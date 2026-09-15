"""Reusable numerical maps for the bounded spin-subspace campaign."""
from fractions import Fraction as F
from pathlib import Path
import json
import time

import cvxpy as cp
import numpy as np
from scipy import sparse

from experiments.marginal_symbolic import canonical, mono, multiplier_basis, number_shift, product
from research.certificate_scaling.adaptive_block_discovery import partition
from research.certificate_scaling.direct_sparse_discovery import sparse_columns
from research.certificate_scaling.polynomial_gram_contraction import prepare, contract, whitening_transform
from research.joint_patterns_20260913.core import prepare_anticommutators
from research.joint_patterns_20260913.discovery import export_compact
from research.joint_patterns_20260913.spin_diagnostic import spin_generator
from research.molecular_collective_20260913.core import digest, extract, factor_operators, retained_polynomial, tail_replay
from research.spin_subspace_20260913.core import combine, replay

DIRECTION_DENOMINATOR = 10**10


class Model:
    def __init__(self, data, tail):
        start = time.monotonic(); self.data = data; self.tail = tail
        tail_replay(data, tail); p = extract(data, tail['center_number'])
        if p['modes'] != 12 or p['particles'] != 6 or len(tail['factors']) != 10:
            raise ValueError('This campaign targets the frozen ten-pattern H6 model')
        self.p = p; self.h = retained_polynomial(p, tail)
        base, signature, symmetry = partition(self.h, p['modes'], 'quadratic', True)
        self.symmetry = symmetry
        self.base = [{'name': g['name'], 'kind': 'square', 'polynomials': [canonical(mono(w)) for w in g['words']]} for g in base]
        patterns = [q for _, q in factor_operators(p, tail)]; parts = {}
        for k, spin in [(-1, -1)]+[(k, spin) for k in range(len(patterns)) for spin in range(2)]:
            for mode in range(p['modes']):
                ref = [k, spin, mode]; q = spin_generator(patterns, ref, p['modes'])
                keys = {signature(w) for w in q}
                if not q or len(keys) != 1:
                    raise ValueError('Empty or symmetry-mixed spin generator')
                key = next(iter(keys)); group = parts.setdefault(key, {'name': 'spin_'+str(key), 'kind': 'anticommutator', 'generators': [], 'polynomials': []})
                group['generators'].append(ref); group['polynomials'].append(q)
        self.frames = [parts[key] for key in sorted(parts)]
        zero = signature(())
        self.basis = [q for q in multiplier_basis(p['modes'], max_body=1) if all(signature(w) == zero for w in q)]
        free_polys = [mono(())]+[product(number_shift(p['modes'], p['particles']), q) for q in self.basis]
        base_prepared, words, base_work = prepare(self.base)
        self.spin_prepared, spin_words, spin_work = prepare_anticommutators(self.frames)
        allwords = words|spin_words|set(self.h)|{w for q in free_polys for w in q}
        if any(len(w) > 4 or signature(w) != zero or sum(2*c-1 for c, _ in w) for w in allwords):
            raise AssertionError('Unexpected coefficient row')
        self.rows = sorted(allwords, key=lambda w: (len(w), w)); self.lookup = {w: i for i, w in enumerate(self.rows)}
        self.free = sparse_columns(free_polys, self.lookup)
        self.rhs = np.array([float(self.h.get(w, 0)) for w in self.rows])
        self.base_maps = [contract(block, g['polynomials'], self.lookup) for g, block in zip(self.base, base_prepared)]
        self.price_maps = []; self.frame_coefficients = []
        for frame, block in zip(self.frames, self.spin_prepared):
            words = block['words']; s = len(words); terms = block['terms']
            self.price_maps.append(sparse.csr_matrix(([float(c) for _, _, c in terms],
                ([index for _, index, _ in terms], [self.lookup[w] for w, _, _ in terms])), shape=(s*s, len(self.rows))))
            self.frame_coefficients.append(np.array([[float(poly.get(w, 0)) for poly in frame['polynomials']] for w in words]))
        self.construction = {'construction_seconds': time.monotonic()-start,
            'coefficient_rows': len(self.rows), 'base_preparation': base_work, 'spin_preparation': spin_work,
            'base_map_nonzeros': sum(matrix.nnz for matrix, _, _ in self.base_maps),
            'pricing_map_nonzeros': sum(matrix.nnz for matrix in self.price_maps),
            'spin_frame_dimensions': [len(g['polynomials']) for g in self.frames],
            'base_Gram_dimensions': [len(g['polynomials']) for g in self.base],
            'fixture_sha256': digest(data), 'tail_sha256': digest(tail), 'symmetry': symmetry,
            'many_body_states_enumerated': 0}

    def full_span(self):
        return [{'group': gid, 'vector': [DIRECTION_DENOMINATOR*int(i == j) for i in range(len(frame['polynomials']))]}
            for gid, frame in enumerate(self.frames) for j in range(len(frame['polynomials']))]

    def groups(self, span):
        groups = []; ids = []
        for gid, frame in enumerate(self.frames):
            directions = [entry['vector'] for entry in span if entry['group'] == gid]
            if not directions:
                continue
            if len(directions) > len(frame['polynomials']) or any(len(v) != len(frame['polynomials']) or any(type(c) is not int for c in v) or not any(v) for v in directions):
                raise ValueError('Invalid proposed subspace')
            groups.append({**frame, 'directions': directions, 'direction_denominator': DIRECTION_DENOMINATOR,
                'polynomials': [combine(frame['polynomials'], v, DIRECTION_DENOMINATOR) for v in directions]})
            ids.append(gid)
        if sum(len(g['directions']) for g in groups) != len(span):
            raise ValueError('Invalid proposed frame index')
        return groups, ids

    def price(self, dual, span):
        start = time.monotonic()
        rows = [tuple(tuple(letter) for letter in w) for w in dual['rows']]
        if rows != self.rows or len(dual['values']) != len(rows):
            raise ValueError('Dual proposal uses a different coefficient map')
        y = np.asarray(dual['values'], dtype=float)
        if not np.all(np.isfinite(y)):
            raise ValueError('Nonfinite pricing moments')
        proposals = []; diagnostics = []
        for gid, (frame, mapping, C) in enumerate(zip(self.frames, self.price_maps, self.frame_coefficients)):
            M = (mapping@y).reshape(C.shape[0], C.shape[0]); gram = C.T@M@C
            values, vectors = np.linalg.eigh((gram+gram.T)/2)
            old = [C@(np.asarray(e['vector'], dtype=float)/DIRECTION_DENOMINATOR) for e in span if e['group'] == gid]
            accepted = []
            for index, value in enumerate(values):
                if value >= -1e-6 or len(accepted) >= 2:
                    break
                z = np.rint(vectors[:, index]*DIRECTION_DENOMINATOR).astype(np.int64)
                column = C@(z.astype(float)/DIRECTION_DENOMINATOR)
                if old:
                    basis = np.column_stack(old)
                    residual = column-basis@np.linalg.lstsq(basis, column, rcond=1e-12)[0]
                    if np.linalg.norm(residual) < 1e-7*max(1., np.linalg.norm(column)):
                        continue
                if np.linalg.norm(column) < 1e-10:
                    continue
                old.append(column); accepted.append({'group': gid, 'vector': list(map(int, z)), 'proposed_moment': float(value)})
            proposals.append(accepted)
            diagnostics.append({'group': gid, 'dimension': len(values), 'minimum_proposed_moment': float(values[0]),
                'negative_eigenvalues': int(np.sum(values < -1e-6)), 'new_independent_candidates': len(accepted)})
        return proposals, {'seconds': time.monotonic()-start, 'frames': diagnostics,
            'scope': 'Floating pricing proposal only; candidate lower bounds require exact replay.'}

    def solve(self, span, out, remaining, separate=False, solver_cap=30.):
        start = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
        anti, frame_ids = self.groups(span); groups = self.base+anti
        maps = list(self.base_maps); transforms = [None]*len(self.base); conditioning = []
        for g, gid in zip(anti, frame_ids):
            W, stats = whitening_transform(g['polynomials']) if not separate else (None, None)
            maps.append(contract(self.spin_prepared[gid], g['polynomials'], self.lookup, transform=W))
            transforms.append(W); conditioning.append(stats)
        x = cp.Variable(self.free.shape[1]); rem = cp.Variable(len(self.rows)); expr = self.free@x
        variables = []; squared = np.asarray(self.free.power(2).sum(axis=1)).ravel(); used_nnz = 0
        for g, (mapping, ix, _) in zip(groups, maps):
            k = len(g['polynomials'])
            if separate and g['kind'] == 'anticommutator':
                positions = [j for j, index in enumerate(ix) if index//k == index % k]
                mapping = mapping[:, positions]; q = cp.Variable(k, nonneg=True); expr += mapping@q
            else:
                q = cp.Variable((k, k), PSD=True); expr += mapping@cp.reshape(q, (k*k,), order='C')[ix]
            variables.append(q); squared += np.asarray(mapping.power(2).sum(axis=1)).ravel(); used_nnz += mapping.nnz
        row_scale = 1/np.maximum(1., np.sqrt(squared))
        equality = cp.multiply(row_scale, expr+rem-self.rhs) == 0
        problem = cp.Problem(cp.Minimize(cp.norm1(rem)-x[0]), [equality, rem[self.lookup[()]] == 0])
        built = time.monotonic(); solver_seconds = min(solver_cap, remaining-(built-start)-20.)
        if solver_seconds < 1:
            raise RuntimeError('No solve/export/replay reserve remains')
        metadata = {'directions': len(span), 'separate': separate, 'anti_dimensions': [len(g['polynomials']) for g in anti],
            'Gram_entries': sum(len(g['polynomials'])**2 for g in self.base)+sum(len(g['polynomials']) if separate else len(g['polynomials'])**2 for g in anti),
            'constructed_map_nonzeros': sum(matrix.nnz for matrix, _, _ in maps), 'used_map_nonzeros': used_nnz,
            'construction_seconds': built-start, 'solver_seconds_limit': solver_seconds,
            'conditioning': conditioning, 'contraction_work': [stats for _, _, stats in maps[len(self.base):]]}
        (out/'construction.json').write_text(json.dumps(metadata, indent=2)+'\n')
        print(json.dumps({'stage': 'solve', 'path': str(out), 'directions': len(span), 'seconds_limit': solver_seconds}), flush=True)
        problem.solve(solver='CLARABEL', tol_gap_abs=1e-8, tol_gap_rel=1e-8, tol_feas=1e-8, max_iter=1000, time_limit=solver_seconds)
        solved = time.monotonic()
        if x.value is None or equality.dual_value is None or any(q.value is None for q in variables):
            raise RuntimeError('No exportable solver proposal')
        grams = [np.diag(q.value) if q.value.ndim == 1 else q.value for q in variables]
        cert, clipped = export_compact(self.data, self.tail, groups, grams, transforms, x.value, self.basis)
        cert['kind'] = 'spin_pattern_subspace_v1'
        for block, g in zip(cert['anti_blocks'], anti):
            block.update({'directions': g['directions'], 'direction_denominator': DIRECTION_DENOMINATOR})
        dual = {'rows': self.rows, 'values': list(map(float, equality.dual_value*row_scale)), 'scope': 'Numerical proposal, not an accepted obstruction.'}
        (out/'dual_proposal.json').write_text(json.dumps(dual, separators=(',', ':'))+'\n')
        (out/'span.json').write_text(json.dumps(span, separators=(',', ':'))+'\n')
        exported = time.monotonic(); exact = replay(self.data, self.tail, cert)
        text = json.dumps(cert, separators=(',', ':'))+'\n'; (out/'certificate.json').write_text(text)
        receipt = {**metadata, 'status': problem.status, 'floating_lower_Ha': -float(problem.value),
            'clipped_negative_Gram_mass': clipped, 'solve_seconds': solved-built, 'export_seconds': exported-solved,
            'accept_seconds': time.monotonic()-exported, 'wall_seconds': time.monotonic()-start,
            'certificate_bytes': len(text.encode()), 'accepted': exact}
        (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
        print(json.dumps({'stage': 'accepted', 'directions': len(span), 'lower_Ha': exact['original_lower_float_Ha'], 'wall_seconds': receipt['wall_seconds']}), flush=True)
        return receipt, dual
