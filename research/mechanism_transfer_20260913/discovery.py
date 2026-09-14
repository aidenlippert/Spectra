"""Spin completion reuses pricing, selection, and CAR contraction machinery.

Construction and solve adapt the preceding frozen Model to the new exact
frame and acceptance schema. The preceding module and artifacts are unchanged.
"""
from fractions import Fraction as F
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
from research.molecular_collective_20260913.core import digest, extract, factor_operators, retained_polynomial, tail_replay
from research.spin_subspace_20260913.discovery import Model as PreviousModel, DIRECTION_DENOMINATOR
from research.mechanism_transfer_20260913.core import frames, replay


class Model(PreviousModel):
    def __init__(self, data, tail):
        start = time.monotonic(); self.data = data; self.tail = tail
        tail_replay(data, tail); p = extract(data, tail['center_number'])
        if not 4 <= p['modes'] <= 16 or p['particles'] != p['spatial'] or len(tail['factors']) != 2*p['spatial']-2:
            raise ValueError('Transfer requires a half-filled model with 2*s-2 retained patterns, up to sixteen spin orbitals')
        self.p = p; self.h = retained_polynomial(p, tail)
        base, signature, symmetry = partition(self.h, p['modes'], 'quadratic', True)
        self.symmetry = symmetry
        self.base = [{'name': g['name'], 'kind': 'square', 'polynomials': [canonical(mono(w)) for w in g['words']]} for g in base]
        patterns = [q for _, q in factor_operators(p, tail)]
        self.frames = frames(patterns, p['modes'], signature)
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


    def solve(self, span, out, remaining, separate=False, solver_cap=40., whiten=True):
        start = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
        anti, frame_ids = self.groups(span); groups = self.base+anti
        maps = list(self.base_maps); transforms = [None]*len(self.base); conditioning = []
        for g, gid in zip(anti, frame_ids):
            W, stats = whitening_transform(g['polynomials']) if whiten and not separate else (None, None)
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
        built = time.monotonic(); solver_seconds = min(solver_cap, remaining-(built-start)-25.)
        if solver_seconds < 1:
            raise RuntimeError('No solve/export/replay reserve remains')
        metadata = {'directions': len(span), 'separate': separate, 'whiten': whiten, 'anti_dimensions': [len(g['polynomials']) for g in anti],
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
        cert['kind'] = 'mechanism_spin_subspace_v1'
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
