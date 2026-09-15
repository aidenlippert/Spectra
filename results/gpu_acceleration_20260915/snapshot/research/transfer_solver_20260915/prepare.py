"""Construct full coefficient equations and sparse highest-weight maps afresh."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsqr, splu
from experiments.marginal_coefficient import gram_map
from experiments.marginal_symbolic import encode, mono, product, number_shift, adj, canonical
from research.collective_completion_20260914.spin_replay import setup, alpha_shift
from research.collective_completion_20260914.spin_screen import spin_squared, ladder_ideal
from research.collective_completion_20260914.spin_rows import build as spin_rows
from research.certificate_scaling.direct_sparse_discovery import sparse_columns
from research.reconstruction_compression_20260914.moments import Proposal
from research.reconstruction_compression_20260914.inputs import sha
from research.transfer_solver_20260915.algebra import frame, highest_weights
from research.transfer_solver_20260915.budget import dump


def build(case):
    started = time.monotonic()
    folder = case/'prepared'
    folder.mkdir(exist_ok=False)
    data = json.loads((case/'fixture.json').read_text())
    state = json.loads((case/'mps/state.json').read_text())
    h, hs, delta = setup(data)
    m, n = data['modes'], data['particles']
    groups, rows, number_basis, stats = frame(hs, m)
    specs = highest_weights(groups)
    # This is an explicit computational cap, not a claim of an inadequate cone.
    forecast = {'rows': len(rows), 'full_word_pairs': sum(len(g['words'])**2 for g in groups),
                'retained_Gram_entries': sum(V.shape[1]**2 for _, V, _ in specs),
                'largest_Gram_block': max(V.shape[1] for _, V, _ in specs),
                'number_ideal_columns': len(number_basis), 'stats': stats}
    dump(folder/'forecast.json', forecast)
    print(json.dumps({'stage': 'forecast', **forecast}), flush=True)
    if len(rows) > 180000 or forecast['retained_Gram_entries'] > 5000000:
        raise RuntimeError('Declared preparation size budget exceeded before map allocation')
    lookup = {w: i for i, w in enumerate(rows)}
    T, selected, projection = spin_rows(rows)
    Ts = T[selected, :].tocsc()
    sparse.save_npz(folder/'twirl.npz', T)
    np.save(folder/'selected.npy', selected)
    print(json.dumps({'stage': 'spin_projection', 'seconds': time.monotonic()-started, **projection}), flush=True)
    spin_basis = [p for p in number_basis if max(map(len, p), default=0) <= 2]
    polynomials = [mono(())]+[product(number_shift(m, n), p) for p in number_basis]
    polynomials += [product(alpha_shift(m, n), p) for p in spin_basis]+[spin_squared(m)]
    ladder_basis = []
    for i in range(1, m, 2):
        for j in range(0, m, 2):
            p = mono(((1, i), (0, j)))
            q = ladder_ideal(m, p)
            if all(w in lookup or all(v in lookup for v in canonical(adj(mono(w)))) for w in q):
                ladder_basis.append(p)
                polynomials.append(q)
    free = Ts@sparse_columns([{w: c for w, c in p.items() if w in lookup} for p in polynomials], lookup)
    rhs = Ts@np.array([float(hs.get(w, 0)) for w in rows])
    maps, meta, bases = [], [], {}
    physical_pairs = 0
    for k, (g, V, details) in enumerate(specs):
        M = (Ts@gram_map(groups[g]['words'], lookup)).tocsc()
        K = sparse.csc_matrix(V)
        S = (M@sparse.kron(K, K, format='csc')).tocsc()
        r = V.shape[1]
        transposed = (np.arange(r)[None, :]*r+np.arange(r)[:, None]).ravel()
        S = ((S+S[:, transposed])*.5).tocsc()
        S.eliminate_zeros()
        maps.append(S)
        bases[f'V_{k}'] = V
        meta.append({'physical_group': g, 'dimension': r, **details})
        physical_pairs += len(groups[g]['words'])**2
        if k % 8 == 0:
            print(json.dumps({'stage': 'maps', 'block': k, 'seconds': time.monotonic()-started}), flush=True)
    norms = np.asarray(free.power(2).sum(axis=1)).ravel()
    for S in maps:
        norms += np.asarray(S.power(2).sum(axis=1)).ravel()
    row_scale = 1/np.sqrt(np.maximum(1, norms))
    D = sparse.diags(row_scale)
    free = (D@free).tocsc()
    rhs = row_scale*rhs
    G = (free@free.T).tocsc()
    for k, S in enumerate(maps):
        S = (D@S).tocsc()
        maps[k] = S
        sparse.save_npz(folder/f'map_{k}.npz', S)
        G += S@S.T
    G = ((G+G.T)*.5).tocsc()
    G.eliminate_zeros()
    sparse.save_npz(folder/'normal.npz', G)
    sparse.save_npz(folder/'free.npz', free)
    np.save(folder/'rhs.npy', rhs)
    np.save(folder/'scale.npy', row_scale)
    np.savez_compressed(folder/'bases.npz', **bases)
    map_end = time.monotonic()
    oracle = Proposal(data, state)
    oracle.fill(rows)
    weights = np.array([1 if tuple(i for c, i in w if c) == tuple(i for c, i in w if not c) else 2 for w in rows])
    physical = weights*np.array([oracle.cache[w] for w in rows])
    pullback = lsqr(T[selected].T, T.T@physical, atol=1e-13, btol=1e-13, iter_lim=1000)
    np.save(folder/'physical_dual.npy', pullback[0]/row_scale)
    np.save(folder/'weights.npy', weights)
    np.save(folder/'moments.npy', physical)
    rec = {'kind': 'input_derived_highest_weight_preparation_v1', 'fixture': str(case/'fixture.json'),
           'fixture_sha256': sha(case/'fixture.json'), 'state': str(case/'mps/state.json'),
           'state_sha256': sha(case/'mps/state.json'), 'modes': m, 'particles': n,
           'groups': groups, 'rows': rows, 'blocks': meta,
           'number_basis': [encode(p) for p in number_basis], 'spin_basis': [encode(p) for p in spin_basis],
           'ladder_basis': [encode(p) for p in ladder_basis], 'spin_projection': projection,
           'upper_Ha': json.loads((case/'upper.json').read_text())['upper_Ha'],
           'original_H_spin_defect_Ha': str(delta), 'forecast': forecast,
           'map_build_seconds': map_end-started, 'physical_word_pairs_constructed': physical_pairs,
           'coefficient_nonzeros': sum(S.nnz for S in maps), 'normal_nonzeros': G.nnz,
           'moment_seconds': time.monotonic()-map_end, 'moment_count': len(rows),
           'moment_transfer_steps': oracle.steps, 'moment_pullback_residual': pullback[3],
           'initial_primal': 'zero Gram matrices and zero ideal coordinates; no previous solution',
           'teacher_certificate_used': False, 'old_checkpoint_used': False,
           'many_body_determinants_enumerated': 0, 'total_seconds': time.monotonic()-started}
    dump(folder/'frame.json', rec)
    print(json.dumps({k: v for k, v in rec.items() if k not in ('groups', 'rows', 'blocks', 'number_basis', 'spin_basis', 'ladder_basis')}, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    build(p.parse_args().case.resolve())
