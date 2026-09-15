"""Build only declared local/collective products and their exact coefficient rows."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsqr
from experiments.marginal_coefficient import gram_map, dagger
from experiments.marginal_symbolic import encode, mono, product, number_shift, word_product
from research.collective_completion_20260914.spin_replay import setup, alpha_shift
from research.collective_completion_20260914.spin_screen import spin_squared, ladder_ideal
from research.collective_completion_20260914 import spin_rows
from research.interacting_scaling_20260915.spin_patterns import twirl
from research.certificate_scaling.direct_sparse_discovery import sparse_columns
from research.reconstruction_compression_20260914.moments import Proposal
from research.reconstruction_compression_20260914.inputs import sha
from research.interacting_scaling_20260915.dictionary import clustered_frame, ideal_basis, close_rows
from research.interacting_scaling_20260915.budget import dump


def build(case):
    start = time.monotonic()
    folder = case/'prepared'
    folder.mkdir(exist_ok=False)
    data = json.loads((case/'fixture.json').read_text())
    from research.nvidia_followup_20260915.strict_replay import require_supported_sector
    require_supported_sector(data)
    design = json.loads((case/'design.json').read_text())
    magnetic = design.get('magnetization', 0) == 1
    if design.get('magnetization', 0) not in (0, 1):
        raise ValueError('Only the declared singlet and M_S=1 constructions are supported')
    state = None if magnetic else json.loads((case/'mps/state.json').read_text())
    h, hs, delta = setup(data)
    m, n = data['modes'], data['particles']
    groups, specs = clustered_frame(m, design['clusters'], design['collective_pairs'], design.get('complete', False), design.get('collective_clusters', ()))
    if magnetic:
        specs = [(g, np.eye(len(group['words'])), {'kind': 'full_magnetic', 'declared_block': group['name']}) for g, group in enumerate(groups)]
    number_basis = ideal_basis(m, design['clusters'], design.get('complete', False), design.get('number_body', 2))
    spin_basis = [p for p in number_basis if max(map(len, p), default=0) <= 2]
    polynomials = [mono(())]+[product(number_shift(m, n), p) for p in number_basis]
    from experiments.marginal_symbolic import add
    shift = add(alpha_shift(m, n), mono((), -int(magnetic)))
    polynomials += [product(shift, p) for p in spin_basis]+[{} if magnetic else spin_squared(m)]
    ladder_basis = [] if magnetic else [mono(((1, i), (0, j))) for i in range(1, m, 2) for j in range(0, m, 2)]
    polynomials += [ladder_ideal(m, p) for p in ladder_basis]
    forecast = {'retained_Gram_entries': sum(V.shape[1]**2 for g, V, d in specs),
        'largest_Gram_block': max(V.shape[1] for g, V, d in specs),
        'declared_word_pairs': sum(len(groups[g]['words'])**2 for g, V, d in specs),
        'number_ideal_columns': len(number_basis), 'blocks': len(specs),
        'global_cubic_map_constructed': bool(design.get('complete', False)),
        'full_fixed_N_determinants_enumerated': 0, 'design': design}
    dump(folder/'forecast.json', forecast)
    print(json.dumps({'stage': 'forecast', **forecast}), flush=True)
    if forecast['retained_Gram_entries'] > design.get('max_Gram_entries', 2000000):
        raise RuntimeError('Declared Gram allocation envelope exceeded before coefficient maps')
    # This is the actual union of generated products, not a global row table.
    words = set(hs)
    for p in polynomials:
        words.update(p)
    for g, V, details in specs:
        entries = groups[g]['words']
        for left in entries:
            dl = dagger(left)
            for right in entries:
                words.update(w for w, c in word_product(dl, right) if c)
    if magnetic:
        from research.interacting_scaling_20260915.dictionary import representative
        rows = sorted({representative(w) for w in words}, key=lambda w: (len(w), w))
    else:
        rows = close_rows(words)
    if len(rows) > design.get('max_coefficient_rows', 180000):
        raise RuntimeError('Declared coefficient row envelope exceeded before map allocation')
    lookup = {w: i for i, w in enumerate(rows)}
    if magnetic:
        T, selected = sparse.eye(len(rows), format='csr'), np.arange(len(rows))
        projection = {'identity_no_spin_average': True, 'original_rows': len(rows), 'independent_invariant_rows': len(rows)}
    else:
        spin_rows.twirl = twirl
        T, selected, projection = spin_rows.build(rows)
    Ts = T[selected, :].tocsc()
    sparse.save_npz(folder/'twirl.npz', T)
    np.save(folder/'selected.npy', selected)
    print(json.dumps({'stage': 'spin_projection', 'seconds': time.monotonic()-start, **projection}), flush=True)
    free = Ts@sparse_columns([{w: c for w, c in p.items() if w in lookup} for p in polynomials], lookup)
    rhs = Ts@np.array([float(hs.get(w, 0)) for w in rows])
    maps, meta, bases = [], [], {}
    for k, (g, V, details) in enumerate(specs):
        M = (Ts@gram_map(groups[g]['words'], lookup)).tocsc()
        K = sparse.csc_matrix(V)
        S = (M@sparse.kron(K, K, format='csc')).tocsc()
        r = V.shape[1]
        order = (np.arange(r)[None, :]*r+np.arange(r)[:, None]).ravel()
        S = ((S+S[:, order])*.5).tocsc()
        S.eliminate_zeros()
        maps.append(S)
        bases[f'V_{k}'] = V
        meta.append({'physical_group': g, 'dimension': r, **details})
        if k % 8 == 0:
            print(json.dumps({'stage': 'map', 'block': k, 'seconds': time.monotonic()-start}), flush=True)
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
    weights = np.array([1 if tuple(i for c, i in w if c) == tuple(i for c, i in w if not c) else 2 for w in rows])
    if magnetic:
        from math import comb
        from research.interacting_scaling_20260915.singlet_trace import magnetic_trace
        alpha, beta = n//2+1, n//2-1
        dimension = comb(m//2, alpha)*comb(m//2, beta)
        physical = weights*np.array([float(magnetic_trace(mono(w), m, alpha, beta)/dimension) for w in rows])
        transfer_steps = 0
    else:
        oracle = Proposal(data, state)
        oracle.fill(rows)
        physical = weights*np.array([oracle.cache[w] for w in rows])
        transfer_steps = oracle.steps
    pullback = lsqr(T[selected].T, T.T@physical, atol=1e-13, btol=1e-13, iter_lim=1000)
    np.save(folder/'physical_dual.npy', pullback[0]/row_scale)
    np.save(folder/'weights.npy', weights)
    np.save(folder/'moments.npy', physical)
    rec = {'kind': 'direct_local_collective_preparation_v1', 'fixture': str(case/'fixture.json'),
        'fixture_sha256': sha(case/'fixture.json'), 'state': None if magnetic else str(case/'mps/state.json'),
        'state_sha256': None if magnetic else sha(case/'mps/state.json'), 'modes': m, 'particles': n,
        'magnetization': int(magnetic), 'moment_seed': 'Exact magnetic-sector occupation trace' if magnetic else 'MPS',
        'groups': groups, 'rows': rows, 'blocks': meta,
        'number_basis': [encode(p) for p in number_basis], 'spin_basis': [encode(p) for p in spin_basis],
        'ladder_basis': [encode(p) for p in ladder_basis], 'spin_projection': projection,
        'upper_Ha': json.loads((case/'upper.json').read_text())['upper_Ha'],
        'original_H_spin_defect_Ha': str(delta), 'forecast': forecast,
        'map_build_seconds': map_end-start, 'coefficient_rows': len(rows),
        'physical_word_pairs_constructed': forecast['declared_word_pairs'],
        'coefficient_nonzeros': sum(S.nnz for S in maps), 'normal_nonzeros': G.nnz,
        'moment_seconds': time.monotonic()-map_end, 'moment_count': len(rows),
        'moment_transfer_steps': transfer_steps, 'moment_pullback_residual': float(pullback[3]),
        'teacher_certificate_used': False, 'old_checkpoint_used': False,
        'many_body_determinants_enumerated': 0, 'total_seconds': time.monotonic()-start,
        'constructor_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dump(folder/'frame.json', rec)
    print(json.dumps({k: v for k, v in rec.items() if k not in ('groups', 'rows', 'blocks', 'number_basis', 'spin_basis', 'ladder_basis')}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    build(p.parse_args().case.resolve())
