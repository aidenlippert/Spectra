"""Fresh full-degree CAR maps, retaining the old compact operator spans.

Only fixed orbital operator dictionaries and the old compact V maps are read.
The strong full certificate and three-body guiding moments are not inputs.
All degree-six equations survive, including the quartic number multiplier.
"""
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from experiments.marginal_symbolic import encode, mono, product, number_shift, add, adj, canonical
from experiments.marginal_coefficient import gram_map, dagger
from research.reconstruction_compression_20260914.reduced import frame
from research.reconstruction_compression_20260914.inputs import dump, sha
from research.collective_completion_20260914.spin_rows import build as spin_rows
from research.collective_completion_20260914.spin_replay import setup, alpha_shift
from research.collective_completion_20260914.spin_screen import spin_squared, ladder_ideal
from research.certificate_scaling.direct_sparse_discovery import sparse_columns
from research.sector_quotient_20260914.budget import ROOT, OUT


def run():
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Sealed campaign')
    started = time.monotonic()
    folder = OUT/'prepared'
    folder.mkdir(exist_ok=False)
    frozen = json.loads((OUT/'frozen_inputs.json').read_text())['h8']
    fixture = ROOT/frozen['fixture']
    if sha(fixture) != frozen['sha256']['fixture']:
        raise ValueError('Fixture changed')
    data = json.loads(fixture.read_text())
    h, groups, rows, lookup, number_basis, _, stats = frame(data)
    m, n = data['modes'], data['particles']
    _, hs, delta = setup(data)
    T, selected, projection_receipt = spin_rows(rows)
    Ts = T[selected, :].tocsc()
    sparse.save_npz(folder/'twirl.npz', T)
    np.save(folder/'selected.npy', selected)
    print(json.dumps({'stage': 'exact_spin_projection', 'seconds': time.monotonic()-started, **projection_receipt}), flush=True)
    spin_basis = [p for p in number_basis if max(map(len, p), default=0) <= 2]
    polynomials = [mono(())] + [product(number_shift(m, n), p) for p in number_basis]
    polynomials += [product(alpha_shift(m, n), p) for p in spin_basis]
    polynomials += [spin_squared(m)]
    ladder_basis = []
    for i in range(1, m, 2):
        for j in range(0, m, 2):
            p = mono(((1, i), (0, j)))
            q = ladder_ideal(m, p)
            if all(w in lookup or all(v in lookup for v in canonical(adj(mono(w)))) for w in q):
                ladder_basis.append(p)
                polynomials.append(q)
    free = sparse_columns([{w: c for w, c in p.items() if w in lookup} for p in polynomials], lookup)
    sparse.save_npz(folder/'free_full.npz', free)
    sparse.save_npz(folder/'free.npz', (Ts@free).tocsc())
    rhs = np.array([float(hs.get(w, 0)) for w in rows])
    np.save(folder/'rhs.npy', Ts@rhs)
    seed = ROOT/'results/collective_completion_20260914/candidates/h8_spin_r32_invariant_dual48/round_0/raw.npz'
    compact = ROOT/'results/collective_completion_20260914/prepared/h8/frame.json'
    meta = json.loads(compact.read_text())
    if json.loads(json.dumps(groups)) != meta['groups']:
        raise ValueError('Operator dictionaries changed')
    raw = np.load(seed)
    maps = {}
    entries = 0
    nonzeros = 0
    for k, pair in enumerate(meta['pairs']):
        i = pair['members'][0]
        V = raw[f'V_{k}']
        if V.shape[0] != len(groups[i]['words']):
            raise ValueError('Compact map dimension')
        maps[f'V_{k}'] = V
        maps[f'Q_{k}'] = raw[f'Q_{k}']
        for j in pair['members']:
            M = gram_map(groups[j]['words'], lookup).tocsc()
            M = (Ts@M).tocsc()
            M.eliminate_zeros()
            sparse.save_npz(folder/f'map_{j}.npz', M)
            entries += len(groups[j]['words'])**2
            nonzeros += M.nnz
        if k % 8 == 0:
            print(json.dumps({'stage': 'unpaired_operator_maps', 'pair': k, 'seconds': time.monotonic()-started}), flush=True)
    np.savez_compressed(folder/'seed.npz', **maps, old_x=raw['x'])
    result = {'kind': 'unpaired_number_quotient_v1', 'fixture': str(fixture), 'fixture_sha256': sha(fixture),
              'groups': groups, 'rows': rows, 'pairs': meta['pairs'],
              'number_basis': [encode(p) for p in number_basis], 'spin_basis': [encode(p) for p in spin_basis],
              'ladder_basis': [encode(p) for p in ladder_basis],
              'seed': str(seed), 'seed_sha256': sha(seed), 'compact_frame': str(compact), 'compact_frame_sha256': sha(compact),
              'spin_projection': projection_receipt, 'number_basis_dimension': len(number_basis),
              'quartic_number_basis_dimension': sum(max(map(len, p), default=0) == 4 for p in number_basis),
              'many_body_states_enumerated': 0, 'higher_MPS_moments_used': False, 'full_Gram_teacher_used': False,
              'full_operator_word_pairs': entries, 'map_nonzeros': nonzeros, 'seconds': time.monotonic()-started}
    dump(folder/'frame.json', result)
    print(json.dumps({k: v for k, v in result.items() if k not in ('groups', 'rows', 'pairs', 'number_basis', 'spin_basis', 'ladder_basis')}, indent=2))


if __name__ == '__main__':
    run()
