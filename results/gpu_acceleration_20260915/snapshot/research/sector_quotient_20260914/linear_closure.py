"""Restore the exact number-dressed linear operators lost under compression.

Nhat*a_i=(N-1)*a_i on N particles, and its adjoint a_i†*Nhat=N*a_i†.
Dropping mixed-block linear words in a full cubic dictionary uses these
identities. A subsequently truncated cubic span need not contain the required
directions. Append those exact directions, without selecting another rank.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
import time
import numpy as np
from experiments.marginal_symbolic import product, number_shift, mono, canonical, add, scale
from research.sector_quotient_20260914.budget import OUT
from research.reconstruction_compression_20260914.inputs import dump, sha


def vectors(words, modes):
    lookup = {}
    for j, word in enumerate(words):
        p = canonical(mono(tuple(tuple(x) for x in word)))
        if len(p) != 1: raise ValueError('Expected a monomial cubic dictionary')
        w, sign = next(iter(p.items())); lookup[w] = (j, sign)
    columns = []; labels = []
    for i in range(modes):
        p = product(number_shift(modes, 0), mono(((0, i),)))
        if not p or not all(w in lookup for w in p): continue
        column = [F(0)]*len(words)
        for w, value in p.items():
            j, sign = lookup[w]; column[j] = value/sign/4
        rebuilt = canonical(add(*(scale(mono(tuple(tuple(x) for x in w)), c) for w, c in zip(words, column) if c)))
        if rebuilt != scale(p, F(1, 4)): raise AssertionError('Number-dressed column identity failed')
        columns.append(column); labels.append(i)
    return columns, labels


def run():
    if (OUT/'manifest.json').exists(): raise RuntimeError('Sealed campaign')
    started = time.monotonic(); source = OUT/'prepared'; destination = OUT/'prepared_linear_closure'
    destination.mkdir(exist_ok=False)
    meta = json.loads((source/'frame.json').read_text())
    data = json.loads(Path(meta['fixture']).read_text()); seed = np.load(source/'seed.npz')
    old_raw = np.load(meta['seed']); new = {key: seed[key] for key in seed.files}; details = []
    for k, pair in enumerate(meta['pairs']):
        i = pair['members'][0]; group = meta['groups'][i]
        if not group['name'].startswith('mixed'): continue
        columns, labels = vectors(group['words'], data['modes'])
        if not columns: continue
        D = np.array(columns, dtype=float).T; V = seed[f'V_{k}']; q = seed[f'Q_{k}']
        distances = np.linalg.norm(D-V@np.linalg.lstsq(V, D, rcond=None)[0], axis=0)/np.linalg.norm(D, axis=0)
        new[f'V_{k}'] = np.column_stack((V, D))
        Q = np.zeros((q.shape[0]+len(labels), q.shape[1]+len(labels))); Q[:q.shape[0], :q.shape[1]] = q
        new[f'Q_{k}'] = Q
        details.append({'pair': k, 'orbitals': labels, 'previous_rank': V.shape[1], 'new_rank': V.shape[1]+len(labels),
                        'relative_distances_from_old_span': distances.tolist(), 'new_columns_are_exact_quarter_integer': True})
    if sum(len(d['orbitals']) for d in details) != 16: raise AssertionError('Unexpected H8 linear completion size')
    # Preserve the original quartic guide dual, solely as a numerical initializer.
    new['dual'] = old_raw['dual']; new['x'] = old_raw['x']
    np.savez_compressed(destination/'seed.npz', **new)
    for path in source.iterdir():
        if path.name not in ('frame.json', 'seed.npz'): shutil.copyfile(path, destination/path.name)
    meta.update(original_seed=meta['seed'], original_seed_sha256=meta['seed_sha256'],
                seed=str(destination/'seed.npz'), seed_sha256=sha(destination/'seed.npz'),
                exact_number_dressed_linear_completion=details, source_prepared_frame_sha256=sha(source/'frame.json'))
    dump(destination/'frame.json', meta)
    rec = {'appended_left_directions': 16, 'independent_partner_directions': 16, 'details': details,
           'full_Gram_teacher_used': False, 'higher_MPS_moments_used': False, 'many_body_states_enumerated': 0,
           'same_old_spans_included_exactly': True, 'seconds': time.monotonic()-started}
    dump(destination/'closure.json', rec); print(json.dumps(rec), flush=True)


if __name__ == '__main__': run()
