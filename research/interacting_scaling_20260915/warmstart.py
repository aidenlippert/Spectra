"""Transport a retained local-family proposal into a nested larger family."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy import linalg
from research.interacting_scaling_20260915.budget import dump


def basis_embedding(old_words, old_basis, new_words, new_basis):
    positions = {tuple(map(tuple, w)): j for j, w in enumerate(new_words)}
    embedded = np.zeros((len(new_words), old_basis.shape[1]))
    for i, word in enumerate(old_words):
        if tuple(map(tuple, word)) not in positions:
            raise ValueError('The new dictionary omits a retained operator')
        embedded[positions[tuple(map(tuple, word))]] = old_basis[i]
    if embedded.shape == new_basis.shape and np.array_equal(embedded, new_basis):
        return np.eye(new_basis.shape[1])
    A = linalg.lstsq(new_basis, embedded, lapack_driver='gelsy')[0]
    old_integer = np.rint(embedded*6).astype(np.int64)
    new_integer = np.rint(new_basis*6).astype(np.int64)
    if np.max(abs(old_integer/6-embedded)) > 1e-14 or np.max(abs(new_integer/6-new_basis)) > 1e-14:
        raise ValueError('Only exact denominator-six bases are supported')
    for denominator in (6, 12, 36, 72):
        integer = np.rint(A*denominator).astype(np.int64)
        overflow_bound = new_integer.shape[1]*int(np.max(abs(new_integer)))*int(np.max(abs(integer)))
        if overflow_bound >= 2**63: raise ValueError('Exact integer embedding check would overflow')
        if np.array_equal(new_integer@integer, denominator*old_integer):
            return integer/denominator
    raise ValueError('No exact checked basis embedding was recovered')


def transport(source, proposal, target):
    started = time.monotonic()
    old = json.loads((source/'prepared/frame.json').read_text())
    new = json.loads((target/'prepared/frame.json').read_text())
    if old['fixture_sha256'] != new['fixture_sha256']:
        raise ValueError('Nested-family continuation requires the exact same Hamiltonian')
    old_bases = np.load(source/'prepared/bases.npz')
    new_bases = np.load(target/'prepared/bases.npz')
    raw = np.load(source/proposal/'export/raw.npz')
    def name(meta, block): return meta['groups'][block['physical_group']]['name']
    indices = {name(new, b): j for j, b in enumerate(new['blocks'])}
    grams = [np.zeros((b['dimension'], b['dimension'])) for b in new['blocks']]
    for k, block in enumerate(old['blocks']):
        if name(old, block) not in indices: raise ValueError('A retained Gram block is missing')
        j = indices[name(old, block)]
        A = basis_embedding(old['groups'][block['physical_group']]['words'], old_bases[f'V_{k}'],
            new['groups'][new['blocks'][j]['physical_group']]['words'], new_bases[f'V_{j}'])
        grams[j] += A@raw[f'Q_{k}']@A.T
    x = np.zeros(2+len(new['number_basis'])+len(new['spin_basis'])+len(new['ladder_basis']))
    x[0] = raw['x'][0]
    old_offset = new_offset = 1
    for field in ('number_basis', 'spin_basis'):
        lookup = {json.dumps(p, sort_keys=True): j for j, p in enumerate(new[field])}
        for k, p in enumerate(old[field]):
            key = json.dumps(p, sort_keys=True)
            if key not in lookup: raise ValueError('A retained ideal multiplier is missing')
            x[new_offset+lookup[key]] += raw['x'][old_offset+k]
        old_offset += len(old[field]); new_offset += len(new[field])
    x[new_offset] = raw['x'][old_offset]
    old_offset += 1; new_offset += 1
    if old['ladder_basis'] != new['ladder_basis']: raise ValueError('Spin ladder domain changed')
    x[new_offset:] = raw['x'][old_offset:]
    destination = target/'nested_start.npz'
    if destination.exists(): raise ValueError('Preserve the earlier transported proposal')
    np.savez_compressed(destination, x=x, y=-np.load(target/'prepared/physical_dual.npy'),
        **{f'Q_{k}': Q for k, Q in enumerate(grams)})
    dump(target/'nested_start.json', {'source_case': str(source), 'source_proposal': proposal,
        'exact_operator_embedding_checked': True, 'global_full_teacher_used': False,
        'source_discovery_cost_additional': True, 'numerical_proposal_requires_complete_exact_replay': True,
        'dual_initializer': 'Target MPS moments', 'seconds': time.monotonic()-started})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('source', type=Path); p.add_argument('proposal'); p.add_argument('target', type=Path)
    a = p.parse_args(); transport(a.source.resolve(), a.proposal, a.target.resolve())
