"""Give the retained canonical full construction the same state and optimizations."""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb
from pathlib import Path
import shutil
import time
from research.interacting_scaling_20260915.budget import OUT, dump


def initialize(original, local, name):
    from research.transfer_followup_20260915 import rotated_upper
    rotated_upper.SOURCE, rotated_upper.CASE = original, local
    rotated_upper.verify()
    case = OUT/'cases'/name
    case.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(original/'fixture.json', case/'fixture.json')
    shutil.copyfile(local/'original_upper.json', case/'upper.json')
    (case/'mps').mkdir()
    # This state retains its truthful LOCAL fixture hash. It is only a source
    # of transported moments here; acceptance must use --rotated and check it.
    shutil.copyfile(local/'mps/state.json', case/'mps/state.json')
    dump(case/'design.json', {'reference_full_canonical_family': True, 'magnetization': 0})
    dump(case/'local_guide.json', {'source': str(local), 'canonical_MPS_claimed': False,
        'acceptance_requires_rotated_upper': True, 'source_discovery_cost_included_in_cold_parent': True})


def prepare(case, local):
    import numpy as np
    from research.transfer_followup_20260915.correlated_guide import compound, canonical_moment
    from research.reconstruction_compression_20260914.moments import Proposal
    from research.molecular_collective_20260913.core import digest
    from research.transfer_solver_20260915 import prepare as reference
    from research.collective_completion_20260914 import spin_rows
    from research.interacting_scaling_20260915.spin_patterns import twirl
    from experiments.marginal_symbolic import decode
    started = time.monotonic()
    data = json.loads((case/'fixture.json').read_text())
    rotated = json.loads((local/'fixture.json').read_text())
    rotation = json.loads((local/'rotation.json').read_text())
    if rotation['original_fixture_sha256'] != digest(data): raise ValueError('Guide belongs to a different original input')
    state = json.loads((local/'mps/state.json').read_text())
    oracle = Proposal(rotated, state)
    R = np.array([[float(F(v, int(rotation['denominator']))) for v in row] for row in rotation['integer_matrix']])
    O = np.kron(R, np.eye(2))
    matrices, indices = {}, {}
    for degree in range(1, 4):
        for alpha in range(degree+1):
            sets = [v for v in combinations(range(data['modes']), degree) if sum(i%2 == 0 for i in v) == alpha]
            pairs = {(i, j): tuple((1, p) for p in a)+tuple((0, q) for q in reversed(b))
                for i, a in enumerate(sets) for j, b in enumerate(sets[:i+1])}
            oracle.fill(pairs.values())
            matrix = np.zeros((len(sets), len(sets)))
            for (i, j), w in pairs.items(): matrix[i, j] = matrix[j, i] = oracle.cache[w]
            W = compound(O, sets)
            matrix = W@matrix@W.T
            expected = comb(state['spin_counts'][0], alpha)*comb(state['spin_counts'][1], degree-alpha)
            if abs(float(np.trace(matrix))-expected) > 1e-7*max(1, expected): raise ValueError('Transported charge trace failed')
            key = degree, alpha
            matrices[key], indices[key] = matrix, {v: i for i, v in enumerate(sets)}
    energy = sum(float(c)*canonical_moment(w, matrices, indices) for w, c in decode(data['hamiltonian'], data['modes'], 4).items())
    if abs(energy-float(F(json.loads((case/'upper.json').read_text())['upper_Ha']))) > 1e-6:
        raise ValueError('Transported moments disagree with the checked original upper')
    transport_seconds = time.monotonic()-started
    class Guide:
        def __init__(self, *_): self.cache, self.steps = {}, oracle.steps
        def fill(self, rows):
            self.cache.update({tuple(map(tuple, w)): canonical_moment(tuple(map(tuple, w)), matrices, indices) for w in rows})
    reference.Proposal = Guide
    spin_rows.twirl = twirl
    reference.build(case)
    dump(case/'reference_preparation.json', {'same_exact_spin_pattern_backend': True, 'same_sparse_QR_solver': True,
        'local_state': str(local/'mps/state.json'), 'moment_transport_seconds': transport_seconds,
        'total_seconds': time.monotonic()-started, 'full_cubic_reference_constructed': True,
        'full_fixed_N_determinants_enumerated': 0, 'transported_energy_Ha': energy,
        'basis': 'Retained canonical dictionary, parity masks and complete old number freedoms'})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); sub = p.add_subparsers(dest='action', required=True)
    i = sub.add_parser('initialize'); i.add_argument('original', type=Path); i.add_argument('local', type=Path); i.add_argument('name')
    b = sub.add_parser('prepare'); b.add_argument('case', type=Path); b.add_argument('local', type=Path)
    a = p.parse_args()
    if a.action == 'initialize': initialize(a.original.resolve(), a.local.resolve(), a.name)
    else: prepare(a.case.resolve(), a.local.resolve())
