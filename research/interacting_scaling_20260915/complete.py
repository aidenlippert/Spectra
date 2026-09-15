"""One complete exact replay, including transfer to the original Hamiltonian."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time
from research.interacting_scaling_20260915.budget import dump
from research.molecular_collective_20260913.core import digest
from research.transfer_followup_20260915.rotated_upper import rotate_and_round
from research.nvidia_followup_20260915.strict_replay import run as exact_replay


def run(original, case, proposal, output):
    start = time.monotonic()
    data = json.loads((original/'fixture.json').read_text())
    rotation = json.loads((case/'rotation.json').read_text())
    if rotation['original_fixture_sha256'] != digest(data):
        raise ValueError('Rotation is bound to another original model')
    recomputed, allowance = rotate_and_round(data, rotation['integer_matrix'], int(rotation['denominator']))
    supplied = json.loads((case/'fixture.json').read_text())
    if supplied != recomputed: raise ValueError('Exact rotated coefficients did not reproduce')
    # This invocation executes every upper/lower check itself. It does not
    # accept a previously supplied success flag or an unverified interval.
    exact_replay(case, proposal, output)
    interval = json.loads((output/'interval.json').read_text())
    if json.loads((output/'fixture.json').read_text()) != supplied:
        raise ValueError('Replayed input changed during original-model transfer')
    L = F(interval['lower_Ha'])-allowance
    U = F(interval['upper_Ha'])+allowance
    if U < L: raise ValueError('Inconsistent exact endpoints')
    witnesses = [output/'fixture.json', output/'certificate.json', output/'nonsinglet.json', case/'mps/state.json']
    dump(case/'original_interval.json', {'lower_Ha': str(L), 'upper_Ha': str(U), 'width_Ha': str(U-L),
        'width_mHa': float(1000*(U-L)), 'target_met': U-L <= F(1, 625),
        'rotation_allowance_each_endpoint_Ha': str(allowance),
        'original_fixture_sha256': digest(data), 'rotated_fixture_sha256': digest(supplied),
        'upper_and_lower_witnesses_rechecked': True, 'fresh_complete_replay_in_this_invocation': True,
        'witness_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in witnesses},
        'full_fixed_N_determinants_enumerated': 0, 'complete_seconds': time.monotonic()-start,
        'argument': 'Exact paired orbital unitary; coefficient allowance charged separately at both endpoints.'})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('original', type=Path); p.add_argument('case', type=Path)
    p.add_argument('proposal'); p.add_argument('output', type=Path); a = p.parse_args()
    run(a.original.resolve(), a.case.resolve(), a.proposal, a.output.resolve())
