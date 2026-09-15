"""Charged sparse response actions and exact variational witness controls."""
from fractions import Fraction as F
import argparse
import json
import math
import sys
import time

from experiments.marginal_symbolic import decode
from research.certificate_scaling.streaming_reference_upper import compile_term, upper
from research.molecular_collective_20260913.core import digest
from research.compact_response_20260913 import program

OUT = program.OUT


def load_case(name):
    if name == 'h6':
        return program.load_case(name)
    if name != 'fresh_h6_1p6':
        raise ValueError('Sparse physical trial is capped at the two H6 controls')
    source = program.ROOT/'results/response_consistency_20260913/fresh_h6_1p6'
    return tuple(json.loads((source/file).read_text()) for file in ('fixture.json', 'tail.json', 'upper.json'))


def combine(*terms):
    result = {}
    for coefficient, vector in terms:
        for state, value in vector.items():
            result[state] = result.get(state, 0.)+coefficient*value
    return {s: x for s, x in result.items() if x != 0.}


class SparseAction:
    """Proposal-only bit action: no stored basis or Hamiltonian action cache."""
    def __init__(self, data):
        if data['modes'] != 12 or data['particles'] != 6:
            raise ValueError('Sparse control is capped at H6')
        self.terms = [t for w, c in decode(data['hamiltonian'], 12, 4).items()
                      if (t := compile_term(w, float(c))) is not None]
        self.calls = self.sources = self.checks = self.maximum_support = 0
        self.distinct_sources = set(); self.destinations = set()

    def __call__(self, vector, project_q=False):
        if len(vector) > 924 or any(type(s) is not int or not 0 <= s < 4096 or s.bit_count() != 6 for s in vector):
            raise ValueError('Invalid bounded six-electron sparse vector')
        result = {}; self.calls += 1; self.sources += len(vector)
        self.checks += len(vector)*len(self.terms); self.distinct_sources.update(vector)
        for state, x in vector.items():
            for required, occupied, flip, parity, c in self.terms:
                if state & required == occupied:
                    target = state ^ flip; self.destinations.add(target)
                    if project_q and target & (3 << 10) != 3 << 10:
                        continue
                    value = c*x*(-1 if (state & parity).bit_count() % 2 else 1)
                    result[target] = result.get(target, 0.)+value
        result = {s: x for s, x in result.items() if x != 0.}
        self.maximum_support = max(self.maximum_support, len(vector), len(result))
        return result

    def receipt(self):
        return {'Hamiltonian_action_calls': self.calls, 'Hamiltonian_action_sources': self.sources,
            'word_state_checks': self.checks, 'distinct_source_states': len(self.distinct_sources),
            'distinct_generated_destination_states_including_projected_out': len(self.destinations),
            'maximum_stored_vector_support': self.maximum_support, 'Hamiltonian_action_cache_states': 0,
            'full_sector_basis_enumerated': False, 'many_body_matrix_entries': 0}


def response_vector(action, W, cert):
    delta = float(F(cert['delta_Ha'])); M = float(F(cert['M_Ha']))
    target = float(F(cert['target_Ha'])); c = (M+delta)/2; a = (M-delta)/2; z = c/a
    previous = {}; current = combine((1/a, W)); tprev = 1.; tcur = z
    for _ in range(1, cert['order']):
        Dq = combine((1., action(current, True)), (-target, current))
        Zq = combine((c/a, current), (-1/a, Dq))
        previous, current = current, combine((2., Zq), (-1., previous), (2*tcur/a, W))
        tprev, tcur = tcur, 2*z*tcur-tprev
    return combine((1/tcur, current))


def propose(data, tail, reference, response):
    import numpy as np
    start = time.monotonic(); program.check(data, tail, response)
    source = reference['independent_upper']; U, source_receipt = upper(data, source)
    if U != F(reference['upper']):
        raise ValueError('Starting upper failed exact replay')
    v = {s: float(a) for s, a in zip(source['states'], source['amplitudes']) if s & (3 << 10) != 3 << 10 and a}
    norm = math.sqrt(sum(x*x for x in v.values()))
    if not norm:
        raise ValueError('No retained component in the supplied trial')
    v = combine((1/norm, v)); action = SparseAction(data)
    W = action(v, True); w = combine((-1., response_vector(action, W, response)))
    wnorm = math.sqrt(sum(x*x for x in w.values()))
    if not wnorm:
        raise ValueError('No nonzero response direction for this trial')
    w = combine((1/wnorm, w)); Hv = action(v); Hw = action(w)
    dot = lambda x, y: sum(a*y.get(s, 0.) for s, a in x.items())
    h2 = np.array([[dot(v, Hv), dot(v, Hw)], [dot(w, Hv), dot(w, Hw)]])
    if abs(h2[0, 1]-h2[1, 0]) > 1e-10:
        raise ValueError('Numerical trial Hamiltonian is not symmetric')
    values, coefficients = np.linalg.eigh(h2)
    vector = combine((coefficients[0, 0], v), (coefficients[1, 0], w))
    integers = {s: int(round(x*10**10)) for s, x in vector.items() if int(round(x*10**10))}
    witness = {'states': sorted(integers), 'amplitudes': [integers[s] for s in sorted(integers)]}
    candidate, receipt = upper(data, witness)
    cert = {'kind': 'response_trial_upper_v1', 'fixture_sha256': digest(data),
        'response_program_sha256': digest(response), 'source_witness_sha256': digest(source),
        'witness': witness, 'upper_Ha': str(candidate)}
    return cert, {'construction_seconds': time.monotonic()-start, 'action_cost': action.receipt(),
        'D_actions': response['order']-1, 'initial_retained_states': len(v),
        'response_direction_states': len(w), 'final_trial_states': len(witness['states']),
        'source_upper_replay': source_receipt, 'candidate_upper_replay': receipt,
        'numerical_2_by_2_minimum_Ha': float(values[0]), 'candidate_upper_Ha': str(candidate),
        'source_upper_Ha': str(U), 'candidate_gain_mHa': float(1000*(U-candidate)),
        'incumbent_retained_if_candidate_worse': candidate >= U}


def check(data, tail, reference, response, cert):
    start = time.monotonic(); response_receipt = program.check(data, tail, response)
    if (cert.get('kind') != 'response_trial_upper_v1' or cert['fixture_sha256'] != digest(data)
            or cert['response_program_sha256'] != digest(response)
            or cert['source_witness_sha256'] != digest(reference['independent_upper'])):
        raise ValueError('Response-trial input binding failed')
    candidate, receipt = upper(data, cert['witness']); U, source_receipt = upper(data, reference['independent_upper'])
    if str(candidate) != cert['upper_Ha'] or U != F(reference['upper']):
        raise ValueError('Response-trial energy did not reproduce')
    return {'candidate_upper_Ha': str(candidate), 'source_upper_Ha': str(U),
        'accepted_upper_Ha': str(min(candidate, U)), 'candidate_gain_mHa': float(1000*(U-candidate)),
        'accepted_gain_mHa': float(1000*max(U-candidate, F(0))), 'response_program_replay': response_receipt,
        'source_upper_replay': source_receipt, 'candidate_upper_replay': receipt,
        'replay_seconds': time.monotonic()-start,
        'scope': 'Exact variational upper from explicit amplitudes. The response provenance is a numerical construction; full retained-remainder positivity is not proved here.'}


def run(replay=False):
    start = time.monotonic(); rows = []
    for name in ('h6', 'fresh_h6_1p6'):
        data, tail, reference = load_case(name); response_path = OUT/f'{name}_program.json'
        program_discovery = None
        if name == 'fresh_h6_1p6' and not replay:
            response, program_discovery = program.propose(data, tail, F(reference['upper'])-F(1, 1000))
            response_path.write_text(json.dumps(response, indent=2)+'\n')
        else:
            response = json.loads(response_path.read_text())
        path = OUT/f'{name}_trial.json'; discovery = None
        if replay:
            cert = json.loads(path.read_text())
        else:
            cert, discovery = propose(data, tail, reference, response)
            path.write_text(json.dumps(cert, indent=2)+'\n')
        receipt = check(data, tail, reference, response, cert)
        rows.append({'case': name, 'trial_certificate_bytes': path.stat().st_size,
            'response_program_bytes': response_path.stat().st_size,
            'source_upper_provenance': 'inherited 200-amplitude FCI witness' if name == 'h6' else 'single Hartree-Fock determinant; no FCI',
            'program_discovery': program_discovery, 'trial_discovery': discovery, 'replay': receipt})
    forbidden = [x for x in ('numpy', 'scipy', 'cvxpy', 'pyscf') if x in sys.modules]
    if replay and forbidden:
        raise AssertionError('Numerical import during exact trial replay')
    result = {'cases': rows, 'wall_seconds': time.monotonic()-start, 'numerical_packages_loaded': forbidden}
    (OUT/('trial_replay.json' if replay else 'trial_discovery.json')).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'cases': [(r['case'], r['replay']['candidate_gain_mHa']) for r in rows],
        'wall_seconds': result['wall_seconds']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--replay', action='store_true')
    run(parser.parse_args().replay)
