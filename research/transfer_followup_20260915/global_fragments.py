"""Disjoint H4 composition allowing charge redistribution at fixed total N.

Small enumerated local charge certificates are explicit comparison inputs.
Their shared chemical-potential line avoids enumerating charge assignments or
many-body determinants of the composite system.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import time
from experiments.marginal_symbolic import decode
from research.molecular_collective_20260913.core import digest
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, NUM, PREFIX

BASE = OUT/'fragments_global_charge'
SOURCE = OUT/'cases/h4_control'
read = lambda p: json.loads(p.read_text())


def check_line(lowers, neutral, intercept, slope, modes):
    if set(lowers) != set(range(modes+1)):
        raise ValueError('Every local particle sector is required')
    if any(intercept+slope*(k-neutral) > lower for k, lower in lowers.items()):
        raise ValueError('Chemical-potential line exceeds a certified local lower')


def disjoint(data, fragment, copies):
    m, n = fragment['modes'], fragment['particles']
    if (data['modes'], data['particles']) != (copies*m, copies*n):
        raise ValueError('Composite dimensions do not match')
    wanted = {}
    for k in range(copies):
        for w, c in decode(fragment['hamiltonian'], m, 4).items():
            if sum(2*a-1 for a, i in w):
                raise ValueError('Local number conservation required')
            shifted = tuple((a, i+k*m) for a, i in w)
            wanted[shifted] = wanted.get(shifted, F(0))+c
    if decode(data['hamiltonian'], copies*m, 4) != wanted:
        raise ValueError('Composite is not the exact disjoint operator sum')


def construct():
    from research.transfer_followup_20260915.enumerated_baseline import construct as propose
    fragment = read(SOURCE/'fixture.json')
    if (fragment['modes'], fragment['particles']) != (8, 4):
        raise ValueError('This control is capped at the supplied H4 fragment')
    for charge in range(9):
        if charge == 4:
            continue
        case = BASE/f'charge_{charge}'
        dump(case/'fixture.json', {**fragment, 'particles': charge,
            'kind': 'same_rational_H4_operator_in_another_particle_sector',
            'parent_fixture_sha256': digest(fragment)})
        propose(case)


def verify():
    import sys
    from research.transfer_followup_20260915.enumerated_baseline import verify as replay
    from research.transfer_followup_20260915.fragments import construct as product
    from research.correlated_pair_20260913.mps_exact import check as upper
    from research.collective_completion_20260914.spin_screen import check as lower
    start = time.monotonic()
    fragment, state = read(SOURCE/'fixture.json'), read(SOURCE/'mps/state.json')
    L = F(lower(fragment, read(SOURCE/'certificate.json'), read(SOURCE/'nonsinglet.json'))['lower'])
    U = F(upper(fragment, state)['upper_Ha'])
    lowers, local_labels = {4: L}, 0
    h = decode(fragment['hamiltonian'], 8, 4)
    for charge in range(9):
        if charge == 4:
            continue
        case = BASE/f'charge_{charge}'
        data = read(case/'fixture.json')
        if data['modes'] != 8 or data['particles'] != charge or decode(data['hamiltonian'], 8, 4) != h:
            raise ValueError('Charge certificate changed the local operator')
        replay(case)
        result = read(case/'enumerated_baseline/replay.json')
        lowers[charge] = F(result['lower_Ha'])
        local_labels += result['cost']['full_sector_labels_enumerated']
    left = max((lowers[k]-L)/(k-4) for k in range(4))
    right = min((lowers[k]-L)/(k-4) for k in range(5, 9))
    if left > right:
        raise ValueError('No supporting line at the neutral certified lower')
    mu = (left+right)/2
    check_line(lowers, 4, L, mu, 8)
    for bad, slope in (({k: v for k, v in lowers.items() if k != 0}, mu), (lowers, right+1)):
        try:
            check_line(bad, 4, L, slope, 8)
        except ValueError:
            pass
        else:
            raise AssertionError('Incomplete or invalid supporting proof was accepted')
    rows = []
    for copies in (2, 4, 8):
        data, combined = product(fragment, state, copies)
        data.pop('fragment_particle_constraints')
        data['kind'] = 'disjoint_H4_fragments_total_N_only'
        data['interpretation'] = 'Complete fixed total-N sector; local charge redistribution allowed'
        combined['fixture_sha256'] = digest(data)
        disjoint(data, fragment, copies)
        rec = upper(data, combined)
        if F(rec['upper_Ha']) != copies*U:
            raise ValueError('Composite upper failed additivity')
        dump(BASE/f'{copies}_fixture.json', data)
        dump(BASE/f'{copies}_state.json', combined)
        row = {'copies': copies, 'total_particles': copies*4,
            'lower_Ha': str(copies*L), 'upper_Ha': str(copies*U),
            'width_mHa': float(1000*copies*(U-L)), 'local_particle_constraints': False,
            'composite_determinants_enumerated': 0, 'charge_assignments_enumerated': 0}
        dump(BASE/f'{copies}_receipt.json', row)
        rows.append(row)
    bad = deepcopy(data)
    bad['hamiltonian'].append({'word': [[1, 0], [0, 8]], 'coefficient': '1/1000'})
    try:
        disjoint(bad, fragment, copies)
    except ValueError:
        pass
    else:
        raise AssertionError('Interfragment coupling was accepted')
    forbidden = [n for n in ('numpy', 'scipy', 'pyscf', 'quimb') if n in sys.modules]
    if forbidden:
        raise ValueError(('Numerical accepting import', forbidden))
    dump(BASE/'summary.json', {'local_lowers_Ha': {str(k): str(v) for k, v in sorted(lowers.items())},
        'chemical_potential_Ha': str(mu), 'admissible_slope_interval_Ha': [str(left), str(right)],
        'neutral_lower_Ha': str(L), 'results': rows, 'all_local_charge_sectors_checked': True,
        'extra_local_determinant_labels_enumerated': local_labels,
        'invalid_slope_missing_charge_and_coupling_refused': True,
        'monolithic_SOS_optimization_tested': False, 'seconds': time.monotonic()-start,
        'argument': 'H_i >= L + mu(N_i-4) on every local charge sector. Sum on total N=4k gives kL without fixing individual charges.',
        'numerical_imports': forbidden})


def run():
    BASE.mkdir(exist_ok=False)
    steps = [('construct', 90, [NUM, '-B', '-m', 'research.transfer_followup_20260915.global_fragments', 'construct']),
             ('verify', 180, [STD, '-B', '-S', '-m', 'research.transfer_followup_20260915.global_fragments', 'verify'])]
    dump(BASE/'protocol.json', {'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'local_enumeration_explicit': True, 'no_monolithic_optimizer_claim': True, 'steps': steps})
    outcomes = []
    for name, seconds, command in steps:
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', 'global_fragments_'+name,
            '--seconds', str(seconds), '--']+command, cwd=ROOT)
        outcomes.append({'stage': name, 'exit_code': result.returncode})
        if result.returncode:
            break
    dump(BASE/'execution.json', outcomes)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'construct', 'verify'))
    a = p.parse_args()
    {'run': run, 'construct': construct, 'verify': verify}[a.action]()
