"""Place an already accepted witness inside the unreduced magnetic family."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from experiments.marginal_symbolic import decode
from research.transfer_solver_20260915.budget import OUT, dump


def in_span(target, basis):
    pivots = {}
    for polynomial in basis:
        row = polynomial.copy()
        while row:
            first = min(row)
            if first not in pivots:
                q = row[first]
                pivots[first] = {w: c/q for w, c in row.items()}
                break
            value = row[first]
            for w, c in pivots[first].items():
                row[w] = row.get(w, F(0))-value*c
                if not row[w]:
                    del row[w]
    row = target.copy()
    while row:
        first = min(row)
        if first not in pivots:
            return False
        value = row[first]
        for w, c in pivots[first].items():
            row[w] = row.get(w, F(0))-value*c
            if not row[w]:
                del row[w]
    return True


def run(name):
    source = OUT/'cases'/name
    control = OUT/'controls'/(name+'_magnetic')
    read = lambda p: json.loads(p.read_text())
    original, target = read(source/'prepared/frame.json'), read(control/'prepared/frame.json')
    keys = ('fixture_sha256', 'groups', 'rows', 'number_basis', 'spin_basis', 'ladder_basis')
    if any(original[k] != target[k] for k in keys):
        raise ValueError('Control changed the physical dictionary or ideal spaces')
    cert_path = source/'certificate.json'
    cert = read(cert_path)
    allowed = [g['words'] for g in target['groups']]
    for block in cert['core']['blocks']:
        if block['words'] not in allowed:
            raise ValueError('Accepted factor is outside the unreduced magnetic dictionary')
        if any(len(row) != len(block['words']) for row in block['factor']):
            raise ValueError('Factor coordinates do not match its dictionary')
    m = target['modes']
    multipliers = [('number_multiplier', cert['core']['number_multiplier'], target['number_basis'], 4),
                   ('alpha_multiplier', cert['alpha_multiplier'], target['spin_basis'], 2),
                   ('spin_ladder_multiplier', cert['spin_ladder_multiplier'], target['ladder_basis'], 2)]
    for label, polynomial, basis, degree in multipliers:
        if not in_span(decode(polynomial, m, degree), [decode(p, m, degree) for p in basis]):
            raise ValueError(('Accepted multiplier outside control span', label))
    interval = read(source/'interval.json')
    dump(control/'known_feasible_witness.json', {'source_certificate': str(cert_path),
        'source_certificate_sha256': hashlib.sha256(cert_path.read_bytes()).hexdigest(),
        'source_interval': interval, 'factor_dictionary_membership_exact': True,
        'ideal_multiplier_span_membership_exact': True, 'physical_problem_and_ideals_identical': True,
        'new_energy_replay_performed_here': False, 'energy_replay_dependency': str(source/'lower.json'),
        'consequence': 'This unreduced family contains the already independently accepted source bound. A worse finite control result cannot be attributed to insufficient expressiveness of this family.'})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', choices=('h4_control', 'h8_cold'))
    run(p.parse_args().case)
