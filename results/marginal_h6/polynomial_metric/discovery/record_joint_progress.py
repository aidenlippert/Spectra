"""Record accepted exact replays; retain previous finite results as history."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from math import comb
import json

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / 'results/marginal_h6/polynomial_metric'


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


checks = []
energies = []
proofs = {}
suffixes = ['joint_metric_proof', 'joint_metric_transfer/proof',
            'vanishing_metric_proof', 'vanishing_metric_transfer/proof']
for suffix in suffixes:
    folder = BASE / suffix
    certificate = read(folder / 'certificate.json')
    receipt = read(folder / 'receipt.json')
    assert receipt == read(folder / 'independent_replay.json')
    checks.append(str((folder / 'independent_replay.json').relative_to(ROOT)))
    item = {'target_lower': certificate['target_lower'], 'receipt': receipt,
            'certificate_bytes': (folder / 'certificate.json').stat().st_size,
            'positive_atom_completions': {}}
    for name in ['weight_proof', 'numerator_proof']:
        counts = Counter()
        for atom in certificate[name]['positive_indicators']:
            count = 1
            for spin in (0, 1):
                required = sum(bool(atom['required'] & (1 << i)) for i in range(spin, 12, 2))
                occupied = sum(bool(atom['occupied'] & (1 << i)) for i in range(spin, 12, 2))
                available, needed = 6-required, 3-occupied
                count *= comb(available, needed) if 0 <= needed <= available else 0
            assert count > 1
            counts[str(count)] += 1
        item['positive_atom_completions'][name] = dict(counts)
    proofs[suffix] = item
    energy = read(folder / 'energy_independent_replay.json')
    assert F(energy['width']) > 0
    energies.append({'certificate': str((folder / 'energy_certificate.json').relative_to(ROOT)),
                     'independent': str((folder / 'energy_independent_replay.json').relative_to(ROOT)),
                     'width': energy['width'], 'width_float': float(F(energy['width'])),
                     'retained_dimension': energy['retained_dimension'],
                     'response_dimension': energy['response_dimension'],
                     'witness_support': energy['witness_support'],
                     'scope': 'Existing witness and response reused and exactly rechecked; no new witness discovery.'})
for original, transfer in [('joint_metric_proof', 'joint_metric_transfer/proof'),
                           ('vanishing_metric_proof', 'vanishing_metric_transfer/proof')]:
    a = read(BASE / original / 'certificate.json')
    b = read(BASE / transfer / 'certificate.json')
    assert a['polynomial_metric'] == b['polynomial_metric']
    assert a['weight_proof'] == b['weight_proof']

obstructions = {}
for suffix in ['quadratic_sos_obstruction', 'localized_quadratic_sos_projected_obstruction']:
    folder = BASE / suffix
    receipt = read(folder / 'receipt.json')
    assert receipt == read(folder / 'independent_replay.json')
    obstructions[suffix] = receipt
    checks.append(str((folder / 'independent_replay.json').relative_to(ROOT)))

remaining = [
    'Replace explicit400-state metric/proof discovery and sector-sized number-ideal completion with controlled construction cost.',
    'Control degree, atom count, coefficient bit size and accuracy as systems grow.',
    'Remove explicit full-sector upper-witness and retained-response work.',
    'Transfer accurate certificates across larger molecular families.',
    'General physical marginal representability and thermodynamic, kinetic and synthesis modeling remain open.'
]
progress = {'status': 'Four finite joint/vanishing gap and energy proofs accepted; general goal remains active.',
            'report': 'research/marginal_joint_metric.md',
            'projector_degree_report': 'research/marginal_projector_degree.md',
            'tests_passed': 350, 'elapsed_seconds': 268.051,
            'gap_and_obstruction_replays': checks, 'energy_replays': energies,
            'new_stdlib_replay_count': len(checks)+len(energies),
            'proofs': proofs, 'fixed_old_metric_obstructions': obstructions,
            'scaling_diagnostic': 'results/marginal_h6/polynomial_metric/projector_scaling_diagnostic.json',
            'discovery_scope': 'Joint LPs and transfer LPs in400 explicit sector coordinates; no general compression claim.',
            'remaining': remaining}
save(BASE / 'joint_progress.json', progress)
old = read(BASE / 'progress.json')
old['latest_progress'] = 'results/marginal_h6/polynomial_metric/joint_progress.json'
old['historical_scope'] = 'This file retains the first fixed-metric result. Its singleton obstruction is resolved for other metrics in latest_progress.'
save(BASE / 'progress.json', old)
validation = read(ROOT / 'results/marginal_final_validation.json')
validation.update(tests_passed=350, elapsed_seconds=268.051,
    latest_validation_scope='350 full marginal tests; ten new standard-library replays: four joint/vanishing gaps, four energies using existing witnesses/responses, two fixed-old-metric full-quadratic SOS obstructions. General scaling remains open.',
    joint_metric_report='research/marginal_joint_metric.md',
    joint_metric_progress='results/marginal_h6/polynomial_metric/joint_progress.json',
    joint_metric_stdlib_replays=checks+[e['independent'] for e in energies])
for energy in energies:
    if energy['certificate'] not in validation['stdlib_replays_passed']:
        validation['stdlib_replays_passed'].append(energy['certificate'])
assert len(validation['stdlib_replays_passed']) == 142
save(ROOT / 'results/marginal_final_validation.json', validation)
print(json.dumps({'new_replays': progress['new_stdlib_replay_count'],
                  'energy_ledger': len(validation['stdlib_replays_passed']),
                  'singleton_positive_atoms': 0, 'tests': 350}))
