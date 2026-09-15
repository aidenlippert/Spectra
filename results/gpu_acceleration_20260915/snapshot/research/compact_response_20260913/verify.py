"""Read-only, standard-library replay of the compact-response experiment."""
import hashlib
import json
import sys
import time

from research.compact_response_20260913 import program, closure, trial


def replay():
    start = time.monotonic(); cases = []
    for name in ('h6', 'h8', 'fresh_h6_1p6'):
        data, tail, reference = trial.load_case(name) if name == 'fresh_h6_1p6' else program.load_case(name)
        cert = json.loads((program.OUT/f'{name}_program.json').read_text())
        cases.append({'case': name, 'receipt': program.check(data, tail, cert)})
    data, tail, reference = program.load_case('h6')
    response = json.loads((program.OUT/'h6_program.json').read_text())
    cert = json.loads((program.OUT/'expanded_closure_certificate.json').read_text())
    complete = closure.check(data, tail, response, cert, reference)
    groups, hden, scalar_cost = closure.blocks(data)
    scalar = closure.check_scalar_obstruction(groups, hden, response, json.loads((program.OUT/'scalar_obstruction.json').read_text()))
    trials = []
    for name in ('h6', 'fresh_h6_1p6'):
        data, tail, reference = trial.load_case(name)
        response = json.loads((program.OUT/f'{name}_program.json').read_text())
        cert = json.loads((program.OUT/f'{name}_trial.json').read_text())
        trials.append({'case': name, 'receipt': trial.check(data, tail, reference, response, cert)})
    forbidden = [name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules]
    if forbidden:
        raise AssertionError('Numerical package imported during exact replay')
    return {'programs': cases, 'complete_H6': complete, 'scalar_obstruction': scalar,
        'scalar_obstruction_additional_cost': scalar_cost, 'trials': trials,
        'numerical_packages_loaded': forbidden, 'wall_seconds': time.monotonic()-start}


def check_manifest():
    manifest = program.OUT/'manifest.json'
    entries = json.loads(manifest.read_text())['files']; mismatches = []
    for name, expected in entries.items():
        raw = (program.ROOT/name).read_bytes()
        if len(raw) != expected['bytes'] or hashlib.sha256(raw).hexdigest() != expected['sha256']:
            mismatches.append(name)
    if mismatches:
        raise ValueError('Experiment manifest mismatch: '+', '.join(mismatches))
    return {'files_verified': len(entries), 'manifest_sha256': hashlib.sha256(manifest.read_bytes()).hexdigest()}


if __name__ == '__main__':
    integrity = check_manifest()
    result = replay()
    print(json.dumps({'integrity': integrity, 'full_H6_width_mHa': result['complete_H6']['width_mHa'],
        'trial_upper_gains_mHa': {r['case']: r['receipt']['accepted_gain_mHa'] for r in result['trials']},
        'numerical_packages_loaded': result['numerical_packages_loaded'], 'replay_seconds': result['wall_seconds']}, indent=2))
