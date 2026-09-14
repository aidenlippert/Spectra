"""Collect frozen campaign receipts without rerunning scientific discovery."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/all_angles_20260913'


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    rows, bindings = [], []
    for name in ('h4', 'h6', 'h8', 'h10'):
        path = OUT / f'intervals/{name}.json'
        interval = read(path)
        provenance = interval['upper_provenance']
        lower, upper = Fraction(interval['lower']), Fraction(interval['upper'])
        if upper - lower != Fraction(interval['width']):
            raise ValueError('Interval arithmetic mismatch')
        if not 0 <= upper - lower <= Fraction(16, 10000):
            raise ValueError('Interval misses target')
        for field in ('reference', 'certificate', 'proof'):
            if provenance[field] is None:
                continue
            artifact = ROOT / provenance[field]
            actual = digest(artifact)
            if actual != interval[field + '_sha256']:
                raise ValueError(f'Hash mismatch: {artifact}')
            bindings.append({'path': provenance[field], 'sha256': actual})
        if digest(ROOT / provenance['fixture']) != provenance['fixture_sha256']:
            raise ValueError('Fixture hash mismatch')
        rows.append({
            'system': name.upper(), 'lower': interval['lower'],
            'upper': interval['upper'], 'width': interval['width'],
            'lower_float': float(lower), 'upper_float': float(upper),
            'width_float': float(upper-lower),
            'passes_0_0016_Ha': True,
            'witness_states': interval['upper_replay']['witness_states'],
            'upper_discovery_through_first_pass_seconds': provenance['upper_discovery_through_pass_seconds'],
            'complete_two_sided_replay_seconds': interval['wall_seconds'],
            'interval_path': str(path.relative_to(ROOT)),
            'interval_sha256': digest(path),
            'fixture_sha256': provenance['fixture_sha256'],
        })
    runs = []
    for path in sorted((OUT / 'selected_refinement').glob('*/receipt.json')):
        receipt = read(path)
        runs.append({'path': str(path.relative_to(ROOT)), 'seconds': receipt['wall_seconds'],
                     'final_basis': receipt['rows'][-1]['basis_size'],
                     'final_upper': receipt['rows'][-1]['upper']})
    coverage = read(OUT / 'coverage.json')
    counts = dict(sorted(Counter(r['status'] for r in coverage['routes']).items()))
    if len(coverage['routes']) != 72 or len({r['route'] for r in coverage['routes']}) != 72:
        raise ValueError('Coverage must account for 72 distinct directions')
    spin = []
    for name in ('h6_restricted', 'h6_full'):
        base = OUT / 'spin_matched' / name
        interval, pre, receipt = read(base/'interval.json'), read(base/'pre_solve.json'), read(base/'receipt.json')
        spin.append({'variant': name, 'lower': interval['lower'], 'lower_float': interval['lower_float'],
                     'width_float': interval['width_float'], 'ideal_variables': pre['ideal_variables'],
                     'coefficient_rows': pre['coefficient_rows'], 'receipt': receipt})
    hashes = read(OUT/'control/download_hashes.json')
    for name, host, expected in [
        ('h8', '146.235.200.232', 'a26379834a79244e3104c2e526f4305084a5e7055b31269df35106332900681d'),
        ('h10', '129.159.32.200', '6ae786c06481b1d4dc8f43aee9b248fa2c0c2680e7787743a45e5c9317f3c48b'),
    ]:
        path = OUT/f'intervals/{name}.json'
        if digest(path) != expected:
            raise ValueError('Remote interval download hash mismatch')
        item = {'path': str(path.relative_to(ROOT)), 'sha256': expected, 'host': host}
        if item not in hashes:
            hashes.append(item)
    for item in hashes:
        if digest(ROOT/item['path']) != item['sha256']:
            raise ValueError(f'Download changed: {item["path"]}')
    (OUT/'control/download_hashes.json').write_text(json.dumps(hashes, indent=2)+'\n')
    summary = {
        'scope': 'First broad campaign: finite-H results, constructed examples, diagnostics, and explicitly unexecuted routes.',
        'subagents': 11, 'existing_lambda_instances_used': 2, 'new_instances_launched': 0,
        'rows': rows, 'coverage_counts': counts, 'coverage_routes': 72,
        'selected_refinement_all_runs': runs,
        'selected_refinement_all_runs_sum_seconds': sum(r['seconds'] for r in runs),
        'timing_scope': 'Run sum includes initial unsuccessful H10 and independent extended rerun. It excludes other agent experiments, previous lower discovery, full interval replay, setup, and queue time; it is not total campaign cost.',
        'spin_matched': spin, 'gpu': read(OUT/'gpu/summary.json'),
        'artifact_bindings': bindings, 'remote_local_hash_matches': len(hashes),
        'focused_tests': '6 passed; control/pytest.log. No whole-repository regression claim.',
        'formal_proof': False, 'fixed_accuracy_scaling_proved': False,
    }
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    coverage_lines = ['# Coverage of the 72 brainstorm directions', '', coverage['scope']+'.', '',
                      'Status describes the evidence actually obtained. Related prerequisites do not count as full implementations.', '',
                      '| Route | Direction | Status | Evidence scope |', '|---:|---|---|---|']
    for route in coverage['routes']:
        coverage_lines.append(f'| {route["route"]} | {route["name"]} | {route["status"]} | {route["evidence_scope"]} |')
    (OUT/'COVERAGE.md').write_text('\n'.join(coverage_lines)+'\n')
    print(json.dumps({'intervals_passed': len(rows), 'remote_local_hash_matches': len(hashes),
                      'coverage_counts': counts, 'selected_refinement_all_runs_seconds': summary['selected_refinement_all_runs_sum_seconds']}, indent=2))


if __name__ == '__main__':
    main()
