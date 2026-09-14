"""Fresh new bounds plus a hash-bound comparison to the verified old run."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from research.molecular_collective_20260913.core import digest
from research.spin_completion_20260913.core import replay
from research.spin_enrichment_20260913.replay import run as replay_new

ROOT = Path(__file__).resolve().parents[2]


def best_within(points, count=None, elapsed=None):
    eligible = [p for p in points if (count is None or p['directions'] <= count)
        and (elapsed is None or p['complete_at_seconds'] <= elapsed)]
    return min(eligible, key=lambda p: F(p['width_Ha'])) if eligible else None


def run(results, out):
    start = time.monotonic(); fresh = replay_new(results/'campaign', out)
    old_root = ROOT/'results/spin_enrichment_20260913'
    old_summary = json.loads((old_root/'campaign/continue/summary.json').read_text())
    new_summary = json.loads((results/'campaign/continue/summary.json').read_text())
    if old_summary['source_hashes'] != new_summary['source_hashes'] or old_summary['budget_seconds'] != new_summary['budget_seconds']:
        raise ValueError('Comparison source or budget differs')
    for key in ('fixture_sha256', 'tail_sha256', 'coefficient_rows', 'spin_frame_dimensions', 'base_Gram_dimensions', 'symmetry'):
        if old_summary['model_construction'][key] != new_summary['model_construction'][key]:
            raise ValueError('Comparison changes the Hamiltonian, frame or baseline')
    old_fresh_path = old_root/'fresh_replay.json'; old_fresh = json.loads(old_fresh_path.read_text())
    old_replayed = {r['case']: r for r in old_fresh['rows'] if r['arm'] == 'continue'}
    new_replayed = {r['case']: r for r in fresh['rows'] if r['arm'] == 'continue'}
    U = F(new_summary['upper_Ha']); initial_lower = F(fresh['inherited']['original_lower_Ha'])
    initial = {'case': 'inherited', 'directions': 64, 'complete_at_seconds': 0.,
        'lower_Ha': str(initial_lower), 'width_Ha': str(U-initial_lower), 'width_mHa': float(1000*(U-initial_lower))}
    points = {}
    for name, summary, directory, accepted in (
            ('coefficient', old_summary, old_root/'campaign/continue', old_replayed),
            ('trace', new_summary, results/'campaign/continue', new_replayed)):
        rows = [initial.copy()]
        for row in summary['trials']:
            if not row['accepted'] or row['complete_at_seconds'] > summary['budget_seconds']:
                continue
            record = accepted[row['case']]; path = directory/row['case']/'certificate.json'
            if hashlib.sha256(path.read_bytes()).hexdigest() != record['certificate_sha256'] or F(row['lower_Ha']) != F(record['lower_Ha']) or F(record['upper_Ha']) != U:
                raise ValueError('A comparison point differs from its exact replay')
            rows.append({key: row[key] for key in ('case', 'directions', 'complete_at_seconds', 'lower_Ha')}
                | {'width_Ha': record['width_Ha'], 'width_mHa': record['width_mHa']})
        points[name] = rows
    old_best = best_within(points['coefficient']); old_path = old_root/'campaign/continue'/old_best['case']/'certificate.json'
    prior = ROOT/'results/molecular_collective_20260913/campaign/h6'
    data = json.loads((prior/'fixture.json').read_text()); tail = json.loads((prior/'rank_10/tail.json').read_text())
    control = replay(data, tail, json.loads(old_path.read_text()))
    if F(control['original_lower_Ha']) != F(old_best['lower_Ha']):
        raise ValueError('Frozen control best lower did not reproduce')
    count_table = []
    for count in (64, 80, 96, 112, 128, 144, 160, 176, 192, 214, 256):
        row = {'combination_ceiling': count, **{name: best_within(p, count=count) for name, p in points.items()}}
        row['trace_gain_mHa'] = row['coefficient']['width_mHa']-row['trace']['width_mHa']; count_table.append(row)
    time_table = []
    for elapsed in (60, 120, 180, 240, 300, 360):
        row = {'elapsed_ceiling_seconds': elapsed, **{name: best_within(p, elapsed=elapsed) for name, p in points.items()}}
        row['trace_gain_mHa'] = row['coefficient']['width_mHa']-row['trace']['width_mHa']; time_table.append(row)
    forbidden = [name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules]
    if forbidden:
        raise AssertionError('Numerical import during comparison acceptance')
    fresh.update({'frozen_control_best_fresh': control,
        'frozen_control_replay_sha256': hashlib.sha256(old_fresh_path.read_bytes()).hexdigest(),
        'comparison_points': points, 'by_combination_ceiling': count_table, 'by_elapsed_ceiling': time_table,
        'comparison_wall_seconds': time.monotonic()-start,
        'comparison_scope': 'One historical budget-matched comparison with identical seed, fixed model and proof grammar. Other control checkpoints reuse their previous exact replays with matching hashes; the control best is freshly replayed. No interpolation or randomized benchmark claim.'})
    out.write_text(json.dumps(fresh, indent=2)+'\n'); return fresh


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True); args = parser.parse_args()
    run(args.results, args.out)
