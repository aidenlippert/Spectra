"""Exact contradiction test for treating a timed-out conic value as a ceiling."""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('case_directory', type=Path)
    folder = parser.parse_args().case_directory.resolve()
    history_path = folder / 'global_diagonal/subspace_history.json'
    energy_path = folder / 'limit_refined/range_two_replay.json'
    history = json.loads(history_path.read_text())
    energy = json.loads(energy_path.read_text())
    if not energy['accepted']:
        raise ValueError('Accepted exact lower required')
    for name, expected in energy['source_sha256'].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError('Stale energy source: ' + name)
    lower = Fraction(energy['lower_replay']['periodic_lower_density'])
    contradictions = []
    for row in history:
        reported = Fraction(str(row['model_periodic_upper']))
        if reported < lower:
            if row['full_minimum_slack'] >= 0:
                raise ValueError('Unexpected nonnegative numerical slack; investigate')
            contradictions.append({
                'round': row['round'], 'solver_status': row['status'],
                'reported_value': str(reported), 'accepted_periodic_lower': str(lower),
                'exact_lower_minus_reported_value': str(lower - reported),
                'full_minimum_slack': row['full_minimum_slack'],
            })
    if not contradictions:
        raise ValueError('No contradictory reported ceiling found')
    files = {Path(__file__).resolve(), history_path, energy_path}
    result = {
        'accepted': True, 'contradictions': contradictions,
        'source_sha256': {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(files)
        },
        'scope': 'The displayed decimal solver values cannot be global ceilings for this family because an independently accepted exact periodic lower exceeds them. This refutes using these incomplete numerical solves as optimality evidence; it does not refute convex duality or determine the actual optimum.',
    }
    (folder / 'limit_refined/solver_limit_refutation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'accepted': True, 'contradictory_values': len(contradictions)}))


if __name__ == '__main__':
    main()
