"""Bounded continuation experiments; earlier sealed campaigns stay read-only."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/acceptance_channels_20260915'
SOURCE = ROOT / 'results/interacting_scaling_20260915/cases/h12_main_refined'
CASE = OUT / 'h12'


def dump(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(data, stream, indent=2)
        stream.write('\n')


def receipt(path):
    interval = json.loads((path / 'original_interval.json').read_text())
    lower = json.loads((path / 'exact/lower.json').read_text())
    width = Fraction(interval['upper_Ha']) - Fraction(interval['lower_Ha'])
    if width != Fraction(interval['width_Ha']) or width < 0:
        raise ValueError('Inconsistent stored exact endpoints')
    residual = Fraction(lower['singlet']['residual_l1'])
    return {
        'kind': 'analysis_of_previously_accepted_receipt_not_new_replay',
        'source': str(path), 'width_Ha': str(width),
        'width_mHa': float(1000 * width),
        'singlet_residual_Ha': str(residual),
        'singlet_residual_mHa': float(1000 * residual),
        'residual_share_percent': float(100 * residual / width),
        'hypothetical_width_without_this_allowance_mHa': float(1000 * (width-residual)),
        'target_met': width <= Fraction(1, 625),
    }


def initialize():
    OUT.mkdir(exist_ok=False)
    baseline = receipt(SOURCE)
    dump(OUT / 'baseline.json', baseline)
    CASE.mkdir()
    dependencies = {}
    for name in ('fixture.json', 'upper.json', 'nonsinglet.json', 'rotation.json', 'design.json', 'mps', 'prepared'):
        path = (SOURCE / name).resolve(strict=True)
        (CASE / name).symlink_to(path, target_is_directory=path.is_dir())
        if path.is_file():
            dependencies[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    for name in ('solve/checkpoint.npz', 'solve/export/raw.npz', 'continuation_dependency.json'):
        path = SOURCE / name
        dependencies[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    dump(OUT / 'protocol.json', {
        'declared_UTC': datetime.now(timezone.utc).isoformat(),
        'objective': 'Improve the accepted H12 bound within the same family; audit the supplied compact channel and prototype accepting-cost reduction',
        'fixed_source': str(SOURCE), 'source_sha256': dependencies,
        'original_model': str(ROOT / 'results/interacting_scaling_20260915/models/h12_heldout'),
        'fixed_model_upper_sector_family_and_acceptance': True,
        'initial_numerical_comparisons_seconds_each': 300,
        'initial_fixed_gram_ideal_fit_seconds': 180,
        'exact_replay_cap_seconds': 900,
        'source_preparation_discovery_and_failed_attempts_additional': True,
        'not_a_cold_or_matched_reference_performance_claim': True,
        'source_nonsinglet_proof_preserved': True,
        'new_external_spending': False,
        'github_operations': False,
        'H8_source_dual_bundle': 'Not supplied with the pasted request; its negative value is unverified here',
    })
    print(json.dumps(baseline), flush=True)


def run(name, seconds, command):
    from research.correlated_pair_20260913 import budget
    if not (OUT / 'protocol.json').exists():
        raise ValueError('Declare the experiment before running')
    budget.OUT = OUT / 'runs'
    return budget.run(name, seconds, command)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('initialize', 'receipt', 'run'))
    parser.add_argument('--name')
    parser.add_argument('--seconds', type=int)
    args, command = parser.parse_known_args()
    if command and command[0] == '--':
        command = command[1:]
    if args.action == 'initialize':
        initialize()
    elif args.action == 'receipt':
        data = receipt(SOURCE)
        print(json.dumps(data, indent=2))
        sys.exit(0 if data['target_met'] else 1)
    else:
        sys.exit(run(args.name, args.seconds, command))
