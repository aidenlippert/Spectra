"""Bind the exact cone obstruction to an independent prior physical lower."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from experiments.marginal_symbolic import decode
from research.certificate_scaling.streaming_reference_upper import upper
from research.certificate_scaling.wedge_spectral_bound import extract as extract_residual, replay_residual
from research.joint_patterns_20260913.dual import check as check_dual
from research.joint_patterns_20260913.spin_diagnostic import check_separator

ROOT = Path(__file__).resolve().parents[2]


def run(results, out):
    start = time.monotonic()
    prior = ROOT/'results/molecular_collective_20260913/campaign/h6'
    data = json.loads((prior/'fixture.json').read_text()); tail = json.loads((prior/'rank_10/tail.json').read_text())
    dual_path = results/'dual/witness.json'; dual = json.loads(dual_path.read_text())
    ceiling = check_dual(data, tail, dual)
    interval_path = ROOT/'results/certificate_scaling/cubic_precision/intervals/final_h6.json'
    source_interval = json.loads(interval_path.read_text())
    physical_path = ROOT/source_interval['certificate']; proof_path = ROOT/source_interval['proof']
    if source_interval['method'] != 'spectral' or hashlib.sha256(physical_path.read_bytes()).hexdigest() != source_interval['certificate_sha256'] or hashlib.sha256(proof_path.read_bytes()).hexdigest() != source_interval['proof_sha256']:
        raise ValueError('Prior precision proof hash or method changed')
    physical = json.loads(physical_path.read_text()); m = data['modes']
    if physical['modes'] != m or physical['particles'] != data['particles'] or decode(physical['hamiltonian'], m, 4) != decode(data['hamiltonian'], m, 4):
        raise ValueError('Prior physical lower is for a different Hamiltonian/sector')
    before = time.monotonic(); proof = json.loads(proof_path.read_text())
    if proof['certificate_sha256'] != source_interval['certificate_sha256']:
        raise ValueError('Residual witness is not bound to the prior physical proof')
    residual, _ = extract_residual(physical)
    physical_lower = replay_residual(physical, residual, proof)
    physical_seconds = time.monotonic()-before
    if F(physical_lower['lower']) != F(source_interval['lower']):
        raise ValueError('Prior precision lower did not reproduce exactly')
    reference_path = ROOT/'results/certificate_scaling/active_space_ladder_references_aligned/h6/upper.json'
    reference = json.loads(reference_path.read_text()); upper_value, upper_stats = upper(data, reference['independent_upper'])
    if upper_value != F(reference['upper']) or F(physical_lower['lower']) > upper_value:
        raise ValueError('Independent interval validation failed')
    cap = F(ceiling['original_lower_ceiling_Ha'])
    error_floor = F(physical_lower['lower'])-cap; width_floor = upper_value-cap
    if error_floor <= 0:
        raise ValueError('No physical accuracy obstruction established')
    separator_path = results/'spin_diagnostic/separator.json'
    separator = check_separator(data, tail, dual, json.loads(separator_path.read_text()))
    forbidden = sorted(name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules)
    if forbidden:
        raise AssertionError('Numerical package imported during obstruction audit')
    result = {'dual': ceiling, 'dual_sha256': hashlib.sha256(dual_path.read_bytes()).hexdigest(),
        'independent_physical_lower': physical_lower,
        'independent_physical_source_path': str(physical_path.relative_to(ROOT)),
        'independent_physical_source_sha256': hashlib.sha256(physical_path.read_bytes()).hexdigest(),
        'independent_physical_source_bytes': physical_path.stat().st_size,
        'independent_residual_proof_path': str(proof_path.relative_to(ROOT)),
        'independent_residual_proof_sha256': hashlib.sha256(proof_path.read_bytes()).hexdigest(),
        'independent_residual_proof_bytes': proof_path.stat().st_size,
        'independent_physical_replay_seconds': physical_seconds,
        'reference_upper': upper_stats, 'reference_source_sha256': hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        'minimum_true_ground_error_Ha': str(error_floor), 'minimum_true_ground_error_mHa': float(1000*error_floor),
        'minimum_width_with_frozen_upper_Ha': str(width_floor), 'minimum_width_with_frozen_upper_mHa': float(1000*width_floor),
        'excludes_1_6_mHa_true_error': error_floor > F(16, 10000),
        'spin_separator': separator, 'separator_sha256': hashlib.sha256(separator_path.read_bytes()).hexdigest(),
        'wall_seconds': time.monotonic()-start, 'numerical_packages_loaded': forbidden,
        'scope': 'Exact obstruction for this certificate grammar, using a separately replayed prior physical lower for the identical frozen molecule. Prior proof discovery and exponential upper-witness discovery are not part of this pass lower-bound discovery costs.'}
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); run(args.results, args.out)
