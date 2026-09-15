"""Exact comparison of a completed-spin dual ceiling with the old physical proof."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from experiments.marginal_symbolic import decode
from research.certificate_scaling.cubic_interval_replay import run as physical_interval
from research.spin_enrichment_20260913.full_dual import check

ROOT = Path(__file__).resolve().parents[2]


def run(results, out):
    start = time.monotonic(); prior = ROOT/'results/molecular_collective_20260913/campaign/h6'
    data = json.loads((prior/'fixture.json').read_text()); tail = json.loads((prior/'rank_10/tail.json').read_text())
    dual_path = results/'full_dual/witness.json'
    dual = check(data, tail, json.loads(dual_path.read_text()))
    meta = json.loads((ROOT/'results/certificate_scaling/cubic_precision/intervals/final_h6.json').read_text())
    cert_path = ROOT/meta['certificate']; proof_path = ROOT/meta['proof']; ref_path = ROOT/meta['reference']
    if hashlib.sha256(cert_path.read_bytes()).hexdigest() != meta['certificate_sha256'] or hashlib.sha256(proof_path.read_bytes()).hexdigest() != meta['proof_sha256']:
        raise ValueError('Prior physical proof hashes changed')
    cert = json.loads(cert_path.read_text())
    if cert['modes'] != data['modes'] or cert['particles'] != data['particles'] or decode(cert['hamiltonian'], 12, 4) != decode(data['hamiltonian'], 12, 4):
        raise ValueError('Prior proof concerns a different Hamiltonian or sector')
    comparison_output = out.with_name(out.stem+'_physical_interval.json')
    physical = physical_interval(cert_path, ref_path, comparison_output, method='spectral', proof=proof_path)
    if F(physical['lower']) != F(meta['lower']):
        raise ValueError('Prior exact physical lower did not reproduce')
    ceiling = F(dual['original_lower_ceiling_Ha'])
    gap = F(physical['lower'])-ceiling; width = F(physical['upper'])-ceiling
    forbidden = sorted(name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules)
    if forbidden:
        raise AssertionError('Numerical import during exact obstruction audit')
    result = {'full_completed_spin_dual': dual, 'dual_sha256': hashlib.sha256(dual_path.read_bytes()).hexdigest(),
        'physical_interval': physical, 'minimum_true_ground_error_Ha': str(gap),
        'minimum_true_ground_error_mHa': float(1000*gap), 'excludes_1_6_mHa_true_error': gap > F(16, 10000),
        'minimum_width_with_frozen_upper_Ha': str(width), 'minimum_width_with_frozen_upper_mHa': float(1000*width),
        'wall_seconds': time.monotonic()-start, 'numerical_packages_loaded': forbidden,
        'scope': 'Exact accuracy obstruction only for the complete 492-generator spin-completion cone and its subspaces. The independent physical comparison proof is not a discovery input.'}
    out.write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True); args = parser.parse_args(); run(args.results, args.out)
