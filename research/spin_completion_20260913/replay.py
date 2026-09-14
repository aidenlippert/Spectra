"""Fresh exact interval replay and comparison with the previous cone ceiling."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from research.certificate_scaling.streaming_reference_upper import upper
from research.spin_subspace_20260913.full_dual import check as check_previous_dual
from research.spin_completion_20260913.diagnostic import check_separator
from research.spin_completion_20260913.core import replay

ROOT = Path(__file__).resolve().parents[2]


def run(campaign, out, separate=None):
    start = time.monotonic(); prior = ROOT/'results/molecular_collective_20260913/campaign/h6'
    data = json.loads((prior/'fixture.json').read_text()); tail = json.loads((prior/'rank_10/tail.json').read_text())
    reference_path = ROOT/'results/certificate_scaling/active_space_ladder_references_aligned/h6/upper.json'
    reference = json.loads(reference_path.read_text()); U, upper_receipt = upper(data, reference['independent_upper'])
    if U != F(reference['upper']):
        raise ValueError('Reference upper did not reproduce')
    old_path = ROOT/'results/spin_subspace_20260913/full_dual/witness.json'
    old = check_previous_dual(data, tail, json.loads(old_path.read_text()))
    separators = []
    diagnostic = campaign.parent/'diagnostic'
    for path in sorted(diagnostic.glob('separator_*.json')):
        result = check_separator(data, tail, json.loads(old_path.read_text()), json.loads(path.read_text()))
        separators.append({'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'fresh': result})
    cap = F(old['original_lower_ceiling_Ha']); results = []
    arms = [(name, campaign/name) for name in ('adaptive', 'separate', 'full')]
    if (campaign.parent/'raw_full/summary.json').exists():
        arms.append(('full_raw', campaign.parent/'raw_full'))
    if separate is not None:
        arms.append(('separate_retry', separate))
    for name, directory in arms:
        path = directory/'summary.json'
        if not path.exists():
            continue
        summary = json.loads(path.read_text())
        for row in summary['trials']:
            if not row['accepted'] or row['complete_at_seconds'] > summary['budget_seconds']:
                continue
            cert_path = directory/row['case']/'certificate.json'
            cert = json.loads(cert_path.read_text()); result = replay(data, tail, cert)
            lower = F(result['original_lower_Ha']); width = U-lower
            if lower != F(row['lower_Ha']) or width < 0:
                raise ValueError('Certificate replay mismatch or inconsistent interval')
            record = {'arm': name, 'case': row['case'], 'directions': row['directions'],
                'fresh': result, 'upper_Ha': str(U), 'width_Ha': str(width), 'width_mHa': float(1000*width),
                'meets_1_6_mHa': width <= F(16, 10000),
                'gain_over_previous_cone_ceiling_Ha': str(lower-cap),
                'gain_over_previous_cone_ceiling_mHa': float(1000*(lower-cap)),
                'certificate_sha256': hashlib.sha256(cert_path.read_bytes()).hexdigest()}
            results.append(record)
            print(json.dumps({k: record[k] for k in ('arm', 'case', 'directions', 'width_mHa', 'gain_over_previous_cone_ceiling_mHa')}), flush=True)
    forbidden = sorted(name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules)
    if forbidden:
        raise AssertionError('Numerical package loaded during exact replay')
    result = {'rows': results, 'separators': separators, 'previous_cone': old, 'previous_dual_sha256': hashlib.sha256(old_path.read_bytes()).hexdigest(),
        'reference_upper': upper_receipt, 'reference_sha256': hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        'wall_seconds': time.monotonic()-start, 'numerical_packages_loaded': forbidden,
        'scope': 'Exact interval validation against the frozen molecule and proof that any positive ceiling gain is outside the prior full diagonal-spin cone. Prior upper-discovery cost remains separate.'}
    out.write_text(json.dumps(result, indent=2)+'\n'); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--separate', type=Path)
    args = parser.parse_args(); run(args.campaign, args.out, args.separate)
