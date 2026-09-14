"""Fresh stdlib-only replay of the inherited and all accepted new intervals."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from research.certificate_scaling.streaming_reference_upper import upper
from research.spin_completion_20260913.core import replay

ROOT = Path(__file__).resolve().parents[2]


def run(campaign, out):
    start = time.monotonic(); prior = ROOT/'results/molecular_collective_20260913/campaign/h6'
    data = json.loads((prior/'fixture.json').read_text()); tail = json.loads((prior/'rank_10/tail.json').read_text())
    reference_path = ROOT/'results/certificate_scaling/active_space_ladder_references_aligned/h6/upper.json'
    reference = json.loads(reference_path.read_text()); U, reference_receipt = upper(data, reference['independent_upper'])
    if U != F(reference['upper']):
        raise ValueError('Reference upper did not reproduce')
    first = json.loads((campaign/'continue/summary.json').read_text()); source = Path(first['source'])
    inherited_cert = json.loads((source/'certificate.json').read_text())
    inherited = replay(data, tail, inherited_cert)
    if inherited['original_lower_Ha'] != first['source_replay']['original_lower_Ha']:
        raise ValueError('Inherited exact lower changed')
    inherited_lower = F(inherited['original_lower_Ha']); rows = []
    for arm in ('continue', 'separate'):
        directory = campaign/arm
        if not (directory/'summary.json').exists():
            continue
        summary = json.loads((directory/'summary.json').read_text()); source = Path(summary['source'])
        for filename, expected in summary['source_hashes'].items():
            if hashlib.sha256((source/filename).read_bytes()).hexdigest() != expected:
                raise ValueError('A bound or proposal source changed after the campaign')
        for trial in summary['trials']:
            if not trial['accepted'] or trial['complete_at_seconds'] > summary['budget_seconds']:
                continue
            path = directory/trial['case']/'certificate.json'; cert = json.loads(path.read_text())
            receipt = replay(data, tail, cert); lower = F(receipt['original_lower_Ha']); width = U-lower
            if lower != F(trial['lower_Ha']) or width < 0:
                raise ValueError('Fresh interval mismatch')
            record = {'arm': arm, 'case': trial['case'], 'directions': trial['directions'],
                'lower_Ha': str(lower), 'upper_Ha': str(U), 'width_Ha': str(width), 'width_mHa': float(1000*width),
                'gain_over_inherited_mHa': float(1000*(lower-inherited_lower)), 'meets_1_6_mHa': width <= F(16, 10000),
                'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'fresh': receipt}
            rows.append(record); print(json.dumps({k: record[k] for k in ('arm', 'case', 'directions', 'width_mHa', 'gain_over_inherited_mHa')}), flush=True)
    forbidden = [name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules]
    if forbidden:
        raise AssertionError('Numerical package loaded on exact acceptance')
    result = {'inherited': inherited, 'inherited_width_mHa': float(1000*(U-inherited_lower)),
        'rows': rows, 'reference_upper': reference_receipt,
        'reference_sha256': hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        'wall_seconds': time.monotonic()-start, 'numerical_packages_loaded': forbidden,
        'scope': 'Exact finite-H6 intervals with the same collective tail. No cone optimum or scaling claim.'}
    out.write_text(json.dumps(result, indent=2)+'\n'); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True); args = parser.parse_args()
    run(args.campaign, args.out)
