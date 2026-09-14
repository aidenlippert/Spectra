"""Fresh standard-library replay of lower certificates and reference uppers."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from research.certificate_scaling.streaming_reference_upper import upper
from research.joint_patterns_20260913.core import replay

ROOT = Path(__file__).resolve().parents[2]


def run(campaign, out):
    start = time.monotonic(); summary = json.loads((campaign/'summary.json').read_text())
    rows = []; references = {}
    for case in summary['cases']:
        if not case['within_budget_accepted']:
            continue
        fixture = ROOT/case['fixture_path']; tail_path = ROOT/case['tail_path']
        if hashlib.sha256(fixture.read_bytes()).hexdigest() != case['fixture_file_sha256'] or hashlib.sha256(tail_path.read_bytes()).hexdigest() != case['tail_file_sha256']:
            raise ValueError('Frozen campaign input hash changed')
        data = json.loads(fixture.read_text()); tail = json.loads(tail_path.read_text())
        path = campaign/case['case']/'certificate.json'
        cert = json.loads(path.read_text()); accepted = replay(data, tail, cert)
        if accepted['original_lower_Ha'] != case['original_lower_Ha']:
            raise ValueError('Fresh exact bound disagrees with discovery receipt')
        name = case['system']
        if name not in references:
            reference_path = ROOT/f'results/certificate_scaling/active_space_ladder_references_aligned/{name}/upper.json'
            reference = json.loads(reference_path.read_text())
            value, stats = upper(data, reference['independent_upper'])
            if value != F(reference['upper']):
                raise ValueError('Reference witness does not reproduce its upper bound')
            references[name] = {'upper_Ha': str(value), 'receipt': stats,
                'source_path': str(reference_path.relative_to(ROOT)),
                'source_sha256': hashlib.sha256(reference_path.read_bytes()).hexdigest()}
        width = F(references[name]['upper_Ha'])-F(accepted['original_lower_Ha'])
        if width < 0:
            raise ValueError('Lower certificate exceeded independent exact upper')
        row = {'case': case['case'], 'fresh': accepted,
            'reference_upper_Ha': references[name]['upper_Ha'], 'interval_width_Ha': str(width),
            'interval_width_mHa': float(1000*width), 'meets_1_6_mHa': width <= F(16, 10000),
            'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
        rows.append(row)
        print(case['case'], 'interval_mHa', row['interval_width_mHa'], flush=True)
    forbidden = sorted(name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules)
    if forbidden:
        raise AssertionError('Numerical package loaded on the accepting path')
    result = {'rows': rows, 'references': references, 'fresh_wall_seconds': time.monotonic()-start,
        'numerical_packages_loaded': forbidden, 'scope': 'Exact lower and reference-upper replay. Prior upper-witness discovery is outside the compact lower-bound method.'}
    out.write_text(json.dumps(result, indent=2)+'\n'); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); run(args.campaign, args.out)
