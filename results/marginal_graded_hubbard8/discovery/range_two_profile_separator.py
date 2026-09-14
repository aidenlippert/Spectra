"""Exact free-profile separator in a fixed-profile range-two dual mixture."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]


def charges(s, sites):
    return [((s >> (2*i)) & 3).bit_count()-1 for i in range(sites)]


def five_site_y(s):
    q = charges(s, 5)
    return q[0]*q[2]-q[2]*q[4]


def derivative(s):
    q = charges(s, 6)
    return q[0]*q[2]-q[1]*q[3]-q[2]*q[4]+q[3]*q[5]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    folder = parser.parse_args().directory.resolve()
    certificate_path = folder/'range_two_family_limit_certificate.json'
    receipt_path = folder/'range_two_family_limit_replay.json'
    receipt = json.loads(receipt_path.read_text())
    if not receipt['accepted'] or receipt['source_sha256'][str(certificate_path.relative_to(ROOT))] != hashlib.sha256(certificate_path.read_bytes()).hexdigest():
        raise ValueError('Exact accepted mixture required')
    for s in range(4096):
        if derivative(s) != five_site_y(s & 1023)-five_site_y(s >> 2):
            raise ValueError('Physical profile derivative is not the claimed telescope')
    c = json.loads(certificate_path.read_text())
    expectation = F(0)
    for item in c['mixture']:
        vector = {int(s): a for s, a in item['vector'].items()}
        expectation += F(item['weight'])*F(sum(a*a*derivative(s) for s, a in vector.items()), sum(a*a for a in vector.values()))
    result = {'accepted': True, 'violated': expectation != 0, 'exact_derivative_expectation': str(expectation),
              'derivative_expectation_float': float(expectation), 'identity_states_checked': 4096,
              'five_site_diagonal_nonzero_entries': sum(five_site_y(s) != 0 for s in range(1024)),
              'compact_Y': 'q0*q2-q2*q4', 'six_site_profile_direction': [1, -1, -1, 1],
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(), certificate_path, receipt_path]},
              'scope': 'Exact profile derivative and telescoping identity. A nonzero expectation excludes this local dual from translation-consistent marginals and from a dual that cancels arbitrary mean-correct range-two profiles. The fixed-profile family cap remains valid. No optimized stronger energy bound is claimed.'}
    (folder/'free_profile_separator.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'exact_derivative_expectation')}), flush=True)


if __name__ == '__main__':
    main()
