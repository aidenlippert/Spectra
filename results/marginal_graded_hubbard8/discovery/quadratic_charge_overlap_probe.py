"""Exact compact quadratic charge separators in an accepted free-profile dual."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]


def charges(s, sites):
    return [((s >> (2*i)) & 3).bit_count()-1 for i in range(sites)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    folder = parser.parse_args().directory.resolve()
    certificate_path = folder/'range_two_family_limit_certificate.json'
    receipt_path = folder/'range_two_family_limit_replay.json'
    receipt = json.loads(receipt_path.read_text())
    if not receipt['accepted'] or not receipt.get('free_range_two_profile') or receipt['source_sha256'][str(certificate_path.relative_to(ROOT))] != hashlib.sha256(certificate_path.read_bytes()).hexdigest():
        raise ValueError('An accepted free-profile mixture is required')
    c = json.loads(certificate_path.read_text())
    directions = []
    for i in range(5):
        for j in range(i, 5):
            pair, reflected = (i, j), (4-j, 4-i)
            if pair >= reflected:
                continue
            def y(s):
                q = charges(s, 5)
                return q[i]*q[j]-q[reflected[0]]*q[reflected[1]]
            moment = F(0)
            for item in c['mixture']:
                vector = {int(s): a for s, a in item['vector'].items()}
                moment += F(item['weight'])*F(sum(a*a*(y(s & 1023)-y(s >> 2)) for s, a in vector.items()), sum(a*a for a in vector.values()))
            if pair != (0, 3) and moment:
                raise ValueError('An already-canceled quadratic direction is nonzero')
            directions.append({'pair': pair, 'reflected_pair': reflected,
                               'exact_moment': str(moment), 'moment_float': float(moment),
                               'five_site_nonzero_entries': sum(y(s) != 0 for s in range(1024))})
    result = {'accepted': True, 'directions': directions,
              'violated_directions': sum(F(row['exact_moment']) != 0 for row in directions),
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(), certificate_path, receipt_path]},
              'scope': 'Exact moments of every reflection-odd quadratic charge polynomial on five sites. Each Y_left-Y_right sums to zero under translation. Five directions are implied by existing nearest/free-range-two profiles and cancel exactly. Nonzero remaining moments are compact necessary consistency separators, not accepted energy improvements.'}
    (folder/'quadratic_charge_overlap.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps([{k: v for k, v in row.items() if k != 'exact_moment'} for row in directions]), flush=True)


if __name__ == '__main__':
    main()
