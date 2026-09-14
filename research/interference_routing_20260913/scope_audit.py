"""Exact support obstruction to the local template, in supplied orbital order."""
from pathlib import Path
import argparse
import hashlib
import json
from experiments.marginal_symbolic import decode

ROOT = Path(__file__).resolve().parents[2]


def obstruction(poly):
    outside = [(w, c) for w, c in poly.items() if w and
               max(i for _, i in w)-min(i for _, i in w) > 2]
    outside.sort(key=lambda item: (-abs(item[1]), item[0]))
    if not outside:
        return {'status': 'no_support_obstruction', 'is_template_acceptance': False}
    w, c = outside[0]
    return {'status': 'outside_template_span', 'outside_terms': len(outside),
            'witness_word': [[cr, i] for cr, i in w], 'witness_coefficient': str(c),
            'reason': 'Every declared template and diagonal term has support on at most three consecutive modes.',
            'scope': 'Supplied orbital ordering only; does not exclude another basis or a larger template family.'}


def run():
    records = []
    for sites in (4, 6, 8, 10):
        path = ROOT/f'results/certificate_scaling/active_space_ladder/h{sites}/fixture.json'
        raw = path.read_bytes()
        data = json.loads(raw)
        result = obstruction(decode(data['hamiltonian'], data['modes'], 4))
        records.append({'fixture': str(path), 'fixture_sha256': hashlib.sha256(raw).hexdigest(),
                        'modes': data['modes'], 'particles': data['particles'], **result})
    return records


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    result = run()
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps([{'modes': r['modes'], 'status': r['status'], 'outside_terms': r.get('outside_terms')} for r in result]))
