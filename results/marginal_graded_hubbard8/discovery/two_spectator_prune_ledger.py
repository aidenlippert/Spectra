"""Retain all accepted physical sources and determinant anchors in a proposal ledger."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[3]


def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    parser.add_argument('accepted', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    receipt = args.accepted / 'range_two_family_limit_replay.json'
    cert = args.accepted / 'range_two_family_limit_certificate.json'
    r = json.loads(receipt.read_text())
    if not r.get('accepted') or not r.get('two_spectator_hopping'):
        raise ValueError('Accepted two-spectator family required')
    for name, h in r['source_sha256'].items():
        if digest(ROOT / name) != h:
            raise ValueError('Stale accepted family source')
    c = json.loads(cert.read_text())
    canonical = lambda v: tuple(sorted((int(s), a) for s, a in v.items() if a))
    required = {canonical(item['vector']) for item in c['mixture']}
    path = args.source / 'candidate_ledger.json.gz'
    with gzip.open(path, 'rt') as f:
        ledger = json.load(f)
    count = len(ledger['states'])
    if not 1 <= count <= 20000 or any(len(ledger[k]) != count for k in ('rows', 'energies')):
        raise ValueError('Bounded complete candidate ledger required')
    retained, found, seen = [], set(), set()
    for i, state in enumerate(ledger['states']):
        key = canonical(state['vector'])
        if key in required:
            found.add(key)
        if (key in required or state['sector'] is None) and key not in seen:
            retained.append(i)
            seen.add(key)
    if found != required:
        raise ValueError('Ledger does not contain every accepted physical source')
    if len(retained) > 13000:
        raise ValueError('Pruned ledger exceeds bounded resume input')
    args.output.mkdir(parents=True, exist_ok=False)
    for name in ('profile_joint_r1_2_certificate.json', 'congruence_witnesses.json'):
        shutil.copy2(args.source / name, args.output / name)
    pruned = dict(ledger)
    for k in ('states', 'rows', 'energies'):
        pruned[k] = [ledger[k][i] for i in retained]
    with gzip.open(args.output / 'candidate_ledger.json.gz', 'wt') as f:
        json.dump(pruned, f, separators=(',', ':'))
    files = [Path(__file__), path, receipt, cert, args.source / 'profile_joint_r1_2_certificate.json']
    out = {'accepted': False, 'preservation_checked': True, 'input_candidates': count,
           'retained_candidates': len(retained), 'retained_accepted_sources': len(required),
           'retained_indices': retained,
           'source_sha256': {str(p.resolve().relative_to(ROOT)): digest(p) for p in files},
           'scope': 'Exact vector identity preserves every source in the independently accepted mixture. Numerical rows remain untrusted; this is not an energy or family proof.'}
    (args.output / 'pruning_diagnostic.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: v for k, v in out.items() if k not in ('source_sha256', 'retained_indices')}))


if __name__ == '__main__':
    main()
