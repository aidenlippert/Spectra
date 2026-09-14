"""Independent CAR closure check for an accepted two-spectator family cap."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import sys
from two_spectator_overlap_probe import LABELS,operator

ROOT=Path(__file__).resolve().parents[3]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    folder=parser.parse_args().directory.resolve()
    cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json'
    c=json.loads(cp.read_text());r=json.loads(rp.read_text())
    if not r.get('accepted') or not r.get('two_spectator_hopping'):
        raise ValueError('Accepted enlarged two-spectator family required')
    for name,h in r['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=h:raise ValueError('Stale family source')
    moments={}
    for label in LABELS:
        action=operator(label);moment=F(0)
        for item in c['mixture']:
            v={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in v.values())
            moment+=F(item['weight'])*F(sum(a*b*v.get(t,0) for s,a in v.items() for t,b in action[s].items()),norm)
        if moment:raise ValueError('Independent two-spectator moment is nonzero')
        moments[','.join(map(str,label))]='0'
    files={Path(__file__).resolve(),Path(__file__).with_name('two_spectator_overlap_probe.py'),cp,rp}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name and str(Path(name).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(name).resolve())
    result={'accepted':True,'two_spectator_moments':moments,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
            'scope':'Independent full-Fock CAR reconstruction confirms all30 two-spectator hopping moments vanish in this accepted local mixture. Necessary stationary consistency checks only; no general quantum extension or representability claim.'}
    (folder/'two_spectator_closure.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'accepted':True,'zero_moments':30}))


if __name__=='__main__':main()
