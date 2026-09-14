"""Independent exact closure check for the fourteen spectator moments."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from spectator_hopping_basis import LABELS,operator


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);folder=parser.parse_args().directory.resolve();cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json';r=json.loads(rp.read_text());c=json.loads(cp.read_text())
    if not r['accepted'] or not r.get('spectator_hopping') or r['source_sha256'][str(cp.relative_to(ROOT))]!=hashlib.sha256(cp.read_bytes()).hexdigest():raise ValueError('Accepted spectator family required')
    moments={}
    for label in LABELS:
        action=operator(label);moment=F(0)
        for item in c['mixture']:
            v={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in v.values());moment+=F(item['weight'])*F(sum(a*b*v.get(t,0) for s,a in v.items() for t,b in action[s].items()),norm)
        if moment:raise ValueError('Independent spectator moment does not vanish')
        moments[','.join(map(str,label[:3]))]=str(moment)
    files={Path(__file__).resolve(),Path(__file__).with_name('spectator_hopping_basis.py'),cp,rp}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'spectator_moments':moments,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Independent full-Fock CAR reconstruction confirms all14 nonconstant one-spectator charge-hopping moments vanish. These necessary stationary consistency conditions are not sufficient for general quantum representability.'}
    (folder/'one_spectator_overlap.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'accepted':True,'zero_moments':len(moments)}))


if __name__=='__main__':main()
