"""Fresh standard-library replay of all operator, energy, and rank proofs."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from research.molecular_collective_20260913.core import tail_replay,support_replay,rank_replay

ROOT=Path(__file__).resolve().parents[2]


def run(campaign,out):
    start=time.monotonic();campaign=campaign.resolve();rows=[];counts={'tail':0,'support':0,'rank':0}
    for directory in sorted(p for p in campaign.iterdir() if p.is_dir()):
        data=json.loads((directory/'fixture.json').read_text())
        for path in sorted(directory.glob('rank_*/tail.json'))+[directory/'uncentered_tail.json']:
            cert=json.loads(path.read_text());receipt=tail_replay(data,cert)
            previous=json.loads((path.parent/'receipt.json').read_text())['tail'] if path.name=='tail.json' else json.loads((directory/'uncentered_receipt.json').read_text())
            if receipt['tail_interval_width_Ha']!=previous['tail_interval_width_Ha']:raise ValueError('Tail replay changed')
            rows.append({'path':str(path.relative_to(ROOT)),'kind':'tail','seconds':receipt['replay_seconds'],
                         'width_Ha':receipt['tail_interval_width_Ha'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest()});counts['tail']+=1
            support=path.parent/'support.json'
            if path.name=='tail.json' and support.exists():
                rr=support_replay(data,cert,json.loads(support.read_text()))
                previous=json.loads((path.parent/'support_receipt.json').read_text())['support']
                for key in ('certified_lower_Ha','support_family_ceiling_Ha'):
                    if rr[key]!=previous[key]:raise ValueError('Supporting-bound replay changed')
                rows.append({'path':str(support.relative_to(ROOT)),'kind':'support','seconds':rr['replay_seconds'],
                             'lower_Ha':rr['certified_lower_Ha'],'ceiling_Ha':rr['support_family_ceiling_Ha'],
                             'sha256':hashlib.sha256(support.read_bytes()).hexdigest()});counts['support']+=1
        path=directory/'rank_obstruction.json';rr=rank_replay(data,json.loads(path.read_text()))
        rows.append({'path':str(path.relative_to(ROOT)),'kind':'rank','seconds':rr['replay_seconds'],
                     'minimum_factors':rr['minimum_factors'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest()});counts['rank']+=1
    numerical=sorted(n for n in sys.modules if n.split('.')[0] in ('numpy','scipy','cvxpy','pyscf'))
    if numerical:raise AssertionError('Numerical package loaded by accepting replay')
    sources={}
    for module in tuple(sys.modules.values()):
        filename=getattr(module,'__file__',None)
        if filename:
            path=Path(filename).resolve()
            if path.is_relative_to(ROOT) and path.suffix=='.py':sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    result={'counts':counts,'rows':rows,'wall_seconds':time.monotonic()-start,
            'site_packages_disabled':bool(sys.flags.no_site),'forbidden_imports':numerical,
            'source_sha256':dict(sorted(sources.items()))}
    out.write_text(json.dumps(result,indent=2)+'\n');print(counts,result['wall_seconds']);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--campaign',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();run(args.campaign,args.out)
