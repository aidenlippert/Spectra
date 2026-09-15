"""Fresh-process standard-library replay of all persisted energy certificates."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from research.collective_interference_20260913.collective import replay

ROOT=Path(__file__).resolve().parents[2]


def run(campaign,out):
    start=time.monotonic();summary=json.loads((campaign/'summary.json').read_text());rows=[]
    for original in summary['rows']:
        directory=campaign/original['name']
        model_path=directory/'model.json';certificate_path=directory/'certificate.json'
        result=replay(json.loads(model_path.read_text()),json.loads(certificate_path.read_text()))
        for key in ('lower','upper','width'):
            if result[key]!=original['accepted'][key]:raise ValueError('Persisted endpoint does not replay')
        rows.append({'name':original['name'],'lower':result['lower'],'upper':result['upper'],
                     'width':result['width'],'seconds':result['replay_seconds'],
                     'model_file_bytes':model_path.stat().st_size,'certificate_file_bytes':certificate_path.stat().st_size,
                     'model_file_sha256':hashlib.sha256(model_path.read_bytes()).hexdigest(),
                     'certificate_file_sha256':hashlib.sha256(certificate_path.read_bytes()).hexdigest()})
    forbidden=sorted(name for name in sys.modules if name.split('.')[0] in ('numpy','scipy','cvxpy','pyscf'))
    if forbidden:raise AssertionError('Numerical package imported in accepting replay')
    sources={}
    for module in tuple(sys.modules.values()):
        source=getattr(module,'__file__',None)
        if source:
            path=Path(source).resolve()
            if path.is_relative_to(ROOT) and path.suffix=='.py':
                sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    result={'rows':rows,'accepted':len(rows),'wall_seconds':time.monotonic()-start,
            'site_packages_disabled':bool(sys.flags.no_site),'forbidden_numerical_imports':forbidden,
            'python_version':sys.version,'loaded_local_source_sha256':dict(sorted(sources.items()))}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('accepted','wall_seconds','site_packages_disabled','forbidden_numerical_imports')}))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--campaign',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args();run(args.campaign,args.out)
