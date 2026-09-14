"""Independently replay the two projector moments of an accepted new mixture."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
from full_overlap_telescope import ROOT,build,moment

def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    directory=parser.parse_args().directory.resolve()
    cp=directory/'range_two_family_limit_certificate.json';rp=directory/'range_two_family_limit_replay.json'
    c=json.loads(cp.read_text());r=json.loads(rp.read_text())
    if not r.get('accepted') or not r.get('coherent_projector'):raise ValueError('Accepted enlarged family required')
    for name,h in r['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=h:raise ValueError('Stale family receipt')
    moments={}
    for label,vector in [('358,601,-1',{358:1,601:-1}),('346,613,1',{346:1,613:1})]:
        denominator,_,matrix=build(vector);value=moment(c['mixture'],matrix,denominator)
        if value or F(r['family_replay']['coherent_projector_moments'][label]):raise ValueError('New projector moment did not close')
        moments[label]=str(value)
    files={Path(__file__).resolve(),cp,rp,Path(sys.modules['full_overlap_telescope'].__file__).resolve(),Path(sys.modules['full_overlap_density'].__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name and str(Path(name).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(name).resolve())
    result={'accepted':True,'coherent_projector_moments':moments,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
            'scope':'Both independent exact projector-telescope contractions are zero. This does not establish full RDM overlap agreement or a quantum extension.'}
    (directory/'coherent_projector_closure.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'accepted':True,'moments':moments}))

if __name__=='__main__':main()
