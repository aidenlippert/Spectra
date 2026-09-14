"""Untrusted frozen-coefficient ablation; remove only the full diagonal spin-word telescope."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,sys
from spin_word_frozen_target import ROOT,minimum

def main():
    parser=argparse.ArgumentParser();parser.add_argument('seed',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    c=json.loads(args.seed.read_text())
    if c.get('kind')!='hubbard_projector_extension_v18' or not c.get('spin_word_telescope'):
        raise ValueError('Explicit new-operator source required')
    c.pop('spin_word_telescope');c['kind']='hubbard_projector_extension_v17'
    value=minimum(c);c['penalized_lower']=str(F(math.floor(value*10**7)-1,10**7))
    (args.output/'profile_joint_r1_2_certificate.json').write_text(json.dumps(c,indent=2)+'\n')
    files={Path(__file__).resolve(),args.seed.resolve(),Path(sys.modules['spin_word_frozen_target'].__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name and str(Path(name).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(name).resolve())
    result={'accepted':False,'numerical_minimum':value,'proposed_threshold':c['penalized_lower'],
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
            'scope':'Removed only spin_word_telescope and lowered schema to v17; all other coefficients fixed, scalar threshold reproposed. Requires exact replay.'}
    (args.output/'ablation_proposal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'accepted':False,'proposed_threshold':c['penalized_lower']}))

if __name__=='__main__':main()
