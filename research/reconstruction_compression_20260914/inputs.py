"""Preserve inherited bytes and freeze exact endpoint budgets."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/reconstruction_compression_20260914'
PARENT=ROOT/'results/correlated_pair_20260913'
STATES={'h6':'h6_b48_real','h8':'h8_spatial_warm144'}

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,obj): Path(path).write_text(json.dumps(obj,indent=2)+'\n')

def initialize():
    inherited=json.loads((PARENT/'preservation_before.json').read_text())['files']
    inherited.update(json.loads((PARENT/'manifest.json').read_text())['files'])
    p=PARENT/'manifest.json'
    inherited[str(p.relative_to(ROOT))]={'bytes':p.stat().st_size,'sha256':sha(p)}
    errors=[p for p,v in inherited.items() if not (ROOT/p).is_file() or sha(ROOT/p)!=v['sha256']]
    if errors: raise ValueError(('Inherited seal mismatch',errors))
    dump(OUT/'preservation_before.json',{'files':inherited,'checked_files':len(inherited)})
    cases={}
    for case,stem in STATES.items():
        folder=PARENT/'mps'/stem
        old=json.loads((folder/'interval.json').read_text())
        meta=json.loads((ROOT/f'results/certificate_scaling/cubic_precision/intervals/final_{case}.json').read_text())
        U,L=F(old['upper_Ha']),F(old['lower_Ha'])
        paths={'fixture':ROOT/f'results/certificate_scaling/active_space_ladder/{case}/fixture.json',
               'state':folder/'state.json','upper_receipt':folder/'interval.json',
               'teacher':ROOT/meta['certificate'],'teacher_residual':ROOT/meta['proof']}
        cases[case]={**{k:str(p.relative_to(ROOT)) for k,p in paths.items()},
                     'sha256':{k:sha(p) for k,p in paths.items()},
                     'upper_Ha':str(U),'lower_Ha':str(L),'required_lower_Ha':str(U-F(1,625)),
                     'lower_deterioration_budget_Ha':str(L-U+F(1,625)),
                     'width_mHa':float(1000*(U-L))}
    dump(OUT/'frozen_inputs.json',cases)
    print(json.dumps({'preserved_files':len(inherited),'cases':{k:{'width_mHa':v['width_mHa'],
          'slack_mHa':float(1000*F(v['lower_deterioration_budget_Ha']))} for k,v in cases.items()}}))

if __name__=='__main__': initialize()
