"""Run the exterior residual replay on unchanged local-family certificates."""
import json
from fractions import Fraction as F
from pathlib import Path
from research.certificate_scaling.wedge_spectral_bound import propose, replay_residual
from experiments.marginal_symbolic import verified_residual, verify
ROOT=Path(__file__).resolve().parents[3]
def run(tag):
 p=ROOT/'results/global_response_20260913/joint_factor'/tag/'certificate.json'; c=json.loads(p.read_text())
 residual,_=verified_residual(c)
 if max(map(len,residual))>6: raise ValueError('refuse residual above body order 3')
 witness=propose(residual,c['modes'],c['particles'])
 spectral=replay_residual(c,residual,witness)
 return {'case':tag,'raw_b':str(F(c['b'])),'original_lower':verify(c)['lower'],
         'spectral_lower':spectral['lower'],'residual_terms':len(residual),
         'residual_degrees':sorted(set(map(len,residual))),
         'certificate_unchanged':True,'many_body_states_enumerated':0}
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True); a=p.parse_args()
 rows=[run('frozen_h6'),run('fresh_h6_1p73')]; a.out.write_text(json.dumps(rows,indent=2)+'\n'); print(json.dumps(rows,indent=2))
