"""Exact lower/upper replay for the frozen local-family runs."""
import json, time
from fractions import Fraction as F
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from experiments.marginal_symbolic import verify
from experiments.marginal_determinant_tree import DeterminantOracle

CASES=[
 ('frozen_h6', ROOT/'results/global_response_20260913/joint_factor/frozen_h6/certificate.json', ROOT/'results/certificate_scaling/active_space_ladder_references_aligned/h6/upper.json'),
 ('fresh_h6_1p73', ROOT/'results/global_response_20260913/joint_factor/fresh_h6_1p73/certificate.json', ROOT/'results/global_response_20260913/fresh_h6_1p73/reference_upper.json')]
def main(out):
 rows=[]
 for name,cp,rp in CASES:
  cert=json.loads(cp.read_text()); ref=json.loads(rp.read_text()); upper=ref.get('independent_upper',ref)
  t=time.monotonic(); lo=F(verify(cert)['lower']); hi=DeterminantOracle(cert).upper(upper)
  row={'case':name,'lower':str(lo),'upper':str(hi),'interval_width':str(hi-lo),'interval_width_float':float(hi-lo),
       'upper_support':len(upper['states']),'many_body_states_enumerated_lower':0,
       'lower_replay_seconds':time.monotonic()-t,'certificate':str(cp),'reference':str(rp),
       'family':'local','ideal_body':2,'symmetry':True,'width':4,'rounds':5,
       'terminal_rule_constructed':True,'full_interval_passes_1p6mHa':hi-lo<=F(16,10000)}
  rows.append(row)
 out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(rows,indent=2)+'\n'); return rows
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);print(json.dumps(main(p.parse_args().out),indent=2))
