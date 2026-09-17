#!/usr/bin/env python3
"""Recheck all dynamics certificates and refusal tests. Standard library only."""
import argparse,json,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--skip-tests',action='store_true');a=p.parse_args()
out=Path(a.out).resolve()
if out.exists() or out==ROOT or ROOT in out.parents:raise SystemExit('Use a new output directory outside this bundle.')
out.mkdir(parents=True);start=time.monotonic()
if not a.skip_tests:
 with (out/'tests.log').open('w')as f:
  t=time.monotonic();proc=subprocess.run([sys.executable,'-B','-S',str(ROOT/'code/test_exact_control.py')],stdout=f,stderr=subprocess.STDOUT)
  if proc.returncode:raise SystemExit('Exact tests failed; inspect tests.log')
 tests_time=time.monotonic()-t
else:tests_time=None
sys.path.insert(0,str(ROOT/'code'))
from exact_control import run
r=run(['short24','long64','baseline32','short8_refusal'],out/'dynamics')
required={k:v['target_proved']for k,v in r['results'].items()}
if required!={'short24':True,'long64':True,'baseline32':True,'short8_refusal':False}:raise SystemExit('Unexpected acceptance/refusal outcomes')
res={'tests_seconds':tests_time,'dynamics':r,'total_seconds':time.monotonic()-start,'expected_claims_and_refusal_reproduced':True}
(out/'complete.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps({'targets':required,'total_seconds':res['total_seconds']},indent=2))
