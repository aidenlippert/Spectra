#!/usr/bin/env python3
"""Independent complete replay. Standard library only; no stored energy is trusted.
Run: python -B -S replay.py --out /a/new/output/directory
"""
from pathlib import Path
from fractions import Fraction
import argparse,json,sys,time,hashlib
BASE=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 if a.out.exists():raise FileExistsError('Use a fresh output directory; sealed receipts are never overwritten.')
 a.out.mkdir(parents=True);start=time.monotonic()
 manifest=json.loads((BASE/'SHA256.json').read_text())
 for name,want in manifest.items():
  path=BASE/name
  if hashlib.sha256(path.read_bytes()).hexdigest()!=want:raise ValueError('Changed bundle file: '+name)
 sys.path.insert(0,str(BASE/'source'))
 from research.correlated_pair_20260913.mps_exact import check as check_upper
 from research.collective_completion_20260914.spin_screen import check as check_lower
 read=lambda p:json.loads(p.read_text())
 fixture=read(BASE/'inputs/fixture.json')
 print('Checking the actual rational MPS upper...',flush=True)
 upper=check_upper(fixture,read(BASE/'inputs/state.json'))
 (a.out/'upper.json').write_text(json.dumps(upper,indent=2)+'\n')
 print('Checking new singlet and inherited nonsinglet proofs...',flush=True)
 lower=check_lower(fixture,read(BASE/'certificates/singlet.json'),read(BASE/'certificates/nonsinglet.json'))
 U=Fraction(upper['upper_Ha']);L=Fraction(lower['lower'])
 if U<L:raise ValueError('Inconsistent exact endpoints')
 result={'lower_Ha':str(L),'upper_Ha':str(U),'width_Ha':str(U-L),'width_mHa':float(1000*(U-L)),'target_1p6mHa_met':U-L<=Fraction(1,625),'upper_replay':upper,'lower_replay':lower,'total_replay_seconds':time.monotonic()-start}
 if any(k in sys.modules for k in ('numpy','scipy','cvxpy','quimb','pyscf')):raise AssertionError('Numerical library loaded on accepting path')
 expected=read(BASE/'RESULT.json')
 if str(L)!=expected['lower_Ha'] or str(U)!=expected['upper_Ha']:raise ValueError('Computed endpoints differ from this bundle; inspect the new receipts')
 (a.out/'complete.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k] for k in ('lower_Ha','width_mHa','target_1p6mHa_met','total_replay_seconds')},indent=2))
 if not result['target_1p6mHa_met']:raise SystemExit('Accuracy target not met')
if __name__=='__main__':main()
