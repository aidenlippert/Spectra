#!/usr/bin/env python3
"""Rerun the numerical H8-specific proposal procedure in a new external directory.

Requires existing NumPy, SciPy and Numba; installs nothing. This reuses the
provided rational Hamiltonian and MPS, not molecular integrals. The algorithm is
not a general molecule dispatcher. It enumerates the balanced-spin sector.
"""
from pathlib import Path
import argparse,sys,shutil,subprocess,json,time
p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args()
ROOT=Path(__file__).resolve().parent;out=Path(a.out).resolve()
if out.exists()or out==ROOT or ROOT in out.parents:raise SystemExit('Use a new directory outside the bundle')
(out/'inputs').mkdir(parents=True);(out/'development').mkdir()
shutil.copytree(ROOT/'code',out/'code',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
for f in ['fixture.json','state.json']:shutil.copy2(ROOT/'inputs'/f,out/'inputs'/f)
steps=['propose.py','propose_short.py','adapt.py','long_control_probe.py','long_adapt.py','baseline.py','export_trajectory.py','numerical_validation.py']
ledger=[]
for script in steps:
 t=time.monotonic()
 with (out/'development'/f'{script}.log').open('w')as f:proc=subprocess.run([sys.executable,'-B',str(out/'code'/script)],stdout=f,stderr=subprocess.STDOUT)
 ledger.append({'script':script,'process_seconds':time.monotonic()-t,'returncode':proc.returncode})
 (out/'rediscovery_costs.json').write_text(json.dumps(ledger,indent=2))
 if proc.returncode:raise SystemExit(f'{script} failed; output and costs preserved')
print('Numerical proposals generated. They require independent exact_control.py acceptance; numerical logs are not certificates.')
