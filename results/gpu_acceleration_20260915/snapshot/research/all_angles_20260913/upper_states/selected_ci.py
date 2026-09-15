"""Bounded selected-CI refinement from a Hamiltonian-only CISD seed."""
from pathlib import Path
import argparse
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_determinant_tree import DeterminantOracle
from research.certificate_scaling.streaming_reference_upper import upper

def run(fixture,out,add=200,rounds=1,max_states=None,scale=10**12):
 t0=time.monotonic(); c=json.loads(fixture.read_text()); o=DeterminantOracle(c)
 def de(s): return sum(float(v) for k,v in o.diagonal.items() if s&k==k)
 # Explicit polynomial-size reference: lowest diagonal determinant is found by enumeration only for this bounded fixture.
 allst=[s for s in range(1<<o.modes) if o.valid_state(s)]; hf=min(allst,key=de); occ=[i for i in range(o.modes) if hf>>i&1]; vir=[i for i in range(o.modes) if not hf>>i&1]
 basis={hf}
 for i in occ:
  for a in vir:basis.add(hf^(1<<i)^(1<<a))
 for x,i in enumerate(occ):
  for j in occ[x+1:]:
   for y,a in enumerate(vir):
    for b in vir[y+1:]:basis.add(hf^(1<<i)^(1<<j)^(1<<a)^(1<<b))
 def matdiag(B):
  B=sorted(B); idx={s:i for i,s in enumerate(B)}; M=np.zeros((len(B),len(B)))
  for j,s in enumerate(B):
   for z,v in o.action(s).items():
    if z in idx:M[idx[z],j]=float(v)
  ev,V=np.linalg.eigh(M);return B,ev[0],V[:,0],M
 seed_size=len(basis); total_candidates=0
 for _ in range(rounds):
  B,e,v,M=matdiag(basis); base={s:x for s,x in zip(B,v)}; scores={}
  for s,x in base.items():
   if abs(x)<1e-10:continue
   for z,h in o.action(s).items():
    if z not in base:scores[z]=scores.get(z,0)+abs(float(h*x)/(de(z)-e+1e-8))
  total_candidates+=len(scores); room=(max_states-len(basis)) if max_states else add
  chosen=sorted(scores,key=scores.get,reverse=True)[:min(add,room)]
  if not chosen:break
  basis.update(chosen); o.cache.clear()
  if max_states and len(basis)>=max_states:break
 B,e,v,M=matdiag(basis)
 amps=[int(round(float(x)*scale)) for x in v]; w={'states':[s for s,a in zip(B,amps) if a],'amplitudes':[a for a in amps if a]}; val,rep=upper(c,w)
 out.mkdir(parents=True,exist_ok=True);(out/'upper.json').write_text(json.dumps({'independent_upper':w,'upper':str(val)},indent=2)+'\n')
 rec={'kind':'hamiltonian_only_selected_ci_v2','fixture':str(fixture),'seed_cisd_basis':seed_size,'added':len(B)-seed_size,'rounds':rounds,'final_basis':len(B),'witness_support':len(w['states']),'upper':str(val),'upper_float':float(val),'discovery_seconds':time.monotonic()-t0,'matrix_entries':len(B)**2,'candidate_count':total_candidates,'legacy_reference_scan_states':len(allst),'replay':rep,'scope':'bounded projected diagonalization; dense cost is not claimed scalable','discovery_inputs':['fixture Hamiltonian only']}
 (out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n');return rec
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--fixture',type=Path);p.add_argument('--out',type=Path);p.add_argument('--add',type=int,default=200);p.add_argument('--rounds',type=int,default=1);p.add_argument('--max-states',type=int);a=p.parse_args()
 if a.fixture:
  print(json.dumps(run(a.fixture,a.out or Path('selected_ci_out'),a.add,a.rounds,a.max_states),indent=2))
 else:
  root=ROOT/'results/certificate_scaling/active_space_ladder';out=ROOT/'results/all_angles_20260913/upper_states';rows=[]
  for n in (6,8):rows.append(run(root/f'h{n}/fixture.json',out/f'h{n}_selected200',200))
  print(json.dumps(rows,indent=2))
