import json,time
from pathlib import Path
from itertools import product
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import qr
from experiments.marginal_spin_reduction import SpinZeroOracle
from experiments.marginal_spin_constructor import spin_states
root=Path('/Users/aidenlippert/Documents/Spectra');source=root/'results/marginal_h6/charge_spin/proof/certificate.json';c=json.loads(source.read_text());o=SpinZeroOracle(c);m=o.modes//2;states=[s for s in spin_states(o) if any(((s>>(2*i))&3)==3 for i in range(m))];idx={s:i for i,s in enumerate(states)};x=np.zeros((len(states),len(states)))
for j,s in enumerate(states):
 for t,v in o.action(s).items():
  if t in idx:x[idx[t],j]=float(v)
a=np.diag(np.diag(x)+6.264)-np.abs(x-np.diag(np.diag(x)));charges=np.array([[((s>>(2*i))&3).bit_count()-1 for i in range(m)] for s in states]);results=[]
out=root/'results/marginal_h6/polynomial_metric/with_constant';out.mkdir(parents=True,exist_ok=True)
for degree in [2,3,4,5,6]:
 orbits={}
 for powers in product(range(3),repeat=m):
  if sum(powers)>degree:continue
  canonical=min(powers,powers[::-1]);orbits.setdefault(canonical,[]).append(powers)
 orbits=list(orbits.values());f=np.column_stack([sum(np.prod(charges**np.array(powers),axis=1) for powers in orbit) for orbit in orbits]).astype(float);_,rr,piv=qr(f,mode="economic",pivoting=True);rank=int(np.sum(np.abs(np.diag(rr))>1e-8));selected=sorted(piv[:rank]);f=f[:,selected];orbits=[orbits[j] for j in selected];k=a@f;n=f.shape[1]
 objective=np.r_[np.zeros(n),-1.];ub=np.vstack([np.column_stack([-k,np.ones(len(states))]),np.column_stack([-f,np.zeros(len(states))])]);rhs=np.r_[np.zeros(len(states)),-np.full(len(states),1e-4)]
 start=time.monotonic();r=linprog(objective,A_ub=ub,b_ub=rhs,A_eq=[np.r_[np.mean(f,axis=0),0]],b_eq=[1],bounds=[(None,None)]*(n+1),method='highs')
 entry={'degree':degree,'features':n,'status':r.message,'seconds':time.monotonic()-start}
 if r.success:
  u=f@r.x[:-1];res=a@u;entry.update(raw_margin=float(min(res)),min_weight=float(min(u)),max_weight=float(max(u)),minimum_actual_row=float(min(res/u-6.264)),coefficient_l1=float(sum(abs(r.x[:-1]))))
  if min(res)>1e-7:
   proposal={'orbits':orbits,'coefficients':r.x[:-1].tolist(),'target_lower':'-6.264','charges_enumerated':140,'spin_rows_enumerated':380,'degree':degree,'scope':'Numerical polynomial-metric proposal using explicit Q rows; all reported margins reconstructed directly.'};(out/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n');results.append(entry);print(json.dumps(entry),flush=True);break
 results.append(entry);print(json.dumps(entry),flush=True)
(out/'degree_probe.json').write_text(json.dumps(results,indent=2)+'\n')
