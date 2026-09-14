import json,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import numpy as np
from scipy.sparse import csc_matrix
from scipy.optimize import linprog
from experiments.marginal_polynomial_metric import JointPolynomial
from experiments.marginal_spin_constructor import spin_states
root=Path('/Users/aidenlippert/Documents/Spectra/results/marginal_h6/polynomial_metric');out=root/'cubic_sos_cutting';out.mkdir(exist_ok=True)
c=json.loads((root/'candidate.json').read_text());o=JointPolynomial(c);compiled=json.loads((root/'compiled.json').read_text());states=spin_states(o.o.oracle);sa=np.array(states,dtype=np.int64);d=np.array([sum(((s>>(2*i))&3)==3 for i in range(o.sites))-1 for s in states]);labels=[['b']];ri=list(range(len(states)));ci=[0]*len(states);values=[1]*len(states);start=time.monotonic()
for degree in range(7):
 for bits in combinations(range(o.modes),degree):
  mask=sum(1<<i for i in bits)
  for assignment in range(1<<degree):
   occupied=sum(1<<i for j,i in enumerate(bits) if assignment&(1<<j));active=np.flatnonzero((sa&mask)==occupied)
   if len(active)>1:
    j=len(labels);labels.append(['positive',mask,occupied]);ri.extend(active);ci.extend([j]*len(active));values.extend([1]*len(active))
   if degree<=4:
    nonzero=active[d[active]!=0]
    if len(nonzero):
     j=len(labels);labels.append(['charge',mask,occupied]);ri.extend(nonzero);ci.extend([j]*len(nonzero));values.extend(d[nonzero].tolist())
a=csc_matrix((np.array(values,dtype=float),(ri,ci)),shape=(len(states),len(labels)));exact={m:int(v) for m,v in compiled['numerator']};scale=compiled['receipt']['numerator_scale'];p=[sum(v for m,v in exact.items() if m&s==m) for s in states];objective=np.zeros(len(labels));objective[0]=-1

from scipy.sparse import hstack

seed=json.loads((root/'full_sos/degree3_proposal.json').read_text());basis=seed['basis'];f=np.array([[(s&m)==m for m in basis] for s in states],dtype=float);history=[];squares=[]
def add_square(vector):
 global a
 coefficients=np.round(vector*10**6).astype(np.int64);values=(f@coefficients/10**6)**2;a=hstack([a,csc_matrix(values[:,None])],format='csc');labels.append(['sos',len(squares)]);squares.append({'denominator':10**6,'terms':[{'mask':m,'coefficient':int(v)} for m,v in zip(basis,coefficients) if v]})
_,vectors=np.linalg.eigh(np.array(seed['gram']))
for j in range(vectors.shape[1]):add_square(vectors[:,j])
for iteration in range(8):
 objective=np.zeros(a.shape[1]);objective[0]=-1;start=time.monotonic();r=linprog(objective,A_eq=a,b_eq=np.array([float(F(v,scale)) for v in p]),bounds=[(None,None)]+[(0,None)]*(a.shape[1]-1),method='highs-ipm',options={'time_limit':45})
 entry={'iteration':iteration,'status':r.message,'seconds':time.monotonic()-start,'columns':a.shape[1]}
 if not r.success:history.append(entry);print(json.dumps(entry),flush=True);break
 z=-r.eqlin.marginals;moment=f.T@(z[:,None]*f);eigen,vectors=np.linalg.eigh(moment);entry.update(bound=float(r.x[0]),minimum_moment_eigenvalue=float(eigen[0]),basis_degree=3);history.append(entry);print(json.dumps(entry),flush=True)
 checkpoint={'iteration':iteration,'basis':basis,'squares':squares,'dual_states':states,'dual_weights':z.tolist(),'bound':float(r.x[0])};(out/'checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n');(out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
 if r.x[0]>1e-5:
  (out/'proposal.json').write_text(json.dumps({'labels':[l for l,x in zip(labels,r.x) if abs(x)>1e-12],'values':[float(x) for x in r.x if abs(x)>1e-12],'squares':squares},indent=2)+'\n');break
 if eigen[0]>=-1e-8:break
 for index in np.flatnonzero(eigen<-1e-8)[:16]:add_square(vectors[:,index])
(out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
