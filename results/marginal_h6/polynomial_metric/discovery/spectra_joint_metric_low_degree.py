import json,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import numpy as np
from scipy.sparse import csc_matrix,bmat
from scipy.optimize import linprog
from experiments.marginal_polynomial_metric import JointPolynomial
from experiments.marginal_spin_constructor import spin_states
root=Path('/Users/aidenlippert/Documents/Spectra/results/marginal_h6/polynomial_metric');out=root/'joint_metric_low_degree_cone';out.mkdir(exist_ok=True);data=json.loads((root/'candidate.json').read_text());o=JointPolynomial(data);states=spin_states(o.o.oracle);sa=np.array(states,dtype=np.int64);d=np.array([sum(((s>>(2*i))&3)==3 for i in range(o.sites))-1 for s in states]);q=np.array([[((s>>(2*i))&3).bit_count()-1 for i in range(o.sites)] for s in states]);ionic=d>=0;orbits=json.loads((root/'with_constant/proposal.json').read_text())['orbits'];features=np.array([sum(np.prod(q**np.array(powers),axis=1) for powers in orbit) for orbit in orbits]).T.astype(float);labels=[];ri=[];ci=[];values=[];vcols=[];start=time.monotonic()
for degree in range(7):
 for bits in combinations(range(o.modes),degree):
  mask=sum(1<<i for i in bits)
  for assignment in range(1<<degree):
   occupied=sum(1<<i for j,i in enumerate(bits) if assignment&(1<<j));active=np.flatnonzero((sa&mask)==occupied)
   if len(active)>1 and degree<=4:
    j=len(labels);labels.append(['positive',mask,occupied]);ri.extend(active);ci.extend([j]*len(active));values.extend([1]*len(active))
    if degree<=4:vcols.append(j)
   if degree<=4:
    nonzero=active[d[active]!=0]
    if len(nonzero):
     j=len(labels);labels.append(['charge',mask,occupied]);ri.extend(nonzero);ci.extend([j]*len(nonzero));values.extend(d[nonzero].tolist())
     if degree<=2:vcols.append(j)
a=csc_matrix((np.array(values,dtype=float),(ri,ci)),shape=(len(states),len(labels)));va=a[:,vcols];index={s:i for i,s in enumerate(states)};h=np.zeros((len(states),len(states)))
for j,s in enumerate(states):
 for t,v in o.o.oracle.action(s).items():h[index[t],j]=float(v)
off=np.abs(h-np.diag(np.diag(h)));off[:,~ionic]=0
for gamma in [-6.24,-6.264]:
 b=(np.diag(np.diag(h)-gamma)-off)@features;one=np.ones((len(states),1));mean=np.mean(features[ionic],axis=0)[None,:]
 matrix=bmat([[csc_matrix(features),-va,None,None],[csc_matrix(b),None,-a,csc_matrix(-one)],[csc_matrix(mean),None,None,csc_matrix((1,1))]],format='csc');rhs=np.r_[np.full(len(states),.001),np.zeros(len(states)),1.];cost=np.zeros(matrix.shape[1]);cost[-1]=-1;bounds=[(None,None)]*features.shape[1]+[(0,None)]*(va.shape[1]+a.shape[1])+[(None,None)];start=time.monotonic();print(json.dumps({'phase':'solve','gamma':gamma,'shape':matrix.shape,'nnz':matrix.nnz}),flush=True)
 r=linprog(cost,A_eq=matrix,b_eq=rhs,bounds=bounds,method='highs-ipm',options={'time_limit':90});receipt={'gamma':gamma,'success':r.success,'status':r.message,'seconds':time.monotonic()-start,'metric_features':features.shape[1],'scope':'Joint metric and nonsingleton-cone discovery in explicit400-state polynomial quotient; both metric and numerator receive their own positivity decompositions. Mean Q metric fixed1. No exact certificate yet.'}
 if r.success:
  w=r.x[:features.shape[1]];v=features@w;k=b@w;receipt.update(bound=float(r.x[-1]),min_metric_Q=float(min(v[ionic])),min_numerator_Q=float(min(k[ionic])),max_equation_error=float(max(abs(matrix@r.x-rhs))))
  proposal={'orbits':orbits,'metric_coefficients':w.tolist(),'metric_bound':.001,'numerator_bound':float(r.x[-1]),'weight_labels':[labels[j] for j in vcols],'weight_values':r.x[features.shape[1]:features.shape[1]+va.shape[1]].tolist(),'numerator_labels':labels,'numerator_values':r.x[features.shape[1]+va.shape[1]:-1].tolist(),'gamma':str(gamma)};(out/('proposal_'+str(gamma)+'.json')).write_text(json.dumps(proposal,indent=2)+'\n')
 (out/('receipt_'+str(gamma)+'.json')).write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
 if r.success and r.x[-1]>.00001:break
