"""Numerical positivity discovery in an exact coefficient quotient, without states."""
import argparse,json,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb
from unittest.mock import patch
import numpy as np
from scipy.sparse import csc_matrix
from scipy.optimize import linprog
from experiments.marginal_polynomial_metric import JointPolynomial
from experiments.marginal_number_quotient import NumberSliceQuotient
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--seed',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--time-limit',type=float,default=60);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);c=json.loads(a.source.read_text());seed=json.loads(a.seed.read_text());assert c['polynomial_metric']==seed['polynomial_metric'];start=time.monotonic();ring=JointPolynomial(c)
with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No states')),patch('experiments.marginal_spin_constructor.spin_states',side_effect=AssertionError('No sector generation')):
 _,poly,cost=ring.compile(F(c['target_lower']))
degrees=[min(ring.target,max(5,max(((m&spin).bit_count() for m in poly),default=0),max(((a['required']&spin).bit_count() for a in seed['numerator_proof']['positive_indicators']),default=0))) for spin in ring.spin_masks];qs=[NumberSliceQuotient(list(range(s,ring.modes,2)),ring.target,degrees[s]) for s in (0,1)];basis=[x|y for x in qs[0].monomials if x not in qs[0].rows for y in qs[1].monomials if y not in qs[1].rows];index={m:i for i,m in enumerate(basis)};cache={}
def normal_monomial(mask):
 if mask not in cache:
  alpha,_=qs[0].monomial(mask&ring.spin_masks[0]);beta,_=qs[1].monomial(mask&ring.spin_masks[1]);cache[mask]={x|y:u*v for x,u in alpha.items() for y,v in beta.items()}
 return cache[mask]
def normal(poly):
 out={}
 for m,v in poly.items():
  if not v:continue
  for k,x in normal_monomial(m).items():out[k]=out.get(k,F(0))+x*v
 return {m:v for m,v in out.items() if v}
def indicator(required,occupied):
 result={occupied:1}
 for i in range(ring.modes):
  if (required^occupied)&(1<<i):result=ring.multiply(result,{0:1,1<<i:-1})
 return result
labels=[['b',0,0]];ri=[index[0]];ci=[0];values=[1.];known=set();charge={0:-1,**{3<<(2*i):1 for i in range(ring.sites)}}
def include(family,required,occupied):
 key=(family,required,occupied)
 if key in known:return
 known.add(key);term=indicator(required,occupied)
 if family=='charge':term=ring.multiply(term,charge)
 reduced=normal(term)
 if not reduced:return
 j=len(labels);labels.append(list(key))
 for m,v in reduced.items():ri.append(index[m]);ci.append(j);values.append(float(v))
# Seed directions are inherited from the ORIGINAL-H proof only; their weights are discarded.
for item in seed['numerator_proof']['positive_indicators']:include('positive',item['required'],item['occupied'])
for degree in range(5):
 for chosen in combinations(range(ring.modes),degree):
  required=sum(1<<i for i in chosen)
  for assignment in range(1<<degree):
   occupied=sum(1<<i for j,i in enumerate(chosen) if assignment&(1<<j));count=1
   for spin in ring.spin_masks:
    available=ring.sites-(required&spin).bit_count();needed=ring.target-(occupied&spin).bit_count();count*=comb(available,needed) if 0<=needed<=available else 0
   if not count:continue
   if count>1:include('positive',required,occupied)
   include('charge',required,occupied)
rhs=normal({m:F(v,cost['numerator_scale']) for m,v in poly.items()});b=np.zeros(len(basis))
for m,v in rhs.items():b[index[m]]=float(v)
A=csc_matrix((values,(ri,ci)),shape=(len(basis),len(labels)));objective=np.zeros(len(labels));objective[0]=-1;prepared=time.monotonic();print(json.dumps({'phase':'coefficient_matrix','rows':A.shape[0],'columns':A.shape[1],'nonzeros':A.nnz,'preparation_seconds':prepared-start}),flush=True)
r=linprog(objective,A_eq=A,b_eq=b,bounds=[(None,None)]+[(0,None)]*(len(labels)-1),method='highs-ipm',options={'time_limit':a.time_limit})
receipt={'success':bool(r.success),'status':r.message,'solve_seconds':time.monotonic()-prepared,'preparation_seconds':prepared-start,'rows':A.shape[0],'columns':A.shape[1],'nonzeros':A.nnz,'quotients':[q.stats for q in qs],'configuration_evaluations':0,'seed':str(a.seed),'source':str(a.source),'scope':'Numerical coefficient-space proposal using inherited original-H positive directions plus all degree<=4 indicators/localizers. No physical state evaluation; metric inherited; exact export is required for acceptance.'}
if r.success:
 receipt['bound']=float(r.x[0]);receipt['minimum_coefficient']=float(min(r.x[1:]));receipt['max_equation_residual']=float(max(abs(A@r.x-b)));proposal={'bound':float(r.x[0]),'labels':[l for l,x in zip(labels[1:],r.x[1:]) if abs(x)>1e-13],'values':[float(x) for x in r.x[1:] if abs(x)>1e-13]};(a.out/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
(a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
