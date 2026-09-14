import json,time
from pathlib import Path
from itertools import combinations
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix
from experiments.marginal_polynomial_metric import JointPolynomial
root=Path('/Users/aidenlippert/Documents/Spectra/results/marginal_h6/polynomial_metric');data=json.loads((root/'candidate.json').read_text());o=JointPolynomial(data);compiled=json.loads((root/'compiled.json').read_text());modes=o.modes;degree=6
out=Path('/tmp/audit_polynomial_nine_completion');out.mkdir(exist_ok=True)
charge={0:-1,**{3<<(2*i):1 for i in range(o.sites)}};shifts=[{0:-o.target,**{1<<i:1 for i in range(s,modes,2)}} for s in (0,1)]
columns=[{0:1}];labels=[['b']];bounds=[(None,None)];masks=[sum(1<<i for i in inds) for k in range(degree+1) for inds in combinations(range(modes),k)]
start=time.monotonic()
for mask in masks:
 bits=[1<<i for i in range(modes) if mask&(1<<i)]
 for assignment in range(1<<len(bits)):
  occ=sum(bit for j,bit in enumerate(bits) if assignment&(1<<j));poly={occ:1}
  if any((occ&s).bit_count()>o.target for s in o.spin_masks):continue
  for bit in bits:
   if not occ&bit:poly=o.multiply(poly,{0:1,bit:-1})
  if not poly:continue
  if len(bits)<=4 or (len(bits)==6 and sum(bool(mask&(1<<i)) for i in range(0,12,2))==3 and sum(bool(mask&(1<<i)) for i in range(1,12,2))==3 and sum(bool(occ&(1<<i)) for i in range(0,12,2)) in (1,2) and sum(bool(occ&(1<<i)) for i in range(1,12,2)) in (1,2)) :
   columns.append(poly);labels.append(['positive',mask,occ]);bounds.append((0,None))
  if len(bits)<=degree-2:columns.append(o.multiply(charge,poly));labels.append(['charge',mask,occ]);bounds.append((0,None))
 if len(bits)<=degree-1 and all((mask&s).bit_count()<=o.target for s in o.spin_masks):
  for spin,shift in enumerate(shifts):columns.append(o.multiply(shift,{mask:1}));labels.append(['ideal',spin,mask]);bounds.append((None,None))
rows=sorted(set().union(*(set(c) for c in columns)));index={m:i for i,m in enumerate(rows)};ri=[];ci=[];values=[]
for j,col in enumerate(columns):
 for mask,value in col.items():ri.append(index[mask]);ci.append(j);values.append(value)
a=csc_matrix((np.asarray(values,dtype=float),(ri,ci)),shape=(len(rows),len(columns)));objective=np.array([0 if label[0] in ['b','ideal'] else (100 if label[0]=='positive' and label[1].bit_count()==6 else 1) for label in labels],dtype=float)
print(json.dumps({'phase':'matrix','rows':len(rows),'columns':len(columns),'nnz':len(values),'seconds':time.monotonic()-start}),flush=True)
for name,scale in [('numerator',compiled['receipt']['numerator_scale'])]:
 poly={mask:int(value)/scale for mask,value in compiled[name]};bounds[0]=((0.01,0.01) if name=='weight' else (0.001,0.001));start=time.monotonic();r=linprog(objective,A_eq=a,b_eq=[poly.get(m,0) for m in rows],bounds=bounds,method='highs-ipm',options={'time_limit':120})
 receipt={'success':r.success,'status':r.message,'bound':float(r.x[0]) if r.x is not None else None,'seconds':time.monotonic()-start,'rows':len(rows),'columns':len(columns),'nonzero_proposal':int(np.count_nonzero(abs(r.x)>1e-10)) if r.x is not None else None,'degree':degree,'scope':'Numerical global fixed-spin,D>=1 positivity proposal; coefficient space, no occupation row constraints.'}
 (out/(name+'_lp_receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(phase=name,**receipt)),flush=True)
 if r.success:(out/(name+'_lp_solution.json')).write_text(json.dumps({'labels':[label for label,x in zip(labels,r.x) if abs(x)>1e-12],'values':[float(x) for x in r.x if abs(x)>1e-12]},indent=2)+'\n')
