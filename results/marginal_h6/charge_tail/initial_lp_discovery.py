import json,time
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix
from experiments.marginal_coherent_tree import CoherentCharge
from experiments.marginal_charge_polynomial import add,scale,multiply,indicator
root=Path('/Users/aidenlippert/Documents/Spectra');data=json.loads((root/'results/marginal_h6/charge_spin/proof/certificate.json').read_text());o=CoherentCharge(data);poly=dict(o.oracle.diagonal);ratios=[]
for (c,a),amp in o.groups.items():
 closed=o.close(c|a,a)
 if closed is None:continue
 lo,hi=o.amplitude_range(amp,*closed)
 if lo>=0:sign=1
 elif hi<=0:sign=-1
 else:raise ValueError('unstable sign')
 delta=tuple(((c>>(2*i))&3).bit_count()-((a>>(2*i))&3).bit_count() for i in range(o.sites))
 ratio=o.ratio_upper(delta,*closed);ratios.append(float(ratio))
 poly=add(poly,scale(multiply(indicator(c|a,a),amp),-sign*ratio))
print('compiled',len(poly),'max ratio',max(ratios),flush=True)
shifts=[{0:F(-o.target),**{1<<i:F(1) for i in range(spin,o.modes,2)}} for spin in (0,1)]
charge={0:F(-2),**{3<<(2*i):F(1) for i in range(o.sites)}}
masks=[sum(1<<i for i in inds) for k in range(5) for inds in combinations(range(o.modes),k)]
columns=[{0:F(1)}];bounds=[(None,None)];labels=[['b']]
for mask in masks:
 bits=[1<<i for i in range(o.modes) if mask&(1<<i)]
 for assignment in range(1<<len(bits)):
  occupied=sum(bit for j,bit in enumerate(bits) if assignment&(1<<j));ip=indicator(mask,occupied)
  columns.append(ip);bounds.append((0,None));labels.append(['positive',mask,occupied])
  if len(bits)<=2:columns.append(multiply(charge,ip));bounds.append((0,None));labels.append(['charge',mask,occupied])
 if len(bits)<=3:
  for spin,shift in enumerate(shifts):columns.append(multiply(shift,{mask:F(1)}));bounds.append((None,None));labels.append(['ideal',spin,mask])
rows=sorted(set(poly).union(*(set(x) for x in columns)));idx={m:i for i,m in enumerate(rows)};ri=[];ci=[];values=[]
for j,col in enumerate(columns):
 for mask,value in col.items():ri.append(idx[mask]);ci.append(j);values.append(float(value))
a=csc_matrix((values,(ri,ci)),shape=(len(rows),len(columns)));obj=np.zeros(len(columns));obj[0]=-1;t=time.monotonic()
r=linprog(obj,A_eq=a,b_eq=[float(poly.get(m,0)) for m in rows],bounds=bounds,method='highs',options={'time_limit':120})
out=root/'results/marginal_h6/charge_tail';out.mkdir(exist_ok=True)
receipt={'success':r.success,'status':r.message,'bound':float(r.x[0]) if r.x is not None else None,'rows':len(rows),'columns':len(columns),'seconds':time.monotonic()-t,'minimum_doublons':2,'scope':'Numerical proposal only: grouped fixed-sign weighted row envelope using independent metric-ratio upper bounds; degree4 Boolean positivity LP on D>=2.'}
(out/'initial_probe.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
if r.success:(out/'initial_lp_solution.json').write_text(json.dumps({'labels':labels,'values':r.x.tolist()},indent=2)+'\n')
