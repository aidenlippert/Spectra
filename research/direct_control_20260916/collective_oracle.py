"""Post-construction numerical validation ONLY: explicitly enumerates 4900.

Never imported by the exact constructor/replayer. This cost is separately
charged and is not evidence of cheap general dynamics.
"""
import json,time,sys
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
base=Path('results/direct_control_20260916');root=base/'imported/Spectra_control_reduction'
sys.path.insert(0,str(root/'code'))
from exact_control import prepare_input
start=time.monotonic();data=prepare_input(root)
cert=json.loads((base/'collective_control_bound64.json').read_text());p=cert['attempts'][-1]
n=len(data['D'])
def sparse(rows,den):
 indices=[];values=[];indptr=[0]
 for row in rows:
  for j,c in row.items():indices.append(j);values.append(c/den)
  indptr.append(len(indices))
 return csr_matrix((values,indices,indptr),shape=(n,n))
h=sparse(data['H'],data['den']);w=sparse(data['W'],1)
initial=np.array([float(F(x,data['pd'])) for x in data['p']]);initial/=np.linalg.norm(initial)
v=float(F(p['v_Ha']));T=float(F(p['T_au']))
generator=-1j*T*(h+v*w)
final=expm_multiply(generator,initial,traceA=generator.diagonal().sum())
D=float(np.sum(np.abs(final)**2*np.array(data['D'])));norm=float(np.linalg.norm(final))
lo,hi=map(lambda x:float(F(x)),p['D_interval'])
if not lo<=D<=hi or abs(norm-1)>1e-10:raise AssertionError('Numerical oracle disagrees')
rec={'status':'independent_full_sector_numerical_comparison_only','final_D':D,'norm':norm,'inside_certified_interval':True,
 'enumerated_configurations':n,'H_nonzeros':int(h.nnz),'seconds':time.monotonic()-start,
 'not_used_by_constructor_or_checker':True}
out=base/'collective_oracle.json'
if out.exists():raise FileExistsError(out)
out.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec))
