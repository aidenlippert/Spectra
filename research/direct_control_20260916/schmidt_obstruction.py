"""Checked Schmidt-tail lower bounds without physical-state enumeration.

The lower bound concerns approximation in state-vector norm in a declared
orbital basis/order. It is not a lower bound on observable error or a universal
obstruction to control/reduced models.
"""
import argparse,json,time,math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.linalg import svd
from .validated_tensor import *
from research.correlated_pair_20260913.mps_exact import State
from research.intervention_reduction_20260916.exact import sqrt_up


def polar_error(a,columns=True):
 g,e=mm(a.conj().T,a) if columns else mm(a,a.conj().T)
 defect=plus(difference_norm(g,np.eye(g.shape[0])),e)
 if defect>=1:raise ArithmeticError('Cannot certify isometry correction')
 denominator=math.nextafter(1+math.nextafter(math.sqrt(math.nextafter(1-defect,-math.inf)),-math.inf),-math.inf)
 return up(defect/denominator)


def spectrum_bound(ts,cut):
 _,tails=tensor_bounds(ts)
 carry=np.ones((1,1),complex);prefix=1.;transformation_error=0.;left_polar=0.
 for i in range(cut):
  a=ts[i];raw,e1=mm(carry,a.reshape(a.shape[0],-1));raw=raw.reshape(carry.shape[0]*a.shape[1],a.shape[2])
  q,_=np.linalg.qr(raw,mode='reduced');nextcarry=q.conj().T@raw
  reconstructed,e2=mm(q,nextcarry)
  defect=plus(e1,difference_norm(raw,reconstructed),e2)
  transformation_error=plus(transformation_error,times(prefix,defect,tails[i+1]))
  on=opnorm(q);left_polar=plus(times(left_polar,on),polar_error(q));prefix=times(prefix,on);carry=nextcarry
 right_polar=0.;right_norm=1.
 for a in reversed(ts[cut:]):
  matrix=a.reshape(a.shape[0],-1);on=opnorm(matrix)
  right_polar=plus(times(right_polar,on),polar_error(matrix,False));right_norm=times(right_norm,on)
 state_to_center=plus(transformation_error,times(left_polar,frob(carry),right_norm),times(right_polar,frob(carry)))
 u,s,vh=svd(carry,full_matrices=False,check_finite=False)
 us=u*s[None,:];reconstructed,e=mm(us,vh)
 es=times(gamma(8),float(max(s)),frob(u))
 matrix_error=plus(difference_norm(carry,reconstructed),e,times(es,opnorm(vh)),
                   times(polar_error(u),frob(np.diag(s)),opnorm(vh)),times(polar_error(vh,False),frob(np.diag(s))))
 return s,plus(state_to_center,matrix_error),{'cut':cut,'center_shape':list(carry.shape),'canonical_map_error':state_to_center,'center_SVD_enclosure':matrix_error}


def lower_sqrt(x,bits=80):
 q=sqrt_up(x,scale=1<<bits)
 return max(F(),q-F(1,1<<bits))


def run(ts,source_error,norm2,ranks):
 result={r:F(0) for r in ranks};winning={};details=[]
 normup=sqrt_up(norm2,scale=1<<80)
 for cut in range(1,len(ts)):
  s,err,info=spectrum_bound(ts,cut);info['tails']={}
  allowance=F.from_float(plus(source_error,err))
  for rank in ranks:
   tail2=sum((F.from_float(float(v))**2 for v in s[rank:]),F())
   bound=max(F(),(lower_sqrt(tail2)-allowance)/normup)
   info['tails'][rank]=float(bound)
   if bound>result[rank]:result[rank]=bound;winning[rank]=cut
  details.append(info)
 return {'lower_bounds':{str(r):{'state_norm_error_at_least':str(result[r]),'float':float(result[r]),'witness_cut':winning.get(r)} for r in ranks},
         'details':details,'enumerated_configurations':0,
         'scope':'All MPS with every bond <= r in this declared orbital basis/order; normalized original input. State norm only. Not a control or observable impossibility claim.',
         'arithmetic':'Outward IEEE binary64 bounds plus exact rational final endpoints; BLAS rounding model assumptions apply.'}


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--basis',choices=['canonical','local'],required=True);p.add_argument('--out',required=True);args=p.parse_args()
 base=Path('results/direct_control_20260916');inputs=base/'imported/Spectra_control_reduction/inputs';start=time.monotonic()
 data=json.loads((inputs/'fixture.json').read_text());cert=json.loads((inputs/'state.json').read_text());state=State(data,cert)
 norm2=F(state.norm_integer,state.den**(2*state.m))
 if args.basis=='canonical':
  ts,error=load_mps(cert);ts,e=right_canonicalize(ts);error=plus(error,e)
 else:
  z=np.load(base/'local_basis_retry/state.npz');ts=[z[f'a{i}'] for i in range(data['modes'])]
  error=json.loads((base/'local_basis_retry/record.json').read_text())['source_to_rotated_state_error_bound']
 rec=run(ts,error,norm2,[8,16,24,32,48,64,96,128]);rec['basis']=args.basis;rec['seconds']=time.monotonic()-start
 out=Path(args.out)
 if out.exists():raise FileExistsError(out)
 out.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps({k:v for k,v in rec.items() if k!='details'}))
