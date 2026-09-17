"""Exact rational orbital circuit applied directly to the supplied MPS."""
from fractions import Fraction as F
from itertools import permutations
from math import lcm,gcd
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import svd
from research.transfer_followup_20260915.rotated_upper import check_rotation
from .validated_tensor import *
from .orbital_order import group_spatial,split_spatial


def determinant(a):
 n=len(a)
 if not n:return F(1)
 total=F()
 for p in permutations(range(n)):
  inv=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
  v=F((-1)**inv)
  for i,j in enumerate(p):v*=a[i][j]
  total+=v
 return total


def local_fock_gate(c,s):
 # G rows new, cols old. mode order alpha-left,beta-left,alpha-right,beta-right.
 u=[[F(0)]*4 for _ in range(4)]
 for p,q in ((0,2),(1,3)):
  u[p][p]=u[q][q]=c;u[p][q]=s;u[q][p]=-s
 # C-order physical indices: |n0,n1,n2,n3>.
 occupied=[[k for k in range(4) if idx>>(3-k)&1] for idx in range(16)]
 gate=[[F(0)]*16 for _ in range(16)]
 for i,a in enumerate(occupied):
  for j,b in enumerate(occupied):
   if len(a)==len(b):gate[i][j]=determinant([[u[p][q] for q in b] for p in a])
 return gate


def givens_design(target):
 A=np.array(target,float);n=len(A)
 R=[[F(int(i==j)) for j in range(n)] for i in range(n)];gates=[]
 for j in range(n-1):
  for k in range(n-1,j,-1):
   a,b=A[k-1,j],A[k,j];h=math.hypot(a,b)
   if h<1e-15:continue
   c,s=a/h,b/h
   if c<0:c,s=-c,-s
   t=F(float(s/(1+c))).limit_denominator(10**8)
   c=(1-t*t)/(1+t*t);s=2*t/(1+t*t)
   r0=R[k-1][:];r1=R[k][:]
   R[k-1]=[c*x+s*y for x,y in zip(r0,r1)];R[k]=[-s*x+c*y for x,y in zip(r0,r1)]
   v0=A[k-1].copy();v1=A[k].copy()
   A[k-1]=float(c)*v0+float(s)*v1;A[k]=-float(s)*v0+float(c)*v1
   gates.append((k-1,c,s))
 signs=[1 if A[j,j]>=0 else -1 for j in range(n)]
 U=[[R[j][i]*signs[j] for j in range(n)] for i in range(n)]
 den=lcm(*(x.denominator for row in U for x in row));Z=[[int(x*den) for x in row] for row in U]
 common=gcd(den,*(abs(x) for row in Z for x in row));den//=common;Z=[[x//common for x in row] for row in Z]
 check_rotation(Z,den)
 return gates,signs,Z,den


def apply_gate(ts,i,gate):
 pref=prefix_gram_bounds(ts);_,suf=tensor_bounds(ts)
 a,b=ts[i:i+2];l=a.shape[0];r=b.shape[2]
 z,e1=mm(a.reshape(l*4,-1),b.reshape(b.shape[0],4*r))
 z=z.reshape(l,16,r).transpose(1,0,2).reshape(16,l*r)
 g=np.zeros((16,16),complex);er2=0.
 for x in range(16):
  for y in range(16):
   g[x,y],err=rational_float(gate[x][y]);er2=plus(er2,times(err,err))
 gg,e2=mm(g,z)
 arith=plus(times(opnorm(g),e1),e2,times(root(er2),plus(frob(z),e1)))
 matrix=gg.reshape(16,l,r).transpose(1,0,2).reshape(l*4,4*r)
 u,s,vh=svd(matrix,full_matrices=False,check_finite=False);v=u.conj().T@matrix
 rec,re=mm(u,v);de=plus(arith,difference_norm(matrix,rec),re)
 out=list(ts);out[i]=u.reshape(l,4,-1);out[i+1]=v.reshape(u.shape[1],4,r)
 return out,times(pref[i],de,suf[i+2])


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',required=True);args=p.parse_args()
 folder=Path(args.out);folder.mkdir(parents=True,exist_ok=False)
 start=time.monotonic();base=Path('results/direct_control_20260916')
 target=json.loads((base/'local_basis_design.json').read_text())['target']
 gates,signs,Z,den=givens_design(target)
 rotation={'integer_matrix':Z,'denominator':str(den),'gates':[[i,str(c),str(s)] for i,c,s in gates],'signs':signs,
           'exact_orthogonality_checked':True,'design':'Symmetric orthogonalized atomic orbital proposal; original rational coefficients authoritative'}
 (folder/'rotation.json').write_text(json.dumps(rotation,separators=(',',':'))+'\n')
 cert=json.loads((base/'imported/Spectra_control_reduction/inputs/state.json').read_text())
 ts,error=load_mps(cert);ts,e=right_canonicalize(ts);error=plus(error,e);ts,e=group_spatial(ts);error=plus(error,e)
 for j,(i,c,s) in enumerate(gates):
  ts,e=apply_gate(ts,i,local_fock_gate(c,s));error=plus(error,e)
  if j%7==6:print(json.dumps({'gate':j+1,'error':error,'seconds':time.monotonic()-start}),flush=True)
 for i,sign in enumerate(signs):
  ts[i]=ts[i].copy();ts[i][:,1:3,:]*=sign
 ts,e=split_spatial(ts);error=plus(error,e);ts,e=right_canonicalize(ts);error=plus(error,e)
 np.savez_compressed(folder/'state.npz',**{f'a{i}':x for i,x in enumerate(ts)})
 rec={'source_to_rotated_state_error_bound':error,'gates':len(gates),'seconds_before_compression':time.monotonic()-start,'compressions':[]}
 for bond in [8,16,24,32,48,64,96]:
  bs,e,_=compress(ts,bond);r={'bond':bond,'compression_error_bound':e,'total_error_bound':plus(e,error),'norm':expectation(bs)}
  rec['compressions'].append(r);print(json.dumps(r),flush=True)
 rec['total_seconds']=time.monotonic()-start
 (folder/'record.json').write_text(json.dumps(rec,indent=2)+'\n')
