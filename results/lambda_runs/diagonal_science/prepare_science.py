"""Current physical tensors and LP-directed extra-overlap candidates for GPU screening."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from results.marginal_graded_hubbard8.discovery.joint_profile_numeric import build
from experiments.marginal_charged_projectors import charged_vectors
BASE=ROOT/'results/marginal_graded_hubbard8/joint_projector/signed_density';OUT=Path(__file__).resolve().parent


def main():
 start=time.monotonic();path=BASE/'symmetric_diagonals_8/profile_joint_r1_2_certificate.json';c=json.loads(path.read_text());local=c['local_window']
 x=[F(local[k][i]) for k,i in [('onsite_profile',0),('onsite_profile',1),('hopping_profile',0),('hopping_profile',1),('density_profile',0),('density_profile',1)]]
 ss=build(x);pert=[]
 for j in range(6):
  y=x.copy();y[j]+=F(1,1000);pert.append(build(y))
 half={int(s):v for s,v in c['vector'].items()};hn=sum(v*v for v in half.values());cv,cn=charged_vectors(c['joint']['vector'])
 Ys=[{int(s):v for s,v in item['diagonal'].items()} for item in json.loads((BASE/'symmetry_diagonal_candidates.json').read_text())['candidates']]
 actual={int(s):F(v) for s,v in c['telescoping_diagonal'].items()};allowed=[j for j in range(8,len(Ys)) if len(set(actual)|set(Ys[j]))<=64]
 alpha,beta=F(c['penalty']),F(c['joint']['penalty']);ratio=F(c['joint']['ratio']);th=F(c['projector_sum_ceiling'])/c['windows'];tj=F(c['joint']['projector_sum_ceiling'])/c['joint']['windows']
 mats=[];raw=[];eigen=[];arrays={}
 for i,(key,(_,K,cols)) in enumerate(ss.items()):
  sc=np.sqrt([sum(v*v for v in col.values()) for col in cols]);A=np.array(K,float)/sc[:,None]/sc[None,:]
  Ds=np.array([(np.array(p[key][1],float)/sc[:,None]/sc[None,:]-A)*1000 for p in pert])
  def P(v,n):
   z=np.array([sum(a*v.get(s,0) for s,a in col.items()) for col in cols],float)/sc/np.sqrt(float(n));return np.outer(z,z)
  ph=P(half,hn);pc=sum((P(v,cn) for v in cv),np.zeros_like(A))
  T=np.array([[Y.get(next(iter(col))&1023,0)-Y.get(next(iter(col))>>2,0) for col in cols] for Y in Ys],float)
  diag=np.array([float(actual.get(next(iter(col))&1023,0)-actual.get(next(iter(col))>>2,0)) for col in cols])
  A+=float(alpha+beta)*ph+float(beta*ratio)*pc+np.diag(diag)
  D=np.concatenate((Ds,ph[None],(ph+float(ratio)*pc)[None]))
  vals,vec=eigh(A,subset_by_index=[0,min(2,len(A)-1)])
  mats.append((key,A,D,T));eigen.append((vals,vec));raw.append(float(vals[0]))
  arrays[f'A{i}']=A;arrays[f'D{i}']=D;arrays[f'T{i}']=T
 floor=min(raw);gradient=[];gaps=[]
 for (_,A,D,T),(vals,vecs) in zip(mats,eigen):
  for j,e in enumerate(vals):
   if e<=floor+1e-4:
    w=vecs[:,j];gradient.append(np.r_[[w@d@w for d in D],T@(w*w)]);gaps.append(float(e-floor))
 gradients=np.array(gradient);points_dense=[np.zeros(8)];points_diag=[np.zeros(len(Ys))];labels=[{'extra':None,'step':0}]
 bounds=[(-.02,.02),(-.08,.08),(-.02,.02),(-.02,.02),(-.05,.05),(-.05,.05),(-min(.02,float(alpha)),.02),(-min(.02,float(beta)),.02)]+[(-.05,.05)]*9+[(None,None)]
 cost=np.r_[np.zeros(6),float(th),float(tj),np.zeros(9),-1.]
 for extra in allowed:
  columns=list(range(8))+list(range(8,16))+[8+extra]
  G=gradients[:,columns];res=linprog(cost,A_ub=np.c_[-G,np.ones(len(G))],b_ub=np.array(gaps),bounds=bounds,method='highs')
  if not res.success:continue
  for step in (.01,.03,.1,.2,.3,.5,1.):
   dense=step*res.x[:8];diagonal=np.zeros(len(Ys));diagonal[:8]=step*res.x[8:16];diagonal[extra]=step*res.x[16]
   if x[0]+x[1]+sum(dense[:2])>10 or x[2]+x[3]+sum(dense[2:4])>2.5 or min(np.array(x[:4],float)+dense[:4])<0:continue
   points_dense.append(dense);points_diag.append(diagonal);labels.append({'extra':extra,'step':step,'lp_predicted_density_gain':float(-res.fun/5)})
 if len(points_dense)>2048:raise ValueError('Screening exceeds2048 candidate cap')
 arrays['points_dense']=np.array(points_dense);arrays['points_diag']=np.array(points_diag);arrays['theta']=np.array([float(th),float(tj)]);arrays['base_penalties']=np.array([float(alpha),float(beta)])
 # Fresh physical matrix check for a perturbed candidate; not only the affine formula.
 check=min(4,len(points_dense)-1);newx=[v+F(round(float(d)*10**9),10**9) for v,d in zip(x,points_dense[check][:6])]
 # Physical profile entries have denominator1e6 caps; use a separate exactly representable probe.
 offsets=[F(i-2,1000000) for i in range(6)];fresh=build([v+d for v,d in zip(x,offsets)])
 discrepancy=0.
 for key,(_,K,cols) in fresh.items():
  sc=np.sqrt([sum(a*a for a in col.values()) for col in cols]);oldK=ss[key][1]
  difference=(np.array(K,float)-np.array(oldK,float))/sc[:,None]/sc[None,:]
  item=next(item for item in mats if item[0]==key)
  discrepancy=max(discrepancy,float(np.max(abs(difference-np.einsum('i,ijk->jk',np.array(offsets,float),item[2][:6])))))
 if discrepancy>1e-10:raise ValueError('Physical affine assembly disagrees')
 np.savez_compressed(OUT/'science.npz',**arrays)
 meta={'certificate':str(path.relative_to(ROOT)),'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'data_sha256':hashlib.sha256((OUT/'science.npz').read_bytes()).hexdigest(),'sectors':len(mats),'basis_shapes':len(Ys),'eligible_extra_shapes':len(allowed),'candidates':len(points_dense),'active_eigenvectors':len(gaps),'baseline_minimum':floor,'baseline_periodic_density':(floor-float(alpha*th+beta*tj))/5,'fresh_matrix_max_error':discrepancy,'preparation_seconds':time.monotonic()-start,'labels':labels,'scope':'Numerical LP-directed screening of one extra symmetric diagonal overlap shape, retaining at most64 nonzero five-site diagonal entries. Current v5 certificate and physical matrices rebuilt locally.'}
 (OUT/'preparation.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps({k:v for k,v in meta.items() if k!='labels'}),flush=True)


if __name__=='__main__':main()
