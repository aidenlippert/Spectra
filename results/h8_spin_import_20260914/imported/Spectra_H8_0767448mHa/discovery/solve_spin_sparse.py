"""Boundary-point discovery in a sparse highest-spin-weight representation.
Physical dictionaries, Hamiltonian and exact accepting checker are unchanged.
"""
from common import *
from scipy.sparse.linalg import splu
import argparse,resource
class SpinOperator:
    def __init__(self):
        old=load_operator();self.meta=old.meta;self.prepared=old.prepared;self.rhs=old.rhs;self.free=old.free;self.scale=old.scale
        p=OUT/'spin_sparse';rep=json.loads((p/'representation.json').read_text());vz=np.load(p/'representation.npz');z=np.load(p/'seed.npz');self.info=rep
        self.V=[vz[f'V_{i}'] for i in range(len(rep['blocks']))];self.members=[[b['physical_group']] for b in rep['blocks']];self.ids=[0]*len(self.V)
        self.M=[sparse.load_npz(p/f'map_{i}.npz').tocsr() for i in range(len(self.V))];self.MT=[s.T.tocsr() for s in self.M]
        self.Q=[z[f'Q_{i}'] for i in range(len(self.V))];self.x=z['x'];self.y=z['y'];self.G=sparse.load_npz(p/'normal.npz').tocsc()
    def A(self,Q):return sum((s@q.ravel() for s,q in zip(self.M,Q)),np.zeros(len(self.rhs)))
    def AT(self,y):return [np.asarray(s@y).reshape(v.shape[1],v.shape[1]) for s,v in zip(self.MT,self.V)]

def run():
 p=argparse.ArgumentParser();p.add_argument('tag');p.add_argument('--seconds',type=float,default=700);p.add_argument('--mu',type=float,default=.3);p.add_argument('--restart');p.add_argument('--overrelax',type=float,default=1.);a=p.parse_args()
 out=OUT/a.tag;out.mkdir(exist_ok=False);start=time.monotonic();op=SpinOperator();Q=[positive(q) for q in op.Q];x=op.x.copy();y=-np.load(OUT/'physical_moment_dual.npy')
 if a.restart:
  z=np.load(a.restart);Q=[z[f'Q_{i}'] for i in range(len(op.Q))];x=z['x'];y=z['y']
 c=np.zeros(op.free.shape[1]);c[0]=-1.;mu=a.mu;Z=[positive(-z) for z in op.AT(y)]
 lu=splu(op.G+sparse.eye(len(op.rhs),format='csc')*1e-11);pre=LinearOperator(op.G.shape,matvec=lu.solve)
 rng=np.random.default_rng(51);v=rng.normal(size=len(op.rhs));err=np.linalg.norm(op.G@v-(op.A(op.AT(v))+op.free@(op.free.T@v)))/np.linalg.norm(op.G@v);assert err<1e-10
 print('SETUP',json.dumps({'seconds':time.monotonic()-start,'gram_entries':sum(q.size for q in Q),'coefficient_map_nonzeros':sum(s.nnz for s in op.M),'normal_nonzeros':op.G.nnz,'normal_map_error':err,'method':'sparse_highest_weight_boundary_point','full_Gram_teacher_used':False,'many_body_states_enumerated':0,'warm_start':'uploaded compact candidate'}),flush=True)
 predictor=FastResidual(op);hist=[];best={'score':-1e100};last=0;end=time.monotonic()+a.seconds
 # Numerical extraction of the scalar/ideal coordinates optimizes exports, not the cone.
 quo=None
 def record(it,extra=None):
  nonlocal quo
  # remove any drift in ideal variables using the source rank-revealing convention
  if quo is None:quo=Quotient(op)
  xx=quo.recover(Q)
  res=(op.A(Q)+op.free@xx-op.rhs)/op.scale;corr,cc=predictor.evaluate(res);score=float(xx[0])+corr-2.36e-10
  dm=min(float(linalg.eigvalsh(-z,subset_by_index=(0,0),check_finite=False)[0]) for z in op.AT(y));pe=float(np.linalg.norm(op.rhs-op.A(Q)-op.free@x));de=float(np.sqrt(sum(np.sum((z+t)**2) for z,t in zip(op.AT(y),Z))+np.linalg.norm(op.free.T@y-c)**2))
  rec={'iteration':it,'seconds':time.monotonic()-start,'mu':mu,'b':float(xx[0]),'raw_width_mHa':1000*(U-xx[0]),'numerical_residual_lower_Ha':corr,'unverified_width_mHa':1000*(U-score),'selected_l1':float(np.sum(abs(res))),'primal_l2':pe,'dual_l2':de,'untrusted_dual_b':float(-op.rhs@y),'minimum_dual_eigenvalue':dm,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
  if extra:rec.update(extra)
  hist.append(rec);record_json(out/'history.json',hist);np.savez_compressed(out/'checkpoint.npz',y=y,x=x,**{f'Q_{i}':q for i,q in enumerate(Q)})
  if score>best['score']:
   best.update(score=score,record=rec,x=xx.copy(),Q=[q.copy() for q in Q]);np.savez_compressed(out/'best.npz',y=y,x=xx,**{f'Q_{i}':q for i,q in enumerate(Q)})
  print(json.dumps(rec),flush=True);return rec
 record(0)
 for it in range(1,100000):
  if time.monotonic()>end:break
  r=op.rhs-op.A(Q)-op.free@x;target=op.free@c-op.A(Z)+r/mu
  count=[0]
  def cb(v):count[0]+=1
  yn,st=cg(op.G,target,x0=y,M=pre,rtol=1e-10,atol=1e-13,maxiter=20,callback=cb)
  y=yn;AY=op.AT(y);Qn=[];Zn=[]
  for q,ay in zip(Q,AY):
   W=q+mu*ay;e,V=linalg.eigh((W+W.T)/2,check_finite=False);qp=(V*np.maximum(e,0))@V.T;zn=(qp-W)/mu
   Qn.append(qp);Zn.append(zn)
  x=x+mu*(op.free.T@y-c);Q,Z=Qn,Zn
  if time.monotonic()-last>15:
   rr=record(it,{'cg_status':int(st),'cg_iterations':count[0]});last=time.monotonic()
   if rr['unverified_width_mHa']<1.3 and rr['primal_l2']<5e-8:break
 record('final')
 ex=export(op,best['Q'],best['x'],out/'export');record_json(out/'discovery.json',{'best':best['record'],'export':ex,'seconds':time.monotonic()-start,'construction':op.info,'status':'requires_independent_exact_replay'});print('EXPORTED',ex,flush=True)
if __name__=='__main__':run()
