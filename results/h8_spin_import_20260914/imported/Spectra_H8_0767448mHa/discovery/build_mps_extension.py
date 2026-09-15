"""Enlarge current spans using explicitly computed three-body MPS moments.
All added directions are proposal information, not assumed ground-state data.
"""
from common import *
from research.sector_quotient_20260914.normal import projected
import argparse,resource
p=argparse.ArgumentParser();p.add_argument('tag');p.add_argument('--big',type=int,default=48);p.add_argument('--small',type=int,default=16);a=p.parse_args()
out=OUT/a.tag;out.mkdir(exist_ok=False);start=time.monotonic();op=load_operator();Q,y,x=load_source(op);guide=np.load(OUT/'physical_moment_matrices.npz');pos=0;records=[]
for k,pair in enumerate(op.meta['pairs']):
 if len(pair['members'])==2 and len(op.V[pos])>64:
  V=op.V[pos];n,r=V.shape;order=np.array(pair['adjoint_order']);C1=guide[f'C_{pos}'];C2=guide[f'C_{pos+1}'][np.ix_(order,order)]
  e1,E1=linalg.eigh(C1,check_finite=False);e2,E2=linalg.eigh(C2,check_finite=False)
  B=linalg.orth(V);added=[];chosen=[];limit=a.big if n>200 else a.small
  for j in range(n):
   for partner,E,e in ((0,E1,e1),(1,E2,e2)):
    v=E[:,j];v=v-B@(B.T@v);v=v-B@(B.T@v);norm=np.linalg.norm(v)
    if norm>.05:
     v/=norm;added.append(v);chosen.append({'partner':partner,'index':j,'eigenvalue':float(e[j]),'missing_norm':float(norm)});B=np.column_stack((B,v))
    if len(added)>=limit:break
   if len(added)>=limit:break
  W=np.column_stack([V]+added);Wp=np.empty_like(W);Wp[order]=W
  for j,Vnew in ((pos,W),(pos+1,Wp)):
   q=np.zeros((Vnew.shape[1],)*2);q[:r,:r]=Q[j];Q[j]=q;op.V[j]=Vnew;op.Q[j]=q;op.identity[j]=False
  records.append({'pair':k,'physical':n,'old':r,'new':len(Q[pos]),'chosen':chosen})
  print('EXTENSION',k,n,r,len(Q[pos]),flush=True)
 pos+=len(pair['members'])
np.savez_compressed(out/'representation.npz',y=y,x=x,**{f'V_{i}':V for i,V in enumerate(op.V)},**{f'Q_{i}':q for i,q in enumerate(Q)})
record_json(out/'representation.json',{'extra_paired_columns':[],'modified_existing_spans':records,'block_members':op.members,'block_pair_ids':op.ids,'gram_entries':sum(q.size for q in Q),'total_with_nonsinglet':sum(q.size for q in Q)+18128,'source':'Hamiltonian + current compact spans + three-body moments of the supplied rational MPS','full_Gram_teacher_used':False,'many_body_states_enumerated':0})
quo=Quotient(op);G=np.zeros((len(op.rhs),)*2,order='F');pos=0

def gram_add(rows,P,weight=1):
 for st in range(0,len(rows),128):
  rr=rows[st:st+128];idx=np.ix_(rr,rows);G[idx]+=weight*(P[st:st+128]@P.T)
def cross_add(rows,cols,P,L,weight):
 for st in range(0,len(rows),128):
  rr=rows[st:st+128];G[np.ix_(rr,cols)]+=weight*(P[st:st+128]@L.T)
selected=np.load(PREP/'selected.npy');sextic=np.array([len(op.meta['rows'][i])==6 for i in selected])
for k,pair in enumerate(op.meta['pairs']):
 t=time.monotonic();V=op.V[pos];M=op.M[pos];active,P=projected(M,V)
 if len(pair['members'])==1:
  gram_add(active,P);pos+=1
 else:
  n=len(V);order=np.array(pair['adjoint_order']);partner=op.M[pos+1][:,(order[None,:]*n+order[:,None]).ravel()]
  paired=M+partner;paired.eliminate_zeros()
  if paired[sextic].nnz:raise AssertionError('Sextic cancellation error')
  low,L=projected(paired,V);gram_add(active,P,2);cross_add(active,low,P,L,-1);cross_add(low,active,L,P,-1);gram_add(low,L);pos+=2
 print('NORMAL',k,'active',len(active),'cols',P.shape[1],'seconds',time.monotonic()-t,'elapsed',time.monotonic()-start,'rss',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)
 del P
rng=np.random.default_rng(211);errs=[]
for _ in range(3):
 v=rng.normal(size=len(op.rhs));ref=op.A(op.AT(v));errs.append(float(np.linalg.norm(G@v-ref)/np.linalg.norm(ref)))
if max(errs)>1e-8:raise AssertionError('Normal operator mismatch')
print('PROJECT',errs,flush=True)
UU=quo.U;UG=UU.T@G;UGU=UG@UU
for st in range(0,len(G),128):
 rr=slice(st,min(st+128,len(G)));Urr=UU[rr];G[rr]-=Urr@UG;G[rr]-=UG[:,rr].T@UU.T;G[rr]+=(Urr@UGU)@UU.T;G[rr]+=Urr@UU.T
for st in range(0,len(G),128):
 en=min(st+128,len(G));G[st:en,st:en]=(G[st:en,st:en]+G[st:en,st:en].T)/2
 for t in range(en,len(G),128):
  u=min(t+128,len(G));block=(G[st:en,t:u]+G[t:u,st:en].T)/2;G[st:en,t:u]=block;G[t:u,st:en]=block.T
G[np.diag_indices(len(G))]+=1e-10
print('CHOLESKY',time.monotonic()-start,flush=True)
C=linalg.cholesky(G,lower=True,overwrite_a=True,check_finite=False);np.save(out/'cholesky.npy',C)
record_json(out/'normal_receipt.json',{'seconds':time.monotonic()-start,'normal_probe_errors':errs,'gram_entries':sum(q.size for q in Q),'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'many_body_states_enumerated':0,'full_Gram_teacher_used':False})
print('FINISHED',time.monotonic()-start,flush=True)
