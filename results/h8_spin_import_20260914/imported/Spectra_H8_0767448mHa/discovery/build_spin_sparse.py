from spin_irrep import *
from scipy.sparse.linalg import splu
import resource,gc
out=OUT/'spin_sparse';out.mkdir(exist_ok=False);start=time.monotonic()
op=load_operator();Q,y,x=load_source(op);specs=construct(op,Q);pos={v[0]:i for i,v in enumerate(op.members)}
S=[];meta=[];qs=[];vs=[];rhsold=op.A(Q)
for j,(g,V,q,info) in enumerate(specs):
    n=V.shape[1];K=sparse.csc_matrix(V);ss=(op.M[pos[g]]@sparse.kron(K,K,format='csc')).tocsc()
    trans=(np.arange(n)[None,:]*n+np.arange(n)[:,None]).ravel();ss=(ss+ss[:,trans])*.5;ss.eliminate_zeros();S.append(ss)
    sparse.save_npz(out/f'map_{j}.npz',ss);qs.append(q);vs.append(V);meta.append({'physical_group':g,'dimension':n,'nnz':ss.nnz,**info})
print('maps',sum(s.nnz for s in S),'dimensions',[len(q) for q in qs], 'seconds',time.monotonic()-start,flush=True)
reconstructed=sum((s@q.ravel() for s,q in zip(S,qs)),np.zeros(len(op.rhs)))
assert np.max(abs(rhsold-reconstructed))<1e-11
np.savez_compressed(out/'seed.npz',y=y,x=x,**{f'Q_{i}':q for i,q in enumerate(qs)})
np.savez_compressed(out/'representation.npz',**{f'V_{i}':v for i,v in enumerate(vs)})
record_json(out/'representation.json',{'blocks':meta,'gram_entries':sum(q.size for q in qs),'total_including_nonsinglet':sum(q.size for q in qs)+18128,'coefficient_nonzeros':sum(s.nnz for s in S),'many_body_states_enumerated':0,'full_cubic_teacher_used':False,'warm_seed':'uploaded compact best candidate','max_embedding_defect':float(np.max(abs(rhsold-reconstructed)))})
G=(op.free@op.free.T).tocsc()
for j,s in enumerate(S):
    G+=s@s.T
    if j%8==0:print('normal',j,G.nnz,'MiB',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,flush=True)
G=(G+G.T)*.5;G.eliminate_zeros();sparse.save_npz(out/'normal.npz',G)
print('NORMAL_READY',G.shape,G.nnz,'seconds',time.monotonic()-start,flush=True)
record_json(out/'build.json',{'seconds':time.monotonic()-start,'normal_nonzeros':G.nnz,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
# The sparse solve is a proposal preconditioner only; its accuracy is measured.
for reg in (1e-10,):
    st=time.monotonic();lu=splu(G+sparse.eye(G.shape[0],format='csc')*reg)
    rng=np.random.default_rng(912);r=G@rng.normal(size=G.shape[0]);z=lu.solve(r)
    rec={'factor_seconds':time.monotonic()-st,'factor_nonzeros':lu.L.nnz+lu.U.nnz,'relative_residual':float(np.linalg.norm(G@z-r)/np.linalg.norm(r)),'regularization':reg,'total_seconds':time.monotonic()-start,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};print('FACTOR',json.dumps(rec),flush=True);record_json(out/'factor_probe.json',rec)
