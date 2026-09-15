"""Highest-weight spin decomposition of the actual uploaded cubic dictionaries.
Proposal algebra only; all resulting certificates use the unchanged checker.
"""
from common import *
from collections import defaultdict

def normalize_word(w):
    cr=[p for c,p in w if c];an=[p for c,p in w if not c]
    if len(set(cr))!=len(cr) or len(set(an))!=len(an):return None,0
    phase=(-1)**(sum(cr[i]>cr[j] for i in range(len(cr)) for j in range(i+1,len(cr)))+sum(an[i]<an[j] for i in range(len(an)) for j in range(i+1,len(an))))
    return tuple([(1,p) for p in sorted(cr)]+[(0,p) for p in sorted(an,reverse=True)]),phase

def raising(groups):
    look={tuple(map(tuple,w)):(g,j) for g,gg in enumerate(groups) if g>=42 for j,w in enumerate(gg['words'])}
    matrices={};trans={}
    for g in range(42,58):
        by=defaultdict(list)
        for j,w in enumerate(groups[g]['words']):
            for k,(c,p) in enumerate(w):
                if (c==1 and p%2==1) or (c==0 and p%2==0):
                    ww=[list(z) for z in w];ww[k][1]=p-1 if c else p+1
                    word,s=normalize_word(ww)
                    if not s:continue
                    h,i=look[word];by[h].append((i,j,s*(1 if c else -1)))
        if len(by)>1:raise AssertionError(('raising_splits_parity',g,list(by)))
        if not by:continue
        h,rows=next(iter(by.items()));v=sparse.coo_matrix(([r[2] for r in rows],([r[0] for r in rows],[r[1] for r in rows])),shape=(len(groups[h]['words']),len(groups[g]['words'])),dtype=np.int64).tocsr();v.sum_duplicates();v.eliminate_zeros()
        matrices[g]=v;trans[g]=h
    return matrices,trans

def kernel(U):
    gram=(U@U.T).toarray()
    if not np.array_equal(gram,3*np.eye(U.shape[0],dtype=int)):raise AssertionError(('bad_norm',np.unique(gram)))
    if max(np.asarray(U.getnnz(axis=0)).ravel())>1:raise AssertionError('overlapping supports')
    cols=[]
    for j in np.flatnonzero(np.asarray(U.getnnz(axis=0)).ravel()==0):
        v=np.zeros(U.shape[1]);v[j]=1.;cols.append(v)
    for i in range(U.shape[0]):
        js=U.indices[U.indptr[i]:U.indptr[i+1]];s=U.data[U.indptr[i]:U.indptr[i+1]]
        if len(js)!=3:raise AssertionError(('notthree',i,len(js)))
        v=np.zeros(U.shape[1]);v[js[0]]=s[1]/2;v[js[1]]=-s[0]/2;cols.append(v)
        v=np.zeros(U.shape[1]);v[js[0]]=s[0]/3;v[js[1]]=s[1]/3;v[js[2]]=-2*s[2]/3;cols.append(v)
    K=np.array(cols).T
    if np.max(abs(U@K))>1e-15:raise AssertionError('kernel')
    return K

def construct(op,Q):
    groups=op.meta['groups'];R,nextg=raising(groups);pos={m[0]:i for i,m in enumerate(op.members)}
    print('raising',[(g,nextg[g],R[g].shape) for g in R],flush=True)
    high=[g for g in range(42,58) if g not in R];chains=[]
    for h in high:
        chain=[h]
        while any(v==chain[0] for v in nextg.values()):chain.insert(0,next(g for g,v in nextg.items() if v==chain[0]))
        if len(chain)!=4:raise AssertionError(chain)
        chains.append(chain)
    specs=[]
    for i,m in enumerate(op.members):
        if m[0]<42:specs.append((m[0],op.V[i].copy(),Q[i].copy(),{'kind':'unchanged'}))
    for chain in chains:
        low,neg,mid,high=chain;Km=kernel(R[mid]);pinv=Km.T/np.sum(Km*Km,axis=0)[:,None]
        qD=np.zeros((Km.shape[1],)*2);qQ=np.zeros((len(groups[high]['words']),)*2)
        for loc,g in enumerate(chain):
            i=pos[g];D=op.V[i]@Q[i]@op.V[i].T
            T=sparse.eye(len(groups[g]['words']),format='csr');cur=g
            while cur!=high:T=R[cur]@T;cur=nextg[cur]
            qQ+=np.asarray(T@D@T.T)/[36,12,3,1][loc]
            if loc==1:
                t=pinv@R[neg];qD+=t@D@t.T
            elif loc==2:qD+=pinv@D@pinv.T
        specs.append((mid,Km,(qD+qD.T)/2,{'kind':'spin_half_highest','chain':chain,'kernel_entries':int(np.count_nonzero(Km))}))
        specs.append((high,np.eye(len(qQ)),(qQ+qQ.T)/2,{'kind':'spin_three_half_highest','chain':chain}))
    return specs

if __name__=='__main__':
    start=time.monotonic();op=load_operator();Q,_,_=load_source(op);specs=construct(op,Q)
    print('specs',[(g,V.shape,info['kind']) for g,V,q,info in specs]);Aold=op.A(Q);Anew=np.zeros_like(Aold);pos={m[0]:i for i,m in enumerate(op.members)}
    for g,V,q,info in specs:Anew+=op.M[pos[g]]@((V@q@V.T).ravel())
    print('relative_embedding_error',np.linalg.norm(Anew-Aold)/np.linalg.norm(Aold),'max',max(abs(Anew-Aold)),'seconds',time.monotonic()-start)
    record_json(OUT/'spin_irrep_probe.json',{'relative_embedding_error':float(np.linalg.norm(Anew-Aold)/np.linalg.norm(Aold)),'max_embedding_error':float(max(abs(Anew-Aold))),'seconds':time.monotonic()-start,'gram_entries':sum(q.size for _,_,q,_ in specs),'independent_entries':sum(len(q)*(len(q)+1)//2 for _,_,q,_ in specs),'description':'Exact highest-weight basis; floating Gram transport checked against complete uploaded coefficient map. No energy claimed.'})
