"""Primal-dual numerical discovery on the uploaded molecular coefficient maps.
No result is accepted here: export is checked in a separate standard-library process.
"""
from common import *
import argparse,resource
from scipy.sparse.linalg import cg
from research.sector_quotient_20260914.normal import projected

p=argparse.ArgumentParser();p.add_argument('tag');p.add_argument('--seconds',type=float,default=500)
p.add_argument('--add',type=int,default=0);p.add_argument('--restart');p.add_argument('--method',choices=['alm','admm'],default='alm')
p.add_argument('--physical-dual',action='store_true');p.add_argument('--reuse-representation');p.add_argument('--mu',type=float,default=.5);p.add_argument('--pricing',choices=['paired','individual'],default='paired')
a=p.parse_args(); out=OUT/a.tag;out.mkdir(exist_ok=False)
start=time.monotonic();op=load_operator();Q,y,x=load_source(op);quo=Quotient(op)
# Maps are unchanged except for added positive generators. The baseline has 58 blocks.
base_count=len(Q);extras=[];chosen_cache=CACHE
if a.restart:
    raw=np.load(a.restart);Q=[raw[f'Q_{i}'].copy() for i in range(base_count)];y=raw['y'].copy()
y=quo.project(y)
if a.reuse_representation:
    previous=Path(a.reuse_representation)
    rp=json.loads((previous/'representation.json').read_text());vz=np.load(previous/'representation.npz')
    extras=rp['extra_paired_columns']
    if rp.get('modified_existing_spans'):
        for i in range(base_count):
            vnew=vz[f'V_{i}'];nnew=vnew.shape[1]
            qnew=np.zeros((nnew,nnew));old=Q[i];qnew[:old.shape[0],:old.shape[1]]=old
            Q[i]=qnew;op.Q[i]=qnew;op.V[i]=vnew
            op.identity[i]=nnew==len(vnew) and np.array_equal(vnew,np.eye(nnew))
        chosen_cache=previous
    for info in extras:
        k=info['pair'];pos=info['old_block'];pair=op.meta['pairs'][k];n=op.V[pos].shape[0];order=np.array(pair['adjoint_order'])
        M=op.M[pos]+op.M[pos+1][:,(order[None,:]*n+order[:,None]).ravel()]
        M=(M+M[:,(np.arange(n)[None,:]*n+np.arange(n)[:,None]).ravel()])*.5;M.eliminate_zeros()
        j=len(Q);V=vz[f'V_{j}'];zero=np.zeros((V.shape[1],)*2)
        op.M.append(M);op.V.append(V);op.Q.append(zero);op.members.append(pair['members']);op.ids.append(k);op.identity.append(False)
        Q.append(raw[f'Q_{j}'].copy() if a.restart else zero)
if a.add:
    yphys=quo.y0-y; pos=0
    for k,pair in enumerate(op.meta['pairs']):
        n=op.V[pos].shape[0]
        if len(pair['members'])==2 and n>64:
            order=np.array(pair['adjoint_order'])
            M=op.M[pos]+op.M[pos+1][:,(order[None,:]*n+order[:,None]).ravel()]
            M=(M+M[:,(np.arange(n)[None,:]*n+np.arange(n)[:,None]).ravel()])*.5
            M.eliminate_zeros()
            C=np.asarray(M.T@yphys).reshape(n,n);C=(C+C.T)/2
            ev,V=linalg.eigh(C,check_finite=False)
            take=np.flatnonzero(ev < -1e-7)[:a.add]
            if len(take):
                V=V[:,take];zero=np.zeros((len(take),len(take)))
                op.M.append(M);op.V.append(V);op.Q.append(zero);op.members.append(pair['members']);op.ids.append(k);op.identity.append(False)
                Q.append(zero)
                info={'pair':k,'rank':len(take),'minimum_pricing_eigenvalue':float(ev[0]),'included_eigenvalues':ev[take].tolist(),'old_block':pos}
                extras.append(info)
                print('ADDED',info,flush=True)
        pos+=len(pair['members'])
np.savez_compressed(out/'representation.npz',**{f'V_{i}':v for i,v in enumerate(op.V)})
record_json(out/'representation.json',{'source':'uploaded compact spans','extra_paired_columns':extras,'modified_existing_spans':rp.get('modified_existing_spans',[]) if a.reuse_representation else [],'block_members':op.members,'block_pair_ids':op.ids,'full_Gram_teacher_used':False,'many_body_states_enumerated':0,'gram_entries':sum(q.size for q in Q),'total_with_nonsinglet':sum(q.size for q in Q)+18128})
# Build an exact (floating proposal) Woodbury correction to the supplied metric.
Cchol=np.load(chosen_cache/'cholesky.npy')
baseinv=lambda b:linalg.cho_solve((Cchol,True),b,check_finite=False)
if extras:
    sel=np.load(PREP/'selected.npy');active_all=np.array([j for j,i in enumerate(sel) if len(op.meta['rows'][i])<=4])
    lookup={int(j):i for i,j in enumerate(active_all)};G=np.zeros((len(active_all),len(active_all)))
    for i in range(base_count,len(Q)):
        active,P=projected(op.M[i],op.V[i])
        if any(int(j) not in lookup for j in active):raise AssertionError('Paired added map retains higher-degree rows')
        inds=np.array([lookup[int(j)] for j in active]);G[np.ix_(inds,inds)]+=P@P.T
    ev,V=linalg.eigh(G,check_finite=False);keep=ev>1e-13*max(ev[-1],1.)
    W=np.zeros((len(op.rhs),int(np.sum(keep))));W[active_all]=V[:,keep]*np.sqrt(ev[keep])
    W-=quo.U@(quo.U.T@W)
    BW=baseinv(W);S=np.eye(W.shape[1])+W.T@BW
    Sch=linalg.cho_factor((S+S.T)/2,lower=True,check_finite=False)
    def inverse(b):
        z=baseinv(b)
        return z-BW@linalg.cho_solve(Sch,W.T@z,check_finite=False)
    print('Woodbury_rank',W.shape[1],'elapsed',time.monotonic()-start,flush=True)
else: inverse=baseinv
gv=inverse(quo.v);den=float(quo.v@gv)
def minverse(b):
    b=quo.project(b);z=inverse(b)-gv*(gv@b)/den
    return quo.project(z)
normal=LinearOperator((len(op.rhs),)*2,matvec=lambda v:quo.A(quo.AT(v)))
pre=LinearOperator(normal.shape,matvec=minverse)
# Normal check is a numerical construction check, not certificate acceptance.
rng=np.random.default_rng(813);v=quo.project(rng.normal(size=len(op.rhs)));v=normal@v
z=minverse(v);err=np.linalg.norm(normal@z-v)/np.linalg.norm(v)
print('normal_preconditioned_residual',err,flush=True)
print('setup_seconds',time.monotonic()-start,'blocks',len(Q),'entries',sum(q.size for q in Q),flush=True)
C=op.AT(quo.y0);rhs=quo.rhs;mu=a.mu
if a.physical_dual:
    y=quo.project(quo.y0-np.load(OUT/'physical_moment_dual.npy'))
    print('INITIALIZED_FROM_ACTUAL_THREE_BODY_MPS_MOMENTS',flush=True)
residual_predictor=FastResidual(op)
history=[];best={'score':-float('inf')};last_log=0.;deadline=time.monotonic()+a.seconds

def log(Q,y,outer,inner,extra=None):
    x=quo.recover(Q);res=(op.A(Q)+op.free@x-op.rhs)/op.scale
    r1=float(np.sum(abs(res)));dual=op.AT(quo.y0-y)
    dmin=min(float(linalg.eigvalsh(c,subset_by_index=(0,0),check_finite=False)[0]) for c in dual)
    dobj=float(quo.offset-rhs@y)
    rec={'outer':outer,'inner':inner,'seconds':time.monotonic()-start,'b':float(x[0]),'raw_width_mHa':1000*(U-x[0]),'selected_row_l1':r1,'primal_l2':float(np.linalg.norm(quo.A(Q)-rhs)),'untrusted_dual_b':dobj,'minimum_dual_eigenvalue':dmin,'mu':mu,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    if extra:rec.update(extra)
    correction,body_corrections=residual_predictor.evaluate(res*1.)
    rec['numerical_spectral_correction_Ha']=correction
    rec['unverified_spectral_width_mHa']=1000*(U-x[0]-correction+2.36e-10)
    score=float(x[0])+correction-2.36e-10
    if score>best['score']:
        best.update(score=score,Q=[q.copy() for q in Q],x=x.copy(),record=rec)
        np.savez_compressed(out/'best.npz',x=x,y=y,**{f'Q_{i}':q for i,q in enumerate(Q)})
    history.append(rec);record_json(out/'history.json',history)
    np.savez_compressed(out/'checkpoint.npz',x=x,y=y,**{f'Q_{i}':q for i,q in enumerate(Q)})
    print(json.dumps(rec),flush=True)
    return rec

log(Q,y,0,0)
if a.method=='admm':
    Z=[positive(c-ay) for c,ay in zip(C,quo.AT(y))]
    for it in range(100000):
        if time.monotonic()>deadline:break
        res=rhs-quo.A(Q);ry=quo.A([c-z for c,z in zip(C,Z)])+res/mu
        count=[0]
        def cb(v):count[0]+=1
        y,status=cg(normal,ry,x0=y,rtol=1e-9,atol=1e-13,maxiter=60,M=pre,callback=cb);y=quo.project(y)
        AY=quo.AT(y);Qnew=[];Znew=[]
        for q,c,ay in zip(Q,C,AY):
            S=q+mu*(ay-c);e,v=linalg.eigh((S+S.T)/2,check_finite=False);qp=(v*np.maximum(e,0))@v.T
            zp=(qp-S)/mu;Qnew.append(qp);Znew.append(zp)
        Q,Z=Qnew,Znew
        if time.monotonic()-last_log>15:
            r=log(Q,y,it+1,1,{'cg_status':int(status),'cg_iterations':count[0]});last_log=time.monotonic()
            if r['primal_l2']<1e-10 and r['minimum_dual_eigenvalue']>-1e-8:break
else:
    for outer in range(80):
        if time.monotonic()>deadline:break
        Qold=[q.copy() for q in Q]
        def evaluate(y,need_derivative=False):
            AY=quo.AT(y);P=[];spectra=[];norm2=0.
            for q,ay,c in zip(Qold,AY,C):
                S=q+mu*(ay-c);e,V=linalg.eigh((S+S.T)/2,check_finite=False)
                ep=np.maximum(e,0);p=(V*ep)@V.T;P.append(p);norm2+=ep@ep
                if need_derivative:
                    ee=e[:,None]-e[None,:];num=ep[:,None]-ep[None,:]
                    L=np.divide(num,ee,out=np.zeros_like(ee),where=abs(ee)>1e-14)
                    eq=abs(ee)<=1e-14;L[eq]=np.broadcast_to((e>0)[:,None],ee.shape)[eq]
                    spectra.append((V,L))
            grad=quo.A(P)-rhs
            return float(norm2/(2*mu)-rhs@y),grad,P,spectra
        f,g,Q,sp=evaluate(y,True);gn0=np.linalg.norm(g)
        goal=max(3e-11,min(gn0*.04,2e-6/(outer+1)**1.5))
        for inner in range(20):
            gn=np.linalg.norm(g)
            if gn<goal or time.monotonic()>deadline:break
            regularizer=max(.001,min(.03,np.sqrt(gn)))
            def hess(v):
                B=quo.AT(v);K=[]
                for z,(V,L) in zip(B,sp):
                    zz=V.T@z@V;K.append(V@(L*zz)@V.T+regularizer*z)
                return mu*quo.A(K)
            H=LinearOperator(normal.shape,matvec=hess);Ppre=LinearOperator(normal.shape,matvec=lambda v:minverse(v)/mu)
            count=[0]
            def cb(v):count[0]+=1
            step,status=cg(H,-g,M=Ppre,rtol=min(.2,max(.001,np.sqrt(gn))),atol=1e-13,maxiter=60,callback=cb)
            step=quo.project(step);slope=float(g@step)
            if slope>=0 or not np.all(np.isfinite(step)):
                step=-minverse(g)/mu;slope=float(g@step)
            alpha=1.;accepted=False
            for ls in range(25):
                yn=quo.project(y+alpha*step);fn,gnn,Qn,_=evaluate(yn)
                if fn<=f+1e-4*alpha*slope+1e-13 or np.linalg.norm(gnn)<.5*gn:
                    accepted=True;break
                alpha*=.5
            if not accepted:
                print('LINE_SEARCH_STALLED',outer,inner,gn,flush=True);break
            y=yn;f,g,Q,sp=evaluate(y,True)
            if time.monotonic()-last_log>15:
                log(Q,y,outer,inner,{'gradient':float(np.linalg.norm(g)),'cg_iterations':count[0],'cg_status':int(status),'line_step':alpha});last_log=time.monotonic()
        r=log(Q,y,outer+1,inner,{'gradient':float(np.linalg.norm(g)),'inner_goal':goal})
        change=np.sqrt(sum(np.sum((q-p)**2) for q,p in zip(Q,Qold))) / mu
        if r['primal_l2']<1e-10 and change<1e-7:break
        if change>1e-5:mu=min(mu*2,1000.)
log(Q,y,'final',0)
print('EXPORT',flush=True)
e=export(op,best['Q'],best['x'],out/'export')
record_json(out/'discovery.json',{'best':best['record'],'export':e,'seconds':time.monotonic()-start,'extra_paired_columns':extras,'many_body_states_enumerated':0,'full_Gram_teacher_used':False,'inherited_numerical_start':str(SOURCE),'setup_normal_preconditioner':str(CACHE),'status':'unverified_until_independent_replay','gram_entries':sum(q.size for q in Q)})
print('FINISHED',json.dumps(e),flush=True)
