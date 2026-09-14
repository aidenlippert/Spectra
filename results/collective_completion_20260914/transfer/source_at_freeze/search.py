"""Reconstruction-preserving paired PSD search with optional spin lifting.

All cross terms inside each selected operator span are optimized. The trial
state supplies the initial span; a coefficient dual can enlarge it. Numerical
outputs are proposals; spin_replay or the existing checker accepts them.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
import cvxpy as cp
from experiments.marginal_symbolic import encode,mono,add,scale,product,number_shift,multiplier_basis
from experiments.marginal_coefficient import dagger
from research.certificate_scaling.direct_sparse_discovery import sparse_columns
from research.reconstruction_compression_20260914.reduced import frame,project,make_maps
from research.reconstruction_compression_20260914.inputs import sha,dump
from research.collective_completion_20260914.spin_replay import setup,alpha_shift

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/collective_completion_20260914'
PREVIOUS=ROOT/'results/reconstruction_compression_20260914/prepared'

def load_case(case,prepared=None):
    prepared=Path(prepared) if prepared else PREVIOUS/case
    meta=json.loads((prepared/'frame.json').read_text());fixture=Path(meta['fixture'])
    if sha(fixture)!=meta['fixture_sha256'] or sha(meta['state'])!=meta['state_sha256']:raise ValueError('Changed prepared input')
    data=json.loads(fixture.read_text())
    if meta.get('kind') in ('paired_quartic_v1','linear_response_v1'):
        groups=[{**g,'words':[tuple(tuple(x) for x in w) for w in g['words']]} for g in meta['groups']]
        rows=[tuple(tuple(x) for x in w) for w in meta['rows']];lookup=set(rows)
        basis=[p for p in multiplier_basis(data['modes'],max_body=1) if any(w in lookup for w in p)]
        return data,groups,rows,basis,prepared,meta
    h,groups,rows,lookup,basis,free,stats=frame(data)
    if json.loads(json.dumps(groups))!=meta['groups'] or json.loads(json.dumps(rows))!=meta['rows']:raise ValueError('Stale CAR frame')
    return data,groups,rows,basis,prepared,meta

def build_maps(groups,rows,prepared,maps,compact=False):
    active=np.array([i for i,w in enumerate(rows) if len(w)<=4]);projected=[];original=[]
    for k,(members,name) in enumerate(maps):
        if compact:
            M=sparse.load_npz(prepared/f'paired_map_{k}.npz');V=members[0][1];P,ii,jj=project(M,V)
            projected.append((P,ii,jj));original.append(M);continue
        i,V=members[0];n=V.shape[0];M=sparse.load_npz(prepared/f'map_{i}.npz')
        if len(members)==2:
            j,W=members[1];loc={w:k for k,w in enumerate(groups[j]['words'])}
            order=np.array([loc[dagger(w)] for w in groups[i]['words']])
            if not np.array_equal(W[order],V):raise ValueError('Inconsistent adjoint maps')
            M+=sparse.load_npz(prepared/f'map_{j}.npz')[:,(order[None,:]*n+order[:,None]).ravel()]
        M.eliminate_zeros()
        if M[[k for k,w in enumerate(rows) if len(w)==6],:].nnz:raise AssertionError('Uncancelled sextic row')
        M=M[active,:].tocsc();P,ii,jj=project(M,V)
        projected.append((P,ii,jj));original.append(M)
    return active,projected,original

def enrich(maps,groups,original,y,add_rank):
    result=[];details=[]
    for (members,name),M in zip(maps,original):
        i,V=members[0];n,r=V.shape
        if len(members)==1 or n==r:result.append((members,name));continue
        C=np.asarray(M.T@y).reshape(n,n);C=(C+C.T)/2
        ev,W=np.linalg.eigh(C)
        P=np.eye(n)-V@np.linalg.pinv(V)
        B=P@W[:,np.where(ev < -1e-7)[0][:add_rank]]
        if not B.shape[1]:result.append((members,name));continue
        u,s,_=np.linalg.svd(B,full_matrices=False)
        fresh=u[:,s>1e-7];V2=np.column_stack((V,fresh));V2=np.rint(V2*1e8)/1e8
        j,_=members[1];loc={w:k for k,w in enumerate(groups[j]['words'])};order=[loc[dagger(w)] for w in groups[i]['words']]
        W2=np.empty_like(V2);W2[order]=V2;result.append(([(i,V2),(j,W2)],name))
        details.append({'name':name,'lowest_full_dual_eigenvalue':float(ev[0]),'old_rank':r,'new_rank':V2.shape[1]})
    return result,details

def solve(case,rank,seconds=180,spin=False,tag='',rounds=1,add_rank=8,prepared=None,solver='SCS',seed=None,enrich_seed=False):
    total=time.monotonic();data,groups,rows,allbasis,prepared,meta=load_case(case,prepared)
    h,hs,delta=setup(data);target=hs if spin else h;m,n=data['modes'],data['particles']
    basis=[p for p in allbasis if max(map(len,p),default=0)<=2]
    active=[i for i,w in enumerate(rows) if len(w)<=4];smallrows=[rows[i] for i in active];lookup={w:i for i,w in enumerate(smallrows)}
    polynomials=[mono(())]+[product(number_shift(m,n),p) for p in basis]
    if spin:polynomials += [product(alpha_shift(m,n),p) for p in basis]
    free=sparse_columns([{w:c for w,c in p.items() if w in lookup} for p in polynomials],lookup)
    rhs=np.array([float(target.get(w,0)) for w in smallrows]);Cs=np.load(prepared/'moments.npz')
    response=meta.get('kind')=='linear_response_v1'
    compact=meta.get('kind') in ('paired_quartic_v1','linear_response_v1')
    if response:
        from research.collective_completion_20260914.linear_response import select
        maps,guide=select(meta,Cs,rank)
    elif compact:
        from research.collective_completion_20260914.prepare import select
        maps,guide=select(meta,Cs,rank)
    else:maps,guide=make_maps(groups,Cs,'low',rank)
    base=OUT/'candidates'/f'{case}_{"spin" if spin else "number"}_r{rank}{tag}';base.mkdir(parents=True,exist_ok=False)
    if seed:
        if response:raise ValueError('Response seed maps require a matching structured exporter')
        raw=np.load(seed);seed_maps=[]
        for k,(members,name) in enumerate(maps):
            i=members[0][0];V=raw[f'V_{k}'];mem=[(i,V)]
            if V.shape[0]!=len(groups[i]['words']):raise ValueError('Seed dictionary mismatch')
            if len(members)==2:
                j=members[1][0];loc={w:k for k,w in enumerate(groups[j]['words'])};order=[loc[dagger(w)] for w in groups[i]['words']]
                W=np.empty_like(V);W[order]=V;mem.append((j,W))
            seed_maps.append((mem,name))
        maps=seed_maps
        if enrich_seed:
            _,_,original=build_maps(groups,rows,prepared,maps,compact)
            maps,guide=enrich(maps,groups,original,raw['dual'],add_rank)
        dump(base/'seed.json',{'source':str(seed),'sha256':sha(seed),'enriched':enrich_seed,'add_rank':add_rank,'guide':guide})
    for rnd in range(rounds):
        start=time.monotonic();folder=base/f'round_{rnd}';folder.mkdir()
        active,projected,original=build_maps(groups,rows,prepared,maps,compact)
        x=cp.Variable(free.shape[1]);expr=free@x;variables=[]
        for members,(P,ii,jj) in zip(maps,projected):
            r=members[0][0][1].shape[1];q=cp.Variable((r,r),PSD=True);variables.append(q);expr+=P@q[ii,jj]
        norm=np.asarray(free.power(2).sum(axis=1)).ravel()
        for P,ii,jj in projected:norm+=np.asarray(P.power(2).sum(axis=1)).ravel()
        row_scale=1/np.maximum(1,np.sqrt(norm));eq=cp.multiply(row_scale,expr-rhs)==0
        problem=cp.Problem(cp.Maximize(x[0]),[eq]);build_seconds=time.monotonic()-start
        details={'case':case,'round':rnd,'spin_lifting':spin,'rank_cap_initial':rank,'gram_entries':sum(q.shape[0]**2 for q in variables),
            'gram_dimensions':[q.shape[0] for q in variables],'full_gram_entries':sum(len(g['words'])**2 for g in groups),
            'coefficient_equations':len(smallrows),'free_variables':free.shape[1],
            'projected_map_nonzeros':sum(P.nnz for P,_,_ in projected),'build_seconds':build_seconds,
            'prepared_dependency':str(prepared),'fixture':meta['fixture'],'state':meta['state'],
            'full_Gram_teacher_used':False,'many_body_states_enumerated':0,'guide':guide}
        dump(folder/'construction.json',details);print(json.dumps({'stage':'solve',**{k:v for k,v in details.items() if k!='guide'}}),flush=True)
        t=time.monotonic()
        options={'eps':1e-8,'max_iters':100000,'time_limit_secs':seconds,'use_indirect':False} if solver=='SCS' else {'tol_gap_abs':1e-9,'tol_feas':1e-9,'tol_gap_rel':1e-9,'max_iter':300,'time_limit':seconds}
        problem.solve(solver=solver,verbose=False,**options)
        details.update(solve_seconds=time.monotonic()-t,solver=solver,solver_status=problem.status,solver_iterations=problem.solver_stats.num_iters)
        if x.value is None or any(q.value is None for q in variables):
            details['status']='no_candidate';dump(folder/'discovery.json',details);break
        xv=np.array(x.value);y=np.array(eq.dual_value)*row_scale
        np.savez_compressed(folder/'raw.npz',x=xv,dual=y,**{f'Q_{i}':q.value for i,q in enumerate(variables)},**{f'V_{i}':mem[0][1] for i,(mem,name) in enumerate(maps)})
        coeff=np.asarray(free@xv);blocks=[];min_ev=0.;actions=[]
        for k,((members,name),Q,(P,ii,jj)) in enumerate(zip(maps,variables,projected)):
            q=(Q.value+Q.value.T)/2;coeff+=P@q[ii,jj];ev,W=np.linalg.eigh(q);min_ev=min(min_ev,float(ev[0]))
            R=np.sqrt(np.maximum(ev,0))[:,None]*W.T
            C=members[0][1].T@Cs[f'C_{k}']@members[0][1] if compact else sum(V.T@Cs[f'C_{i}']@V for i,V in members)
            actions.append(float(np.trace(q@C)))
            if response:
                J=members[0][1];abstract=np.rint((R@J.T)*1e9).astype(np.int64);pair=meta['pairs'][k]
                for i,slots in zip(pair['members'],pair['slot_indices']):
                    Z=abstract[:,slots];Z=Z[np.any(Z,axis=1)]
                    if len(Z):blocks.append({'name':groups[i]['name'],'words':groups[i]['words'],'factor':Z.tolist()})
                continue
            common=None
            for i,V in members:
                if common is None:Z=np.rint((R@V.T)*1e9).astype(np.int64);common=Z
                else:
                    loc={w:k for k,w in enumerate(groups[i]['words'])};order=[loc[dagger(w)] for w in groups[members[0][0]]['words']]
                    Z=np.empty_like(common);Z[:,order]=common
                Z=Z[np.any(Z,axis=1)]
                if len(Z):blocks.append({'name':groups[i]['name'],'words':groups[i]['words'],'factor':Z.tolist()})
        def poly(values):return add(*(scale(p,F(round(float(v)*1e10),10**10)) for p,v in zip(basis,values)))
        X=poly(xv[1:1+len(basis)]);Y=poly(xv[1+len(basis):]) if spin else {}
        core={'modes':m,'particles':n,'operator_degree':3,'hamiltonian':encode(add(target,scale(product(alpha_shift(m,n),Y),-1))) if spin else encode(h),
              'b':str(F(round(float(xv[0])*1e12),10**12)),'number_multiplier':encode(X),'denominator':10**9,'blocks':blocks}
        cert={'kind':'balanced_spin_sos_v1','modes':m,'particles':n,'hamiltonian':encode(h),'alpha_multiplier':encode(Y),'core':core} if spin else core
        (folder/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
        weights=np.array([1 if tuple(i for c,i in w if c)==tuple(i for c,i in w if not c) else 2 for w in smallrows])
        details.update(status='requires_independent_exact_replay',proposed_b=float(xv[0]),numeric_coefficient_l1=float(weights@abs(coeff-rhs)),
            minimum_PSD_eigenvalue=min_ev,total_state_action_Ha=sum(actions),factor_rows=sum(len(b['factor']) for b in blocks),
            certificate_bytes=(folder/'certificate.json').stat().st_size,round_seconds=time.monotonic()-start,total_seconds=time.monotonic()-total)
        dump(folder/'discovery.json',details);print(json.dumps({k:v for k,v in details.items() if k!='guide'}),flush=True)
        if rnd+1<rounds:
            maps,guide=enrich(maps,groups,original,y,add_rank);dump(folder/'enrichment.json',guide)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('rank',type=int);p.add_argument('--seconds',type=int,default=180);p.add_argument('--spin',action='store_true');p.add_argument('--tag',default='');p.add_argument('--rounds',type=int,default=1);p.add_argument('--add-rank',type=int,default=8);p.add_argument('--prepared');p.add_argument('--solver',choices=['SCS','CLARABEL'],default='SCS');p.add_argument('--seed');p.add_argument('--enrich-seed',action='store_true');a=p.parse_args()
    solve(a.case,a.rank,a.seconds,a.spin,a.tag,a.rounds,a.add_rank,a.prepared,a.solver,a.seed,a.enrich_seed)
