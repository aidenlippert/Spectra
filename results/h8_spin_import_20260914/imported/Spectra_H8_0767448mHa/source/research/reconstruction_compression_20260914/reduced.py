"""Reduced PSD search with every CAR coefficient equation retained.

The shared coefficient rows are unknown global moments in the dual. Guiding
MPS moments never appear as constraints. Ordinary rational square factors
are exported to the unchanged independent checker.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
import cvxpy as cp
from experiments.marginal_symbolic import decode,encode,mono,product,number_shift,multiplier_basis,add,scale
from experiments.marginal_coefficient import gram_map,dagger
from research.certificate_scaling.adaptive_block_discovery import partition
from research.certificate_scaling.direct_sparse_discovery import sparse_columns
from research.reconstruction_compression_20260914.inputs import OUT,sha,dump

def frame(data):
    m,n=data['modes'],data['particles'];h=decode(data['hamiltonian'],m,4)
    groups,sig,stats=partition(h,m,'mixed',True)
    # Match the successful campaigns' declared conventions: H6 retained the
    # mixed linear words; H8 used the N>=2 linear quotient. Both are generated
    # from problem size, without loading a successful Gram certificate.
    if n>6:
        groups=[{**g,'words':[w for w in g['words'] if len(w)!=1]} if g['name'].startswith('mixed') else g for g in groups]
    groups=[g for g in groups if g['words']]
    rows=[tuple((1,i) for i in a)+tuple((0,i) for i in b) for k in range(4)
          for a in combinations(range(m),k) for b in combinations(range(m),k)
          if a<=b and sig(tuple((1,i) for i in a)+tuple((0,i) for i in b))==sig(())]
    lookup={w:i for i,w in enumerate(rows)}
    basis=[p for p in multiplier_basis(m,max_body=2) if all(sig(w)==sig(()) for w in p)]
    free=sparse_columns([{w:c for w,c in p.items() if w in lookup} for p in
                        [mono(())]+[product(number_shift(m,n),q) for q in basis]],lookup)
    return h,groups,rows,lookup,basis,free,stats

def prepare(case,fixture,state):
    start=time.monotonic();data=json.loads(Path(fixture).read_text());cert=json.loads(Path(state).read_text())
    folder=OUT/'prepared'/case;folder.mkdir(parents=True,exist_ok=False)
    h,groups,rows,lookup,basis,free,stats=frame(data)
    sparse.save_npz(folder/'free.npz',free)
    for i,g in enumerate(groups):
        sparse.save_npz(folder/f'map_{i}.npz',gram_map(g['words'],lookup).tocsc())
    map_end=time.monotonic()
    dump(folder/'frame.json',{'fixture':str(fixture),'state':str(state),'fixture_sha256':sha(fixture),'state_sha256':sha(state),
                            'groups':groups,'rows':rows,'stats':stats,'map_seconds':map_end-start})
    print(json.dumps({'stage':'maps','case':case,'seconds':map_end-start,'rows':len(rows),
          'dimensions':[len(g['words']) for g in groups],'entries':sum(len(g['words'])**2 for g in groups)}),flush=True)
    from research.reconstruction_compression_20260914.moments import Proposal
    oracle=Proposal(data,cert)
    # One trie for all requested entries, rather than a separate norm/action
    # contraction for every column. No many-body vector is materialized.
    Cs={}
    for i,g in enumerate(groups):
        Cs[f'C_{i}']=oracle.matrix(g['words'])
        if i%10==0:print(json.dumps({'stage':'moments','block':i,'cache':len(oracle.cache)}),flush=True)
    np.savez_compressed(folder/'moments.npz',**Cs)
    oracle.fill(h)
    energy=sum(float(c)*oracle.cache[w] for w,c in h.items())
    receipt={'case':case,'map_seconds':map_end-start,'moment_seconds':time.monotonic()-map_end,
             'numerical_energy_Ha':energy,'moment_count':len(oracle.cache),'transfer_steps':oracle.steps,
             'minimum_moment_eigenvalue':min(float(np.linalg.eigvalsh(C)[0]) for C in Cs.values()),
             'many_body_states_enumerated':0,'numeric_moments_accept_nothing':True}
    dump(folder/'receipt.json',receipt);print(json.dumps(receipt),flush=True)

def project(gmap,V):
    """Sparse CAR map composed with a fixed congruence in bounded batches."""
    n,r=V.shape;ii,jj=np.triu_indices(r);out=[]
    for a in range(0,len(ii),8):
        i,j=ii[a:a+8],jj[a:a+8]
        X=np.einsum('ib,jb->ijb',V[:,i],V[:,j])
        X+=np.einsum('ib,jb->ijb',V[:,j],V[:,i])*(i!=j)[None,None,:]
        out.append(sparse.csc_matrix(gmap@X.reshape(n*n,-1)))
    return sparse.hstack(out,format='csc'),ii,jj

def make_maps(groups,Cs,mode,rank,teacher=None):
    answer=[];used=set();diagnostics=[]
    teacher_by_words={frozenset(tuple(tuple(x) for x in w) for w in b['words']):b for b in teacher['blocks']} if teacher else {}
    for i,g in enumerate(groups):
        if i in used:continue
        words=g['words'];n=len(words)
        if max(map(len,words))<=2:
            answer.append(([(i,np.eye(n))],g['name']));continue
        if mode=='teacher':
            old=teacher_by_words.get(frozenset(words))
            if old is None:
                # Old exporters omitted exactly zero small blocks.
                if n>4:raise ValueError('Teacher dictionary mismatch')
                answer.append(([(i,np.eye(n))],g['name']));continue
            loc={tuple(tuple(x) for x in w):j for j,w in enumerate(old['words'])}
            L=np.array(old['factor'],dtype=float)[:,[loc[w] for w in words]]/teacher['denominator']
            ev,U=np.linalg.eigh(L.T@L);V=U[:,np.argsort(ev)[-min(rank,n):]]
            answer.append(([(i,V)],g['name']));continue
        if mode=='unpaired':
            ev,V=np.linalg.eigh(Cs[f'C_{i}']);V=np.rint(V[:,:min(rank,n)]*10**8)/10**8
            answer.append(([(i,V)],g['name']))
            diagnostics.append({'block':i,'lowest':float(ev[0]),'highest_retained':float(ev[min(rank,n)-1])})
            continue
        # Pair an operator and its adjoint with the *same* reduced PSD block.
        # Odd anticommutators cancel degree six exactly, but the matching
        # equations include those rows and exact acceptance checks them again.
        partner=None;order=None
        for j,other in enumerate(groups):
            if j==i or j in used or len(other['words'])!=n:continue
            loc={w:k for k,w in enumerate(other['words'])}
            if all(dagger(w) in loc for w in words):partner=j;order=[loc[dagger(w)] for w in words];break
        if partner is None:raise ValueError('Missing exact adjoint partner')
        C=Cs[f'C_{i}']+Cs[f'C_{partner}'][np.ix_(order,order)]
        if mode=='simple':V=np.eye(n)[:,:min(rank,n)]
        elif mode=='low':
            ev,V=np.linalg.eigh(C);V=V[:,:min(rank,n)]
            diagnostics.append({'block':i,'lowest':float(ev[0]),'highest_retained':float(ev[min(rank,n)-1])})
        elif mode=='cuts':
            V,details=cut_maps(words,C,rank,details=True)
            diagnostics.append({'block':i,**details})
        else:raise ValueError('Unknown map rule')
        # Quantized maps are fixed before optimization; independent replay
        # accepts the final expanded integer factors, not these floats.
        V=np.rint(V*10**8)/10**8
        W=np.empty_like(V);W[order,:]=V
        answer.append(([(i,V),(partner,W)],g['name']+' paired'));used.add(partner)
    return answer,diagnostics

def cut_maps(words,C,rank,details=False):
    """Nested prefix-tree maps at orbital cuts; local covariances come from C.

    At each cut, disjoint prefix classes provide overlapping cluster columns.
    A root congruence mixes these columns at an independent proof rank. Both
    levels are fixed maps on the CAR dictionary, not conditions on the state.
    This is a finite two-level positivity construction, not a complete RDM
    renormalization hierarchy.
    """
    n=len(words);m=max(i for w in words for c,i in w)+1
    columns=[np.eye(n)[:,i] for i in range(n)]
    representatives=list(words);stages=[]
    for cut in sorted(set((max(1,m//3),max(1,2*m//3))),reverse=True):
        classes={}
        for j,w in enumerate(representatives):
            prefix=tuple((c,i) for c,i in w if i<cut)
            classes.setdefault(prefix,[]).append(j)
        next_columns=[];next_representatives=[]
        for prefix,ix in sorted(classes.items()):
            B=np.column_stack([columns[j] for j in ix]);ev,W=np.linalg.eigh(B.T@C@B)
            keep=min(2,len(ev));sub=B@W[:,:keep]
            next_columns.extend(sub[:,j] for j in range(keep));next_representatives.extend([prefix]*keep)
        stages.append({'cut':cut,'prefix_nodes':len(classes),'input_dimension':len(columns),'output_dimension':len(next_columns),'per_node_dimension_cap':2})
        columns=next_columns;representatives=next_representatives
    B=np.column_stack(columns)
    # Disjoint support at each prefix node and orthogonal child maps make B
    # an isometry. Check, rather than silently repair, the nested construction.
    if not np.allclose(B.T@B,np.eye(B.shape[1]),atol=1e-10):raise AssertionError('Cut-map isometry')
    ev,W=np.linalg.eigh(B.T@C@B)
    V=B@W[:,:min(rank,len(ev))]
    return (V,{'cut_stages':stages,'root_input_dimension':B.shape[1],'root_output_dimension':V.shape[1]}) if details else V

def solve(case,mode,rank,seconds=60,teacher_path=None,tag='',eliminate=False):
    started=time.monotonic();prepared=OUT/'prepared'/case;meta=json.loads((prepared/'frame.json').read_text())
    fixture=Path(meta['fixture']);state=Path(meta['state'])
    if sha(fixture)!=meta['fixture_sha256'] or sha(state)!=meta['state_sha256']:raise ValueError('Prepared input changed')
    data=json.loads(fixture.read_text());h,groups,rows,lookup,basis,free,stats=frame(data)
    if json.loads(json.dumps(groups))!=meta['groups'] or json.loads(json.dumps(rows))!=meta['rows']:
        raise ValueError('Prepared coefficient frame is stale')
    Cs=np.load(prepared/'moments.npz')
    teacher=json.loads(Path(teacher_path).read_text()) if mode=='teacher' else None
    if mode=='teacher' and teacher is None:raise ValueError('Explicit teacher required')
    maps,diagnostics=make_maps(groups,Cs,mode,rank,teacher)
    folder=OUT/'candidates'/f'{case}_{mode}_{rank}{tag}';folder.mkdir(parents=True,exist_ok=False)
    original_basis_dimension=len(basis);active_rows=np.arange(len(rows))
    if eliminate:
        if mode not in ('low','cuts','simple'):raise ValueError('Exact cancellation is required for elimination')
        from research.reconstruction_compression_20260914.ideal_rank import certify
        proof=certify(data['modes'],data['particles'],basis);dump(folder/'ideal_elimination.json',proof)
        keep=[i for i,p in enumerate(basis) if max(map(len,p),default=0)<=2]
        basis=[basis[i] for i in keep];free=free[:,[0]+[i+1 for i in keep]]
        active_rows=np.array([i for i,w in enumerate(rows) if len(w)<=4]);free=free[active_rows,:]
    x=cp.Variable(free.shape[1]);expr=free@x;variables=[];projected=[];packed=[]
    for members,name in maps:
        r=members[0][1].shape[1];P=None
        if len(members)==2:
            i,V=members[0];j,W=members[1];n=V.shape[0]
            loc={w:k for k,w in enumerate(groups[j]['words'])}
            order=np.array([loc[dagger(w)] for w in groups[i]['words']])
            if not np.array_equal(W[order],V):raise ValueError('Adjoint congruence mismatch')
            gmap=sparse.load_npz(prepared/f'map_{i}.npz')
            # Transpose the adjoint block's pair indices. S is symmetric, so
            # this is the same quadratic form, with each odd anticommutator
            # cancellation now visible before the symmetric projection.
            other=sparse.load_npz(prepared/f'map_{j}.npz')[:,(order[None,:]*n+order[:,None]).ravel()]
            combined=gmap+other;combined.eliminate_zeros()
            # Integer CAR coefficients cancel before dense floating projection.
            if combined[[k for k,w in enumerate(rows) if len(w)==6],:].nnz:
                raise AssertionError('Odd anticommutator left sextic coefficients')
            P,ii,jj=project(combined[active_rows,:],V)
        else:
            i,V=members[0];gmap=sparse.load_npz(prepared/f'map_{i}.npz');P,ii,jj=project(gmap[active_rows,:],V)
        P.eliminate_zeros();Q=cp.Variable((r,r),PSD=True)
        expr+=P@Q[ii,jj];variables.append(Q);projected.append(P);packed.append((ii,jj))
    rhs=np.array([float(h.get(rows[i],0)) for i in active_rows])
    rownorm=np.asarray(free.power(2).sum(axis=1)).ravel()
    for P in projected:rownorm+=np.asarray(P.power(2).sum(axis=1)).ravel()
    row_scale=1/np.maximum(1.,np.sqrt(rownorm))
    equality=cp.multiply(row_scale,expr-rhs)==0
    problem=cp.Problem(cp.Maximize(x[0]),[equality])
    build=time.monotonic()-started
    receipt={'case':case,'mode':mode,'rank_cap':rank,'teacher_used':bool(teacher),'all_CAR_rows':len(rows),
             'degree_six_rows':sum(len(w)==6 for w in rows),'sector_multiplier_dimension':len(basis),
             'full_gram_entries':sum(len(g['words'])**2 for g in groups),'reduced_gram_entries':sum(q.shape[0]**2 for q in variables),
             'active_equations_after_exact_elimination':len(active_rows),'original_sector_multiplier_dimension':original_basis_dimension,
             'exact_quartic_ideal_elimination':eliminate,
             'projected_map_nonzeros':sum(p.nnz for p in projected),'map_build_seconds':build,'diagnostics':diagnostics,
             'independent_preparation_receipt':str(prepared/'receipt.json'),'numerical_equality_tolerance':1e-8,
             'unknown_state_restricted_to_MPS':False,'many_body_states_enumerated':0}
    dump(folder/'construction.json',receipt);print(json.dumps({'stage':'solve',**{k:v for k,v in receipt.items() if k!='diagnostics'}}),flush=True)
    t=time.monotonic();problem.solve(solver='SCS',eps=1e-8,max_iters=30000,time_limit_secs=seconds,verbose=False,use_indirect=False)
    receipt.update(solve_seconds=time.monotonic()-t,solver_status=problem.status,solver_iterations=problem.solver_stats.num_iters)
    if x.value is None or any(q.value is None for q in variables):
        receipt.update(status='no_numerical_candidate',total_seconds=time.monotonic()-started);dump(folder/'discovery.json',receipt);print(json.dumps(receipt));return
    xv=np.asarray(x.value);coeff=np.asarray(free@xv);blocks=[];min_ev=0.;actions=[]
    for (members,name),Q,P,(ii,jj) in zip(maps,variables,projected,packed):
        q=(Q.value+Q.value.T)/2;coeff+=P@q[ii,jj]
        ev,W=np.linalg.eigh(q);min_ev=min(min_ev,float(ev[0]));R=np.sqrt(np.maximum(ev,0))[:,None]*W.T
        C=sum(V.T@Cs[f'C_{i}']@V for i,V in members)
        actions.append({'name':name,'dimension':len(q),'numeric_trace':float(np.trace(q)),
                        'numeric_state_action_Ha':float(np.trace(C@q))})
        common=None
        for i,V in members:
            if common is None:
                Z=np.rint((R@V.T)*10**9).astype(np.int64);common=Z
            else:
                first=members[0][0];loc={w:k for k,w in enumerate(groups[i]['words'])}
                order=[loc[dagger(w)] for w in groups[first]['words']]
                Z=np.empty_like(common);Z[:,order]=common
            Z=Z[np.any(Z,axis=1)]
            if len(Z):blocks.append({'name':groups[i]['name'],'words':groups[i]['words'],'factor':Z.tolist()})
    X=add(*(scale(p,F(round(float(v)*10**10),10**10)) for p,v in zip(basis,xv[1:])))
    cert={'modes':data['modes'],'particles':data['particles'],'hamiltonian':encode(h),'operator_degree':3,
          'b':str(F(round(float(xv[0])*10**12),10**12)),'number_multiplier':encode(X),'denominator':10**9,'blocks':blocks}
    (folder/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    weights=np.array([1 if tuple(i for c,i in rows[k] if c)==tuple(i for c,i in rows[k] if not c) else 2 for k in active_rows])
    receipt.update(status='proposal_requires_exact_replay',proposed_b=float(xv[0]),minimum_numeric_PSD_eigenvalue=min_ev,
                   full_numeric_coefficient_l1=float(weights@abs(coeff-rhs)),factor_rows=sum(len(b['factor']) for b in blocks),
                   certificate_bytes=(folder/'certificate.json').stat().st_size,total_seconds=time.monotonic()-started)
    receipt['gram_state_actions']=actions
    dump(folder/'discovery.json',receipt);print(json.dumps({k:v for k,v in receipt.items() if k not in ('diagnostics','gram_state_actions')}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='cmd',required=True)
    a=s.add_parser('prepare');a.add_argument('case');a.add_argument('fixture');a.add_argument('state')
    a=s.add_parser('solve');a.add_argument('case');a.add_argument('mode',choices=['teacher','simple','low','cuts','unpaired']);a.add_argument('rank',type=int);a.add_argument('--seconds',type=int,default=60);a.add_argument('--teacher');a.add_argument('--tag',default='');a.add_argument('--eliminate',action='store_true')
    a=p.parse_args()
    if a.cmd=='prepare':prepare(a.case,a.fixture,a.state)
    else:solve(a.case,a.mode,a.rank,a.seconds,a.teacher,a.tag,a.eliminate)
