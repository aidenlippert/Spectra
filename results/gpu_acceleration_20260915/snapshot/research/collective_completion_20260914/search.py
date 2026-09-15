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
from experiments.marginal_symbolic import encode,mono,add,scale,product,number_shift,multiplier_basis,canonical,adj
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

def build_maps(groups,rows,prepared,maps,compact=False,map_ids=None):
    active=np.array([i for i,w in enumerate(rows) if len(w)<=4]);projected=[];original=[];map_cache={}
    for k,(members,name) in enumerate(maps):
        if compact:
            source_id=map_ids[k] if map_ids is not None else k
            if source_id not in map_cache:map_cache[source_id]=sparse.load_npz(prepared/f'paired_map_{source_id}.npz')
            M=map_cache[source_id];V=members[0][1];P,ii,jj=project(M,V)
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
        if n==r:result.append((members,name));continue
        C=np.asarray(M.T@y).reshape(n,n);C=(C+C.T)/2
        ev,W=np.linalg.eigh(C)
        P=np.eye(n)-V@np.linalg.pinv(V)
        B=P@W[:,np.where(ev < -1e-7)[0][:add_rank]]
        if not B.shape[1]:result.append((members,name));continue
        u,s,_=np.linalg.svd(B,full_matrices=False)
        fresh=u[:,s>1e-7];V2=np.column_stack((V,fresh));V2=np.rint(V2*1e8)/1e8
        mem=[(i,V2)]
        if len(members)==2:
            j,_=members[1];loc={w:k for k,w in enumerate(groups[j]['words'])};order=[loc[dagger(w)] for w in groups[i]['words']]
            W2=np.empty_like(V2);W2[order]=V2;mem.append((j,W2))
        result.append((mem,name))
        details.append({'name':name,'lowest_full_dual_eigenvalue':float(ev[0]),'old_rank':r,'new_rank':V2.shape[1]})
    return result,details

def solve(case,rank,seconds=180,spin=False,tag='',rounds=1,add_rank=8,prepared=None,solver='SCS',seed=None,enrich_seed=False,mag=0,singlet=False,ladder=False,seed_eigen_cutoff=0.,threads=1,spin_twirl=False,share_spin=False,moment_dual=False,verbose=False,atoms=False):
    if (OUT/'manifest.json').exists():raise RuntimeError('Sealed campaign: set a fresh output root before new discovery')
    total=time.monotonic();data,groups,rows,allbasis,prepared,meta=load_case(case,prepared)
    h,hs,delta=setup(data);target=hs if spin else h;m,n=data['modes'],data['particles']
    if (mag or singlet) and not spin:raise ValueError('Spin-sector proof requires exact spin setup')
    if mag not in (0,1) or (singlet and mag):raise ValueError('Unsupported spin sector')
    if ladder and not singlet:raise ValueError('Spin-ladder relations require singlet sector')
    if spin_twirl and not singlet:raise ValueError('Spin averaging currently requires the screened singlet construction')
    Zspin=add(alpha_shift(m,n),mono((),-mag))
    if singlet:
        from research.collective_completion_20260914.spin_screen import spin_squared
        S2=spin_squared(m)
    else:S2={}
    basis=[p for p in allbasis if max(map(len,p),default=0)<=2]
    active=[i for i,w in enumerate(rows) if len(w)<=4];smallrows=[rows[i] for i in active];lookup={w:i for i,w in enumerate(smallrows)}
    polynomials=[mono(())]+[product(number_shift(m,n),p) for p in basis]
    if spin:polynomials += [product(Zspin,p) for p in basis]
    if singlet:polynomials += [S2]
    ladder_basis=[]
    if ladder:
        from research.collective_completion_20260914.spin_screen import ladder_ideal
        for i in range(1,m,2):
            for j in range(0,m,2):
                p=mono(((1,i),(0,j)));q=ladder_ideal(m,p)
                if all(w in lookup or all(v in lookup for v in canonical(adj(mono(w)))) for w in q):
                    ladder_basis.append(p);polynomials.append(q)
    free=sparse_columns([{w:c for w,c in p.items() if w in lookup} for p in polynomials],lookup)
    rhs=np.array([float(target.get(w,0)) for w in smallrows]);Cs=np.load(prepared/'moments.npz')
    fullfree=free;fullrhs=rhs.copy();T=None;spin_row_receipt=None
    if spin_twirl:
        from research.collective_completion_20260914.spin_rows import build as spin_rows
        T,selected,spin_row_receipt=spin_rows(smallrows);Ts=T[selected,:]
        free=(Ts@free).tocsc();rhs=np.asarray(Ts@rhs)
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
        raw=np.load(seed);seed_maps=[]
        seed_layout=json.loads((Path(seed).parent/'construction.json').read_text())
        if seed_layout.get('coherent_atoms') or seed_layout.get('shared_spin_orbits'):
            raise ValueError('This seed loader requires the original unshared Gram-map layout')
        if response:
            from research.collective_completion_20260914.linear_response import embed_paired_seed
            seed_details=json.loads((Path(seed).parent/'construction.json').read_text())
            seed_meta=json.loads((Path(seed_details['prepared_dependency'])/'frame.json').read_text())
            if seed_meta.get('kind')!='paired_quartic_v1' or seed_details.get('shared_spin_orbits'):raise ValueError('Response seeding requires unshared paired quartic maps')
            if seed_meta['rows']!=meta['rows'] or seed_meta['fixture_sha256']!=meta['fixture_sha256']:raise ValueError('Response seed does not match this coefficient frame')
            if seed_eigen_cutoff:raise ValueError('Response seed inclusion does not truncate the old span')
        for k,(members,name) in enumerate(maps):
            i=members[0][0];V=raw[f'V_{k}'];mem=[(i,V)]
            if response:
                old_id=seed_meta['pairs'][k]['members'][0]
                old_words=[tuple(tuple(x) for x in w) for w in seed_meta['groups'][old_id]['words']]
                V=embed_paired_seed(old_words,V,groups[i]['words'],meta['pairs'][k]);seed_maps.append(([(i,V)],name));continue
            if V.shape[0]!=len(groups[i]['words']):
                if meta.get('global_triples') and groups[i]['name'].startswith('global-triples'):
                    seed_maps.append((members,name));continue
                raise ValueError('Seed dictionary mismatch')
            if seed_eigen_cutoff>0 and len(members)==2 and V.shape[0]>32:
                q=raw[f'Q_{k}'];ev,W=np.linalg.eigh((q+q.T)/2);keep=max(2,int((ev>seed_eigen_cutoff).sum()))
                V=np.rint((V@W[:,-keep:])*1e8)/1e8;mem=[(i,V)]
            if len(members)==2:
                j=members[1][0];loc={w:k for k,w in enumerate(groups[j]['words'])};order=[loc[dagger(w)] for w in groups[i]['words']]
                W=np.empty_like(V);W[order]=V;mem.append((j,W))
            seed_maps.append((mem,name))
        maps=seed_maps
        if enrich_seed and not atoms:
            _,_,original=build_maps(groups,rows,prepared,maps,compact)
            maps,guide=enrich(maps,groups,original,raw['dual'],add_rank)
        dump(base/'seed.json',{'source':str(seed),'sha256':sha(seed),'enriched':enrich_seed,'add_rank':add_rank,'seed_eigen_cutoff':seed_eigen_cutoff,'guide':guide,
             'truncation_is_untrusted_proposal':True,'retained_certificate_not_overwritten':True})
    map_ids=list(range(len(maps)))
    if atoms:
        if response or not compact or share_spin:raise ValueError('Atoms require paired quartic maps')
        from research.collective_completion_20260914 import atoms as atom_tools
        maps,map_ids,atom_receipt=atom_tools.initial(meta,groups,maps,seed) if seed else atom_tools.from_guide(groups,maps)
        dump(base/'atoms.json',{'source':str(Path(seed).parent/'certificate.json') if seed else meta['state'],'source_kind':'compact_certificate' if seed else 'MPS_quartic_guide','initial':atom_receipt,'quadratic_blocks_remain_full':True,'maximum_atoms':4000})
        if enrich_seed and seed:maps,map_ids,guide=atom_tools.add(meta,groups,maps,map_ids,prepared,raw['dual'],add_rank)
    if share_spin:
        if not spin_twirl or not compact or response or rounds!=1:raise ValueError('Spin sharing requires one spin-averaged paired solve')
        from research.collective_completion_20260914.spin_share import share
        maps,map_ids,sharing=share(groups,maps);dump(base/'spin_sharing.json',sharing)
    for rnd in range(rounds):
        start=time.monotonic();folder=base/f'round_{rnd}';folder.mkdir()
        active,projected,original=build_maps(groups,rows,prepared,maps,compact,map_ids)
        if spin_twirl:projected=[((Ts@P).tocsc(),ii,jj) for P,ii,jj in projected]
        x=cp.Variable(free.shape[1]);expr=free@x;variables=[]
        for members,(P,ii,jj) in zip(maps,projected):
            r=members[0][0][1].shape[1];q=cp.Variable((r,r),PSD=True);variables.append(q);expr+=P@q[ii,jj]
        norm=np.asarray(free.power(2).sum(axis=1)).ravel()
        for P,ii,jj in projected:norm+=np.asarray(P.power(2).sum(axis=1)).ravel()
        row_scale=1/np.maximum(1,np.sqrt(norm));eq=cp.multiply(row_scale,expr-rhs)==0
        problem=cp.Problem(cp.Maximize(x[0]),[eq])
        if moment_dual:
            from research.collective_completion_20260914.moment_dual import problem as dual_problem
            problem,dual_eq,dual_y,dual_cones=dual_problem(free,rhs,projected,row_scale)
        build_seconds=time.monotonic()-start
        details={'case':case,'round':rnd,'spin_lifting':spin,'rank_cap_initial':rank,'gram_entries':sum(q.shape[0]**2 for q in variables),
            'gram_dimensions':[q.shape[0] for q in variables],'full_gram_entries':sum(len(g['words'])**2 for g in groups),
            'coefficient_equations':len(smallrows),'free_variables':free.shape[1],
            'projected_map_nonzeros':sum(P.nnz for P,_,_ in projected),'build_seconds':build_seconds,
            'prepared_dependency':str(prepared),'fixture':meta['fixture'],'state':meta['state'],
            'full_Gram_teacher_used':False,'many_body_states_enumerated':0,'guide':guide,'magnetization':mag,'singlet':singlet,'ladder_ideal_dimension':len(ladder_basis)}
        if spin_twirl:details.update(spin_row_projection=spin_row_receipt,coefficient_equations=len(selected))
        details.update(shared_spin_orbits=share_spin,map_ids=map_ids,moment_dual_formulation=moment_dual,coherent_atoms=atoms)
        dump(folder/'construction.json',details);print(json.dumps({'stage':'solve',**{k:v for k,v in details.items() if k not in ('guide','gram_dimensions','map_ids')}}),flush=True)
        t=time.monotonic()
        options={'eps':1e-8,'max_iters':100000,'time_limit_secs':seconds,'use_indirect':False} if solver=='SCS' else {'tol_gap_abs':1e-9,'tol_feas':1e-9,'tol_gap_rel':1e-9,'max_iter':300,'time_limit':seconds,'max_threads':threads}
        problem.solve(solver=solver,verbose=verbose,**options)
        details.update(solve_seconds=time.monotonic()-t,solver=solver,solver_status=problem.status,solver_iterations=problem.solver_stats.num_iters,solver_max_threads=threads if solver=='CLARABEL' else None)
        x_value=-np.asarray(dual_eq.dual_value) if moment_dual and dual_eq.dual_value is not None else (None if moment_dual else x.value)
        q_values=[q.dual_value for q in dual_cones] if moment_dual else [q.value for q in variables]
        if x_value is None or any(q is None for q in q_values):
            details['status']='no_candidate';dump(folder/'discovery.json',details);break
        xv=np.array(x_value);y=np.array(dual_y.value) if moment_dual else np.array(eq.dual_value)*row_scale
        if spin_twirl:y=np.asarray(Ts.T@y)
        np.savez_compressed(folder/'raw.npz',x=xv,dual=y,**{f'Q_{i}':q for i,q in enumerate(q_values)},**{f'V_{i}':mem[0][1] for i,(mem,name) in enumerate(maps)})
        coeff=np.asarray(fullfree@xv);blocks=[];min_ev=0.;actions=[]
        for k,((members,name),Q,(P,ii,jj)) in enumerate(zip(maps,q_values,projected)):
            q=(Q+Q.T)/2
            J=members[0][1];coeff+=original[k]@(J@q@J.T).ravel()
            ev,W=np.linalg.eigh(q);min_ev=min(min_ev,float(ev[0]))
            R=np.sqrt(np.maximum(ev,0))[:,None]*W.T
            C=members[0][1].T@Cs[f'C_{map_ids[k]}']@members[0][1] if compact else sum(V.T@Cs[f'C_{i}']@V for i,V in members)
            actions.append(float(np.trace(q@C)))
            if response:
                J=members[0][1];abstract=np.rint((R@J.T)*1e9).astype(np.int64);pair=meta['pairs'][map_ids[k]]
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
        combined={}
        for block in blocks:
            key=tuple(block['words'])
            if key not in combined:combined[key]={**block,'factor':[]}
            combined[key]['factor'].extend(block['factor'])
        blocks=list(combined.values())
        def poly(values):return add(*(scale(p,F(round(float(v)*1e10),10**10)) for p,v in zip(basis,values)))
        X=poly(xv[1:1+len(basis)]);Y=poly(xv[1+len(basis):1+2*len(basis)]) if spin else {}
        a=F(round(float(xv[1+2*len(basis)])*1e10),10**10) if singlet else F(0)
        ladder_poly=add(*(scale(p,F(round(float(v)*1e10),10**10)) for p,v in zip(ladder_basis,xv[2+2*len(basis):]))) if ladder else {}
        ladder_term=ladder_ideal(m,ladder_poly) if ladder else {}
        core={'modes':m,'particles':n,'operator_degree':3,'hamiltonian':encode(add(target,scale(product(Zspin,Y),-1),scale(S2,-a),scale(ladder_term,-1))) if spin else encode(h),
              'b':str(F(round(float(xv[0])*1e12),10**12)),'number_multiplier':encode(X),'denominator':10**9,'blocks':blocks}
        cert={'kind':'balanced_spin_sos_v1','modes':m,'particles':n,'hamiltonian':encode(h),'alpha_multiplier':encode(Y),'core':core} if spin else core
        if mag or singlet:cert.update(kind='spin_sector_sos_v1',magnetization=mag,singlet=singlet,casimir_multiplier=str(a))
        if ladder:cert['spin_ladder_multiplier']=encode(ladder_poly)
        if spin_twirl:cert['spin_twirl']=True
        (folder/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
        weights=np.array([1 if tuple(i for c,i in w if c)==tuple(i for c,i in w if not c) else 2 for w in smallrows])
        difference=coeff-fullrhs
        if spin_twirl:difference=T@difference
        details.update(status='requires_independent_exact_replay',proposed_b=float(xv[0]),numeric_coefficient_l1=float(weights@abs(difference)),
            minimum_PSD_eigenvalue=min_ev,total_state_action_Ha=sum(actions),factor_rows=sum(len(b['factor']) for b in blocks),
            certificate_bytes=(folder/'certificate.json').stat().st_size,round_seconds=time.monotonic()-start,total_seconds=time.monotonic()-total)
        dump(folder/'discovery.json',details);print(json.dumps({k:v for k,v in details.items() if k not in ('guide','gram_dimensions','map_ids')}),flush=True)
        if rnd+1<rounds:
            if atoms:maps,map_ids,guide=atom_tools.add(meta,groups,maps,map_ids,prepared,y,add_rank)
            else:maps,guide=enrich(maps,groups,original,y,add_rank)
            dump(folder/'enrichment.json',guide)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('rank',type=int);p.add_argument('--seconds',type=int,default=180);p.add_argument('--spin',action='store_true');p.add_argument('--tag',default='');p.add_argument('--rounds',type=int,default=1);p.add_argument('--add-rank',type=int,default=8);p.add_argument('--prepared');p.add_argument('--solver',choices=['SCS','CLARABEL'],default='SCS');p.add_argument('--seed');p.add_argument('--enrich-seed',action='store_true');p.add_argument('--mag',type=int,default=0);p.add_argument('--singlet',action='store_true');p.add_argument('--ladder',action='store_true');p.add_argument('--seed-eigen-cutoff',type=float,default=0.);p.add_argument('--threads',type=int,default=1);p.add_argument('--spin-twirl',action='store_true');p.add_argument('--share-spin',action='store_true');p.add_argument('--moment-dual',action='store_true');p.add_argument('--verbose',action='store_true');p.add_argument('--atoms',action='store_true');a=p.parse_args()
    solve(a.case,a.rank,a.seconds,a.spin,a.tag,a.rounds,a.add_rank,a.prepared,a.solver,a.seed,a.enrich_seed,a.mag,a.singlet,a.ladder,a.seed_eigen_cutoff,a.threads,a.spin_twirl,a.share_spin,a.moment_dual,a.verbose,a.atoms)
