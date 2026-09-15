"""Sparse-map PSD block discovery; exact exports, explicit full pricing costs."""
import argparse,json,sys,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from scipy import sparse
import cvxpy as cp
from experiments.marginal_symbolic import canonical,hermitian,decode,encode,mono,add,scale,product,number_shift,multiplier_basis,verify
from experiments.marginal_coefficient import dictionaries,gram_map
from experiments.marginal_molecule_stress import parity_generators
from research.certificate_scaling.direct_sparse_discovery import candidates,sparse_columns
from research.certificate_scaling.adaptive_factor_pricing import pricing_vectors


def partition(h,modes,family,symmetry):
    original=dictionaries(modes,family)
    masks=parity_generators(h,modes) if symmetry else []
    spin=symmetry and all(sum(2*c-1 for c,i in w if i%2==0)==0 for w in h)
    def signature(w):
        return ((sum(2*c-1 for c,i in w if i%2==0),) if spin else ())+tuple(sum((mask>>i)&1 for c,i in w)%2 for mask in masks)
    groups=[]
    for block in original:
        parts={}
        for w in block['words']:parts.setdefault(signature(w),[]).append(w)
        groups.extend({'name':block['name']+str(k),'words':ws} for k,ws in sorted(parts.items()))
    stats={'original_dimensions':[len(b['words']) for b in original],
           'symmetry_dimensions':[len(b['words']) for b in groups],
           'parity_masks':masks,'alpha_charge_quotient':bool(spin)}
    return groups,signature,stats


def run(h,modes,particles,outdir,width=4,rounds=4,seconds=180,denominator=10**7,
        family='quadratic',symmetry=False,full=False,batch=24,ideal_body=1,
        solver='SCS',solver_eps=1e-7,solver_seconds=60.,solver_max_iters=20000,
        save_raw=False,condition='none',free_denominator=10**9,quotient_linear=False,anchor_identity=False,formulation='primal',
        save_scs_state=False,scs_resume=None,scs_stages=1):
    start=time.monotonic();h=canonical(h);outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    if not hermitian(h) or any(len(w)>4 for w in h):raise ValueError('Hermitian two-body H required')
    if width<2 or rounds<1 or seconds<=0:raise ValueError('Invalid budgets')
    if solver not in ('SCS','CLARABEL') or condition not in ('none','row','row-column'):
        raise ValueError('Unsupported solver or conditioning')
    if formulation not in ('primal','dual'):raise ValueError('Unknown formulation')
    if type(scs_stages) is not int or scs_stages<1:raise ValueError('Positive SCS stage count required')
    checkpointing=save_scs_state or scs_resume is not None or scs_stages>1
    if checkpointing and (solver!='SCS' or not full):raise ValueError('SCS checkpoints require a fixed full dictionary')
    if solver_eps<=0 or solver_seconds<=0 or solver_max_iters<1 or denominator<=0 or free_denominator<=0:
        raise ValueError('Invalid solver/export precision')
    groups,sig,stats=partition(h,modes,family,symmetry);identity_sig=sig(())
    if quotient_linear:
        if family!='mixed' or ideal_body!=2 or not full or not 2<=particles<=modes:
            raise ValueError('Linear quotient requires full mixed cone, body-two ideal, and N>=2')
        stats['prequotient_symmetry_dimensions']=stats['symmetry_dimensions']
        groups=[{'name':g['name'],'words':[w for w in g['words'] if len(w)!=1]} for g in groups]
        groups=[g for g in groups if g['words']]
        stats['symmetry_dimensions']=[len(g['words']) for g in groups]
    degree=2 if family=='quadratic' and ideal_body==1 else 3
    rows=[tuple((1,i) for i in l)+tuple((0,i) for i in r) for k in range(degree+1)
          for l in combinations(range(modes),k) for r in combinations(range(modes),k)
          if l<=r and sig(tuple((1,i) for i in l)+tuple((0,i) for i in r))==identity_sig]
    lookup={w:i for i,w in enumerate(rows)}
    weights=np.array([1 if tuple(i for c,i in w if c)==tuple(i for c,i in w if not c) else 2 for w in rows])
    maps=[gram_map(g['words'],lookup).tocsc() for g in groups]
    basis=[p for p in multiplier_basis(modes,max_body=ideal_body) if all(sig(w)==identity_sig for w in p)]
    free=sparse_columns([{w:c for w,c in p.items() if w in lookup} for p in [mono(())]+[product(number_shift(modes,particles),q) for q in basis]],lookup)
    rhs=np.array([float(h.get(w,0)) for w in rows]);supports=[];seen=set()
    row_scale=np.ones(len(rows));free_scale=np.ones(free.shape[1])
    if condition!='none':
        squared=np.asarray(free.power(2).sum(axis=1)).ravel()
        for gmap in maps:squared+=np.asarray(gmap.power(2).sum(axis=1)).ravel()
        row_scale=1/np.maximum(1.,np.sqrt(squared))
    if condition=='row-column':
        free_scale=1/np.maximum(1.,np.sqrt(np.asarray(free.power(2).sum(axis=0)).ravel()))
    def insert(gid,indices):
        key=(gid,tuple(sorted(set(map(int,indices)))))
        if key in seen:return False
        seen.add(key);supports.append(key);return True
    for gid,g in enumerate(groups):
        n=len(g['words'])
        if full:insert(gid,range(n))
        else:
            for i in range(0,n,width):insert(gid,range(i,min(n,i+width)))
    if not full:
        seed,_=candidates(h,modes,512)
        location={w:(gid,i) for gid,g in enumerate(groups) for i,w in enumerate(g['words'])}
        for words,_,_ in seed:
            loc=[location.get(w) for w in words]
            if all(x is not None for x in loc) and len({x[0] for x in loc})==1:insert(loc[0][0],[x[1] for x in loc])
    map_seconds=time.monotonic()-start;history=[];best=None;total_price=0.;total_solve=0.;total_export=0.;stop='round_budget'
    for rnd in range(1 if full else rounds):
        if best and time.monotonic()-start>=seconds:stop='wall_budget';break
        before=time.monotonic();x=cp.Variable(free.shape[1]);physical_x=cp.multiply(free_scale,x)
        rem=cp.Variable(len(rows));expr=free@physical_x;variables=[];active_maps=[]
        for gid,inds in supports:
            n=len(groups[gid]['words']);ix=np.array(inds);cols=(ix[:,None]*n+ix[None,:]).ravel();k=len(ix)
            q=cp.Variable((k,k),PSD=True);variables.append(q)
            active_maps.append(maps[gid][:,cols]);expr += active_maps[-1]@cp.reshape(q,(k*k,),order='C')
        equality=cp.multiply(row_scale,expr+rem-rhs)==0
        constraints=[equality]+([rem[lookup[()]]==0] if anchor_identity else [])
        problem=cp.Problem(cp.Minimize(weights@cp.abs(rem)-physical_x[0]),constraints)
        if formulation=='dual':
            dual_var=cp.Variable(len(rows));physical_y=cp.multiply(row_scale,dual_var)
            unit=np.zeros(free.shape[1]);unit[0]=1.
            free_equality=cp.multiply(free_scale,free.T@physical_y-unit)==0
            # Anchoring rem[I]=0 removes the box on y[I]. Without anchoring,
            # that box is redundant anyway: the free b column fixes y[I]=1.
            bounded=np.array([i for i in range(len(rows)) if i!=lookup[()]])
            box_plus=physical_y[bounded]<=weights[bounded];box_minus=-physical_y[bounded]<=weights[bounded]
            cone_constraints=[]
            for (_,inds),gmap in zip(supports,active_maps):
                matrix=cp.reshape(gmap.T@physical_y,(len(inds),len(inds)),order='C')
                cone_constraints.append((matrix+matrix.T)/2 >> 0)
            problem=cp.Problem(cp.Minimize(rhs@physical_y),[free_equality,box_plus,box_minus]+cone_constraints)
        time_limit=max(1.,min(solver_seconds,seconds-(time.monotonic()-start)))
        options=({'eps':solver_eps,'max_iters':solver_max_iters,'time_limit_secs':time_limit} if solver=='SCS' else
                 {'tol_feas':solver_eps,'tol_gap_abs':solver_eps,'tol_gap_rel':solver_eps,'max_iter':solver_max_iters,'time_limit':time_limit})
        stage_receipts=[]
        if checkpointing:
            from research.certificate_scaling.scs_checkpoint import solve as checkpoint_solve
            resume=scs_resume
            contract={'modes':modes,'particles':particles,'hamiltonian':encode(h),'family':family,
                      'symmetry':symmetry,'ideal_body':ideal_body,'condition':condition,
                      'quotient_linear':quotient_linear,'anchor_identity':anchor_identity,'formulation':formulation}
            for stage in range(scs_stages):
                if stage and time.monotonic()-start>=seconds:break
                options['time_limit_secs']=max(1.,min(solver_seconds,seconds-(time.monotonic()-start)))
                prefix=outdir/f'scs_round_{rnd}_stage_{stage}'
                stage_receipts.append(checkpoint_solve(problem,options,prefix,resume,contract))
                resume=prefix
                (outdir/'scs_stages.json').write_text(json.dumps(stage_receipts,indent=2)+'\n')
        else:
            problem.solve(solver=solver,**options)
        solve_wall=time.monotonic()-before;total_solve+=solve_wall
        if formulation=='primal':
            if x.value is None or equality.dual_value is None or any(q.value is None for q in variables):
                stop='solver_failure';break
            y=np.asarray(equality.dual_value)*row_scale;xvalues=np.asarray(x.value)*free_scale
            raw_grams=[(q.value+q.value.T)/2 for q in variables];remvalues=np.asarray(rem.value)
        else:
            if dual_var.value is None or free_equality.dual_value is None or any(c.dual_value is None for c in cone_constraints):
                stop='solver_failure';break
            y=np.asarray(dual_var.value)*row_scale;xvalues=-np.asarray(free_equality.dual_value)*free_scale
            raw_grams=[(c.dual_value+c.dual_value.T)/2 for c in cone_constraints]
            remvalues=np.zeros(len(rows));remvalues[bounded]=box_minus.dual_value-box_plus.dual_value
        before=time.monotonic();factors=[]
        if save_raw:
            np.savez_compressed(outdir/f'raw_round_{rnd}.npz',x=xvalues,rem=remvalues,dual=y,
                                **{f'gram_{i}':q for i,q in enumerate(raw_grams)})
            (outdir/f'raw_round_{rnd}.json').write_text(json.dumps({'modes':modes,'particles':particles,'hamiltonian':encode(h),
                'family':family,'symmetry':symmetry,'ideal_body':ideal_body,'solver':solver,'solver_options':options,'condition':condition,
                'quotient_linear':quotient_linear,'anchor_identity':anchor_identity,'formulation':formulation,
                'rows':rows,'weights':weights.tolist(),'supports':[{'group':gid,'indices':list(inds),'words':[groups[gid]['words'][i] for i in inds]} for gid,inds in supports]},separators=(',',':'))+'\n')
        raw_coeff=np.asarray(free@xvalues);clip_coeff=raw_coeff.copy();rounded_coeff=raw_coeff.copy()
        min_eigen=0.;negative_mass=0.;negative_count=0;max_gram=0.
        for (gid,inds),q,gmap in zip(supports,raw_grams,active_maps):
            raw_coeff+=gmap@q.ravel()
            val,vec=np.linalg.eigh(q);min_eigen=min(min_eigen,float(val[0]));negative_mass+=float(-val[val<0].sum());negative_count+=int((val<0).sum())
            max_gram=max(max_gram,float(np.max(abs(q))))
            root=np.sqrt(np.maximum(val,0))[:,None]*vec.T
            z=np.rint(root*denominator).astype(np.int64)
            clip_coeff+=gmap@(root.T@root).ravel()
            rounded_coeff+=gmap@((z.astype(float)/denominator).T@(z.astype(float)/denominator)).ravel()
            z=z[np.any(z,axis=1)]
            if len(z):factors.append({'name':str(len(factors)),'words':[groups[gid]['words'][i] for i in inds],'factor':z.tolist()})
        rounded_x=np.rint(xvalues*free_denominator)/free_denominator
        free_delta=np.asarray(free@(rounded_x-xvalues));export_coeff=rounded_coeff+free_delta
        X=add(*(scale(q,F(round(float(v)*free_denominator),free_denominator)) for q,v in zip(basis,xvalues[1:])))
        cert={'modes':modes,'particles':particles,'hamiltonian':encode(h),'b':str(F(round(float(xvalues[0])*free_denominator),free_denominator)),
              'number_multiplier':encode(X),'denominator':denominator,'blocks':factors}
        checked=verify(cert);total_export+=time.monotonic()-before
        if best is None or F(checked['lower'])>F(best[1]['lower']):best=(cert,checked)
        row={'round':rnd,'supports_solved':len(supports),'gram_scalar_entries':sum(len(ix)**2 for _,ix in supports),
             'numeric_lower':float(xvalues[0]-weights@abs(remvalues)),'optimization_objective':float(problem.value),
             'exact_lower':checked['lower_float'],'dual_identity':float(y[lookup[()]]),
             'solver_status':problem.status,'solve_wall_seconds':solve_wall,
             'raw_b':float(xvalues[0]),'solver_residual_l1':float(weights@abs(remvalues)),
             'raw_coefficient_residual_l1':float(weights@abs(rhs-raw_coeff)),
             'equality_violation_l1':float(weights@abs(raw_coeff+remvalues-rhs)),
             'minimum_raw_eigenvalue':min_eigen,'negative_eigenvalue_mass':negative_mass,'negative_eigenvalue_count':negative_count,
             'maximum_gram_entry':max_gram,'maximum_free_coefficient':float(max(abs(xvalues))),
             'clipped_residual_l1':float(weights@abs(rhs-clip_coeff)),
             'clipping_coefficient_delta_l1':float(weights@abs(clip_coeff-raw_coeff)),
             'factor_rounding_delta_l1':float(weights@abs(rounded_coeff-clip_coeff)),
             'free_rounding_delta_l1':float(weights@abs(free_delta)),
             'export_numeric_residual_l1':float(weights@abs(rhs-export_coeff)),
             'exact_residual_l1':checked['residual_l1_float'],'solver_iterations':problem.solver_stats.num_iters,
             'scs_checkpoint_stages':len(stage_receipts),
             'scs_cumulative_iterations':stage_receipts[-1]['cumulative_iterations'] if stage_receipts else None}
        history.append(row)
        (outdir/'certificate.json').write_text(json.dumps(best[0],separators=(',',':'))+'\n')
        (outdir/'history.json').write_text(json.dumps(history,indent=2)+'\n')
        if full or rnd+1==rounds:break
        before=time.monotonic();proposals=[]
        for gid,gmap in enumerate(maps):
            n=len(groups[gid]['words']);mat=np.asarray(gmap.T@y).reshape(n,n);mat=(mat+mat.T)/2
            for value,ix,_ in pricing_vectors(mat,width,batch):proposals.append((value,gid,ix))
        added=0
        for value,gid,ix in sorted(proposals,key=lambda p:p[0]):
            if insert(gid,ix):added+=1
            if added>=batch:break
        total_price+=time.monotonic()-before;row['supports_added']=added
        row['minimum_pricing_value']=min((p[0] for p in proposals),default=None)
        if not added:stop='no_new_heuristic_support';break
    if best is None:raise RuntimeError('No certificate produced')
    cert,rec=best;rec.update(stats);rec.update({'method':'joint_PSD_blocks','full_blocks':full,'symmetry':symmetry,'family':family,'width':width,
        'rounds_solved':len(history),'retained_supports':len(supports),'stop':stop,'coefficient_rows':len(rows),'ideal_body':ideal_body,
        'complete_pricing_matrix_entries':sum(m.shape[1] for m in maps),'complete_pricing_map_nonzeros':sum(m.nnz for m in maps),
        'map_build_seconds':map_seconds,'total_solve_wall_seconds':total_solve,'total_pricing_seconds':total_price,
        'total_exact_export_seconds':total_export,'wall_seconds':time.monotonic()-start,'source_factors_used':False,
        'source_upper_used':False,'many_body_space_enumerated':False,'omitted_family_optimality_proved':False})
    rec.update({'solver':solver,'solver_eps':solver_eps,'solver_seconds':solver_seconds,'solver_max_iters':solver_max_iters,
                'condition':condition,'raw_saved':save_raw,'factor_denominator':denominator,'free_denominator':free_denominator,
                'quotient_linear':quotient_linear,'anchor_identity':anchor_identity,'formulation':formulation,
                'scs_checkpointing':checkpointing,'scs_stages_requested':scs_stages})
    (outdir/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n');(outdir/'history.json').write_text(json.dumps(history,indent=2)+'\n')
    return rec

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fixture',type=Path,required=True);p.add_argument('--outputdir',type=Path,required=True)
    p.add_argument('--width',type=int,default=4);p.add_argument('--rounds',type=int,default=4);p.add_argument('--seconds',type=float,default=180)
    p.add_argument('--family',choices=['quadratic','mixed','local'],default='quadratic');p.add_argument('--ideal-body',type=int,default=1)
    p.add_argument('--symmetry',action='store_true');p.add_argument('--full',action='store_true')
    p.add_argument('--solver',choices=['SCS','CLARABEL'],default='SCS');p.add_argument('--solver-eps',type=float,default=1e-7)
    p.add_argument('--solver-seconds',type=float,default=60);p.add_argument('--solver-max-iters',type=int,default=20000)
    p.add_argument('--save-raw',action='store_true');p.add_argument('--condition',choices=['none','row','row-column'],default='none')
    p.add_argument('--denominator',type=int,default=10**7);p.add_argument('--free-denominator',type=int,default=10**9)
    p.add_argument('--quotient-linear',action='store_true');p.add_argument('--anchor-identity',action='store_true')
    p.add_argument('--formulation',choices=['primal','dual'],default='primal')
    p.add_argument('--save-scs-state',action='store_true');p.add_argument('--scs-resume',type=Path)
    p.add_argument('--scs-stages',type=int,default=1)
    a=p.parse_args();f=json.loads(a.fixture.read_text())
    print(json.dumps(run(decode(f['hamiltonian'],f['modes'],4),f['modes'],f['particles'],a.outputdir,a.width,a.rounds,a.seconds,
        family=a.family,symmetry=a.symmetry,full=a.full,ideal_body=a.ideal_body,solver=a.solver,solver_eps=a.solver_eps,
        solver_seconds=a.solver_seconds,solver_max_iters=a.solver_max_iters,save_raw=a.save_raw,condition=a.condition,
        denominator=a.denominator,free_denominator=a.free_denominator,quotient_linear=a.quotient_linear,anchor_identity=a.anchor_identity,formulation=a.formulation,
        save_scs_state=a.save_scs_state,scs_resume=a.scs_resume,scs_stages=a.scs_stages)))
