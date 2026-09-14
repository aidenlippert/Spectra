"""Fresh bare-H joint metric/atom discovery in an invariant coefficient space.

Free metric variables are represented by signed pairs in a native incremental
phase-I LP. No fixed metric or accepted atom list is imported.
"""
from pathlib import Path
from fractions import Fraction as F
from unittest.mock import patch
import argparse,json,time
import numpy as np
from scipy.sparse import csc_matrix,hstack
from experiments.marginal_fixed_metric_pricing import IncrementalMaster
from experiments.marginal_reynolds_pricing import ReynoldsMomentDictionary
from experiments.marginal_joint_coefficient_constructor import input_digest,export
from experiments.marginal_joint_spinflip import flip_label


def price_joint(dictionary,dual,active,batch,weight_degree):
    if weight_degree==6:return dictionary.price(dual,active,batch)
    if weight_degree!=4:raise ValueError('Weight degree must be4 or6')
    r=dictionary.qrows;original=dictionary.supports
    dictionary.supports=[(degree,indices,[(j,label) for j,label in admitted
        if degree<=(4 if label[0]=='positive' else 2)])
        for degree,indices,admitted in original if degree<=4]
    try:
        weight,wm,wc=dictionary.price(np.r_[dual[:r],np.zeros(r+1)],active,batch)
    finally:dictionary.supports=original
    numerator,nm,nc=dictionary.price(np.r_[np.zeros(r),dual[r:2*r],0.],active,batch)
    result=weight+numerator
    if any(b==0 and label[1].bit_count()>(4 if label[0]=='positive' else 2) for b,label in result):
        raise ValueError('Weight atom degree restriction failed')
    return result,max(wm,nm),wc+nc


def construct(data,gamma,out,time_limit=240,batch=128,hard_normalization=False,method="simplex",audit=False,coordinates="coefficient",weight_degree=6,resume=None,refine=False,allow_coordinate_change=False,crossover=False,weight_floor=.0009,numerator_floor=.0001):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=ReynoldsMomentDictionary(data,F(gamma),feature_degree=2,proof_degree=6,coordinates=coordinates)
    rows,n=d.metric.shape
    if weight_floor<=0 or numerator_floor<=0: raise ValueError('Positive floors required')
    rhs=np.r_[weight_floor*d.one,numerator_floor*d.one,1.]
    solver=IncrementalMaster(rhs)
    if method not in ('simplex','ipm'):raise ValueError('Supported native LP method required')
    if method=='ipm':
        for key,value in [('solver','ipm'),('run_crossover','on' if crossover else 'off')]:
            if solver.solver.setOptionValue(key,value)==solver.highspy.HighsStatus.kError:
                raise ValueError('Native LP option rejected: '+key)
    if refine:
        if method!='ipm':raise ValueError('Refinement requires the interior-point solver')
        for key,value in [('ipm_optimality_tolerance',1e-12),('primal_feasibility_tolerance',1e-10),('dual_feasibility_tolerance',1e-10)]:
            if solver.solver.setOptionValue(key,value)==solver.highspy.HighsStatus.kError:
                raise ValueError('Native refinement tolerance rejected: '+key)
    if hard_normalization:
        # Metric variables can satisfy the mean constraint; forbid artificial
        # residuals in this row so the zero metric is not a phase-I plateau.
        for column in (rows-1,2*rows-1):
            if solver.solver.changeColBounds(column,0.,0.)==solver.highspy.HighsStatus.kError:
                raise ValueError('Hard normalization bound rejected')
    active=[(b,label) for b in (0,1) for label in d.labels(2)]
    signature={'source_sha256':input_digest(data),'gamma':str(F(gamma)),
        'coordinates':coordinates,'weight_degree':weight_degree,'hard_normalization':hard_normalization}
    resume_source=None;previous_rounds=0
    if resume is not None:
        import hashlib
        from math import comb
        raw=Path(resume).read_bytes();checkpoint=json.loads(raw)
        expected=dict(signature);old_signature=checkpoint.get('signature')
        if allow_coordinate_change and isinstance(old_signature,dict) and old_signature.get('coordinates') in ('coefficient','conditional'):
            expected['coordinates']=old_signature['coordinates']
        if old_signature!=expected:raise ValueError('Checkpoint search signature mismatch')
        saved=checkpoint.get('active')
        if type(saved) is not list or not 1<=len(saved)<=50000:raise ValueError('Bounded checkpoint column list required')
        active=[]
        for item in saved:
            if type(item) is not list or len(item)!=2:raise ValueError('Invalid checkpoint column')
            b,label=item
            if type(b) is not int or b not in (0,1) or type(label) is not list or len(label)!=3:raise ValueError('Invalid checkpoint atom')
            family,R,O=label
            if family not in ('positive','charge') or type(R) is not int or type(O) is not int or not 0<=O<=R<1<<d.modes or O&~R:raise ValueError('Invalid checkpoint masks')
            limit=(weight_degree if b==0 else 6)-(2 if family=='charge' else 0)
            if R.bit_count()>limit:raise ValueError('Checkpoint atom degree exceeds block limit')
            count=1
            for spin in d.ring.spin_masks:
                available=d.sites-(R&spin).bit_count();needed=d.ring.target-(O&spin).bit_count()
                count*=comb(available,needed) if 0<=needed<=available else 0
            if count<=(1 if family=='positive' else 0) or tuple(label)!=min(d.atom_members(label)):raise ValueError('Checkpoint atom is not admissible and canonical')
            active.append((b,tuple(label)))
        if len(set(active))!=len(active):raise ValueError('Duplicate checkpoint columns')
        previous_rounds=checkpoint['completed_rounds']
        if type(previous_rounds) is not int or previous_rounds<0:raise ValueError('Invalid checkpoint round count')
        resume_source={'path':str(resume),'sha256':hashlib.sha256(raw).hexdigest(),'columns':len(active),'completed_rounds':previous_rounds,'coordinates_from':old_signature['coordinates'],'coordinates_to':coordinates}
    seen=set(active);built=0
    def save_checkpoint(completed_rounds):
        (out/'checkpoint.json').write_text(json.dumps({'signature':signature,'completed_rounds':previous_rounds+completed_rounds,'active':active,'scope':'Selected atom directions only, from this bare-H search lineage; no imported metric or proof weights. Native LP is rebuilt on resume.'},indent=2)+'\n')
    save_checkpoint(0)
    history=[];accepted=False;reason='Round budget reached';prepared=time.monotonic()
    for iteration in range(100):
        remaining=time_limit-(time.monotonic()-prepared)
        if remaining<=0:reason='Search time budget reached';break
        columns=[]
        for block,label in active[built:]:
            c=np.zeros(rows);c[block*d.qrows:(block+1)*d.qrows]=-d.column(label);columns.append(c)
        fresh=csc_matrix(np.array(columns).T) if columns else csc_matrix((rows,0))
        if iteration==0:fresh=hstack((d.metric,-d.metric,fresh),format='csc')
        result=solver.run(fresh,remaining);built=len(active)
        record={'round':previous_rounds+iteration,'active_atoms':built,'search_seconds':time.monotonic()-prepared,'status':result.message}
        if not result.success:reason='Restricted solver did not finish';history.append(record);break
        record['phase_one_l1']=float(result.fun)
        if audit:
            coefficients=result.x[:n]-result.x[n:2*n]
            info=solver.solver.getInfo()
            record.update(metric_coefficients=coefficients.tolist(),metric_mean=float((d.metric@coefficients)[-1]),nonzero_atom_weights=int(np.sum(result.x[2*n:2*n+built]>1e-9)),simplex_iterations=info.simplex_iteration_count,ipm_iterations=info.ipm_iteration_count,
                primal_infeasibility=info.max_primal_infeasibility,dual_infeasibility=info.max_dual_infeasibility)

        if result.fun<1e-7:
            proposal={'source_sha256':input_digest(data),'gamma':str(F(gamma)),'orbits':d.orbits,'metric_coefficients':(result.x[:n]-result.x[n:2*n]).tolist(),'metric_bound':weight_floor,'numerator_bound':numerator_floor}
            for block,name in enumerate(('weight','numerator')):
                pairs=[]
                for (b,label),value in zip(active,result.x[2*n:2*n+built]):
                    if b!=block or abs(value)<=1e-13:continue
                    members=d.atom_members(label)
                    for member in members:
                        spin=sorted({tuple(member),tuple(flip_label(member,d.sites))})
                        pairs.extend((list(m),float(value)/(len(members)*len(spin))) for m in spin)
                proposal[name+'_labels']=[p[0] for p in pairs];proposal[name+'_values']=[p[1] for p in pairs]
            (out/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
            try:
                precision={'metric_precision':10**12,'proof_denominator':10**14} if refine else {}
                record['exact_receipt']=export(data,proposal,out/'proof',**precision);accepted=True;reason='Exact certificate accepted'
            except ValueError as error:reason='Exact export rejected';record['export_rejection']=str(error)
            history.append(record);print(json.dumps(record),flush=True);break
        fresh,maximum,checked=price_joint(d,result.eqlin.marginals,seen,batch,weight_degree)
        if audit:
            direct=[-float(result.eqlin.marginals[b*d.qrows:(b+1)*d.qrows]@d.column(label)) for b,label in fresh]
            if direct and min(direct)<=1e-8:raise ValueError('Priced atom is not directly improving')
            record.update(minimum_selected_direct_violation=min(direct,default=0),maximum_selected_direct_violation=max(direct,default=0),selected_degrees={str(degree):sum(label[1].bit_count()==degree for b,label in fresh) for degree in range(7)})

        record.update(maximum_dual_violation=maximum,candidates_priced=checked,added_atoms=len(fresh));history.append(record);print(json.dumps(record),flush=True)
        (out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
        if not fresh:reason='No improving column; no infeasibility proof';break
        active.extend(fresh);seen.update(fresh);save_checkpoint(len(history))
    save_checkpoint(len(history))
    receipt=dict(d.stats,source_sha256=input_digest(data),gamma=str(F(gamma)),exact_accepted=accepted,reason=reason,hard_normalization=hard_normalization,method=method,audit=audit,weight_degree=weight_degree,refinement=refine,crossover=crossover,weight_floor=weight_floor,numerator_floor=numerator_floor,resume_source=resume_source,metric_features=n,rounds=len(history),active_atoms=len(active),preparation_seconds=prepared-start,total_seconds=time.monotonic()-start,scope='Fresh bare-H invariant metric and atom discovery. Geometry-only low-degree seeds, direct moment pricing, no imported candidate or atom directions. Exact full original export required; numerical symmetry rank never establishes acceptance.')
    if resume_source:receipt['scope']='Resumed atom directions from an earlier bare-H run with identical Hamiltonian, target, normalization and block settings; any explicit coordinate change is recorded. No imported metric or accepted proof weights; exact original export required.'
    (out/'history.json').write_text(json.dumps(history,indent=2)+'\n');(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--gamma',default='-14');p.add_argument('--time-limit',type=float,default=240);p.add_argument('--out',type=Path,required=True);p.add_argument('--hard-normalization',action='store_true');p.add_argument('--method',choices=('simplex','ipm'),default='simplex');p.add_argument('--audit',action='store_true');p.add_argument('--coordinates',choices=('coefficient','conditional'),default='coefficient');p.add_argument('--weight-degree',choices=(4,6),type=int,default=6);p.add_argument('--resume',type=Path);p.add_argument('--refine',action='store_true');p.add_argument('--allow-coordinate-change',action='store_true');p.add_argument('--crossover',action='store_true');p.add_argument('--weight-floor',type=float,default=.0009);p.add_argument('--numerator-floor',type=float,default=.0001);a=p.parse_args()
    data=json.loads(Path('results/marginal_graded_hubbard8/hamiltonian.json').read_text())
    with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No states')),patch('experiments.marginal_spin_constructor.spin_states',side_effect=AssertionError('No state generation')),patch('experiments.marginal_polynomial_metric.complete_number_ideals',side_effect=AssertionError('No full-population lift')),patch('experiments.marginal_joint_coefficient_constructor.prepare',side_effect=AssertionError('No full atom dictionary')):
        construct(data,a.gamma,a.out,time_limit=a.time_limit,hard_normalization=a.hard_normalization,method=a.method,audit=a.audit,coordinates=a.coordinates,weight_degree=a.weight_degree,resume=a.resume,refine=a.refine,allow_coordinate_change=a.allow_coordinate_change,crossover=a.crossover,weight_floor=a.weight_floor,numerator_floor=a.numerator_floor)
