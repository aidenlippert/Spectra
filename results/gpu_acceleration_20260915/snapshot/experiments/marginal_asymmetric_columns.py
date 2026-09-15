"""Price new asymmetric quartic squares against a restricted LP moment dual."""
from fractions import Fraction as F
import json
from pathlib import Path
import time

import numpy as np
import scipy.sparse as sp
import cvxpy as cp

from experiments.marginal_asymmetric_polish import prepare_polish,solve_weights,export_weights
from experiments.marginal_coefficient import gram_map
from experiments.marginal_dual_exact import full_word_blocks
from experiments.marginal_symbolic import encode,decode

ROOT=Path(__file__).resolve().parents[1]


def prepare_candidates(data):
    lookup={w:i for i,w in enumerate(data['rows'])}
    blocks=full_word_blocks(data['modes'],4)
    return blocks,[data['projection']@gram_map(b['words'],lookup) for b in blocks]


def price(blocks,maps,dual,pool):
    proposals=[]
    for i,(block,g) in enumerate(zip(blocks,maps)):
        n=len(block['words']);matrix=np.asarray(g.T@dual).reshape(n,n)
        values,vectors=np.linalg.eigh((matrix+matrix.T)/2)
        if values[0]<-1e-7:proposals.append((float(values[0]),i,vectors[:,0]))
    proposals.sort(key=lambda item:item[0]);result=[]
    for value,i,v in proposals[:pool]:
        p={w:F(round(float(c)*10**14),10**14) for w,c in zip(blocks[i]['words'],v) if round(float(c)*10**14)}
        rounded=np.array([float(p.get(w,0)) for w in blocks[i]['words']])
        column=np.asarray(maps[i]@np.outer(rounded,rounded).ravel()).ravel()
        result.append({'violation':value,'dictionary':i,'polynomial':p,'column':column})
    return result,len(proposals)


def center_dual(data,columns,objective,tolerance=1e-5,capped=False):
    n=len(data['unit']);matrix=sp.csr_matrix(np.column_stack(columns));y=cp.Variable(n)
    slack=cp.Variable(len(columns),nonneg=True) if capped else np.zeros(len(columns))
    constraints=[data['unit']@y==1,data['a'].T@y==0,matrix.T@y+slack>=0,
                 y[1:]<=data['weights'][1:],y[1:]>=-data['weights'][1:],
                 (data['rhs']-objective*data['unit'])@y+100*cp.sum(slack)<=tolerance]
    problem=cp.Problem(cp.Minimize(cp.sum_squares(cp.multiply(1/data['weights'],y))),constraints)
    started=time.monotonic()
    try:problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9,max_iter=150)
    except cp.error.SolverError:return None,{'status':'solver_error','seconds':time.monotonic()-started}
    if y.value is None or (capped and slack.value is None):return None,{'status':problem.status,'seconds':time.monotonic()-started}
    yy=np.array(y.value);ss=np.array(slack.value) if capped else slack
    if not (np.all(np.isfinite(yy)) and np.all(np.isfinite(ss))):
        return None,{'status':'nonfinite_dual','seconds':time.monotonic()-started,'accepted':False}
    error=max(abs(data['unit']@yy-1),np.max(abs(data['a'].T@yy)),max(0.,float(np.max(-matrix.T@yy-ss))),max(0.,float(-min(ss))),
              max(0.,float(np.max(abs(yy[1:])-data['weights'][1:]))))
    excess=float((data['rhs']-objective*data['unit'])@yy+100*sum(ss))
    status={'status':problem.status,'seconds':time.monotonic()-started,'constraint_error':float(error),'energy_excess':excess,'energy_tolerance':tolerance,'capped_dual':capped}
    accepted=bool(error<=1e-6 and excess<=tolerance+1e-7);status['accepted']=accepted
    return (yy if accepted else None),status


def run(source,iterations=16,pool=8,name='epsilon001',backend='native',matrix_threshold=1e-12,resume=None,center=False):
    if type(iterations) is not int or iterations<0 or type(pool) is not int or pool<1:raise ValueError('Invalid search budget')
    if Path(name).name!=name:raise ValueError('Name must be a directory basename')
    out=ROOT/'results/marginal_asymmetric_columns'/name
    if out.exists():raise ValueError('Preserve previous search; choose another name')
    out.mkdir(parents=True);started=time.monotonic()
    data=prepare_polish(source,average_direct=True);seeds=list(data['seeds']);columns=list(data['columns'])
    blocks,maps=prepare_candidates(data)
    cuts=[]
    if resume is not None:
        resume=Path(resume)
        previous=json.loads((resume.parent/'certificate.json').read_text())
        if previous['hamiltonian']!=data['original']['hamiltonian']:raise ValueError('Checkpoint Hamiltonian mismatch')
        for item in json.loads(resume.read_text()):
            i=item['dictionary'];p=decode(item['polynomial'],data['modes'],4)
            if type(i) is not int or not 0<=i<len(blocks) or not set(p).issubset(blocks[i]['words']):raise ValueError('Invalid checkpoint square')
            v=np.array([float(p.get(w,0)) for w in blocks[i]['words']])
            seeds.append(('exchange_squares',p));columns.append(np.asarray(maps[i]@np.outer(v,v).ravel()).ravel());cuts.append(item)
    structure={'source':str(source),'coefficient_rows':data['projection'].shape[0],'multiplier_rank':len(data['basis']),
               'initial_rays':len(seeds),'quartic_dictionaries':len(blocks),'maximum_pricing_matrix':max(len(b['words']) for b in blocks),
               'iterations_budget':iterations,'pool_size':pool,'maximum_weight':100,'build_seconds':time.monotonic()-started,
               'backend':backend,'matrix_threshold':matrix_threshold,'resumed_cuts':len(cuts),'center_dual':center}
    (out/'structure.json').write_text(json.dumps(structure,indent=2)+'\n');print(json.dumps(structure),flush=True)
    history=[];best=None
    for iteration in range(iterations+1):
        result,status=solve_weights(data,seeds,columns,backend=backend,matrix_threshold=matrix_threshold)
        step=dict(status,iteration=iteration,added_rays=len(cuts))
        if not result.success:
            history.append(step);break
        # HiGHS reports derivatives of the minimized -b objective with respect
        # to rhs. The physical moment functional has the opposite sign.
        dual=-np.asarray(result.eqlin.marginals)
        if center:
            nrays=len(seeds);offset=1+len(data['basis'])
            capped=bool(np.max(result.x[offset:offset+nrays])>=99.999)
            centered,center_status=center_dual(data,columns,status['numerical_objective'],capped=capped)
            step['center']=center_status
            if centered is not None:dual=centered
        proposals,violated=price(blocks,maps,dual,pool)
        step.update(dual_normalization=float(data['unit']@dual),violated_dictionaries=violated,
                    minimum_priced_eigenvalue=proposals[0]['violation'] if proposals else 0)
        if iteration%4==0 or iteration==iterations or not proposals:
            certificate,receipt=export_weights(data,seeds,result)
            receipt.update(iteration=iteration,added_rays=len(cuts),numerical_objective=status['numerical_objective'])
            checkpoint=out/f'iteration_{iteration}';checkpoint.mkdir(exist_ok=True)
            (checkpoint/'certificate.json').write_text(json.dumps(certificate)+'\n');(checkpoint/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
            step['certified_width']=receipt['width_float']
            if best is None or F(receipt['lower'])>F(best['lower']):
                best=receipt
                (out/'certificate.json').write_text(json.dumps(certificate)+'\n');(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        history.append(step);(out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
        print(json.dumps({k:step[k] for k in ('iteration','added_rays','numerical_objective','minimum_priced_eigenvalue','lp_seconds')})+(f" width={step['certified_width']}" if 'certified_width' in step else ''),flush=True)
        if iteration==iterations or not proposals:break
        for p in proposals:
            seeds.append(('exchange_squares',p['polynomial']));columns.append(p['column'])
            cuts.append({'dictionary':p['dictionary'],'polynomial':encode(p['polynomial']),'violation':p['violation']})
        (out/'cuts.json').write_text(json.dumps(cuts)+'\n')
    (out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
    summary={'best_width':best['width_float'] if best else None,'best_iteration':best['iteration'] if best else None,
             'iterations_completed':history[-1]['iteration'],'total_seconds':time.monotonic()-started,
             'stop':'iteration_budget' if iteration==iterations else ('no_detected_violation' if result.success else 'lp_failure')}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
    return summary


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=ROOT/'results/marginal_asymmetric_adapt/1_1000_penalty/polished/certificate.json');parser.add_argument('--iterations',type=int,default=16);parser.add_argument('--pool',type=int,default=8);parser.add_argument('--name',default='epsilon001')
    parser.add_argument('--backend',choices=['scipy','native'],default='native');parser.add_argument('--matrix-threshold',type=float,default=1e-12);parser.add_argument('--resume',type=Path);parser.add_argument('--center',action='store_true')
    args=parser.parse_args();run(args.source,args.iterations,args.pool,args.name,args.backend,args.matrix_threshold,args.resume,args.center)
