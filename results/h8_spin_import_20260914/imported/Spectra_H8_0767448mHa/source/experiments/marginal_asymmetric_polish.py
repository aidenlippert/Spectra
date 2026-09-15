"""Reweight asymmetric square directions with an LP and exact replay."""
from fractions import Fraction as F
import json
from pathlib import Path
import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog,OptimizeResult

from experiments.marginal_asymmetric_adapt import prepare
from experiments.marginal_reynolds import reynolds_matrix
from experiments.marginal_symbolic import decode,encode,product,adj,add,scale,transform
from experiments.marginal_orbit_certificate import average_canonical_polynomial,normalize_scalar_residual
from experiments.marginal_symmetry_transfer import replay

ROOT=Path(__file__).resolve().parents[1]


def prepare_polish(source,average_direct=False):
    started=time.monotonic();source=Path(source);original=json.loads(source.read_text())
    data=prepare(ROOT/'results/marginal_distill_direct/compact_certificate.json')
    modes=original['modes'];particles=original['particles']
    if modes!=data['modes'] or particles!=data['particles'] or original['permutations']!=data['base']['permutations']:
        raise ValueError('Polishing expects the ten-mode matched orbit group and fixed sector')
    rows=data['rows'];projection=data['projection'];full_average=reynolds_matrix(rows,modes)[0]
    seeds=[];columns=[]
    for kind in ('orbit_squares','direct_squares'):
        for item in original.get(kind,[]):
            p=decode(item['polynomial'],modes,4);largest=max(abs(c) for c in p.values());p=scale(p,1/largest)
            square=product(adj(p),p);column=np.array([float(square.get(w,0)) for w in rows])
            if kind=='orbit_squares':column=full_average@column
            output_kind='exchange_squares' if average_direct and kind=='direct_squares' else kind
            seeds.append((output_kind,p));columns.append(projection@column)
    h=decode(original['hamiltonian'],modes,4)
    exchange=[(i+5)%modes for i in range(modes)]
    if transform(h,exchange)!=h:raise ValueError('Hamiltonian breaks the exchange used for coefficient projection')
    data.update(original=original,source=str(source),seeds=seeds,columns=columns,
                rhs=projection@np.array([float(h.get(w,0)) for w in rows]),polish_build_seconds=time.monotonic()-started)
    return data


def native_lp(objective,matrix,rhs,bounds,matrix_threshold):
    import highspy
    h=highspy.Highs()
    for name,value in [('output_flag',False),('threads',1),('small_matrix_value',matrix_threshold),
                       ('primal_feasibility_tolerance',1e-9),('dual_feasibility_tolerance',1e-9)]:
        if h.setOptionValue(name,value)!=highspy.HighsStatus.kOk:raise ValueError(f'HiGHS rejected {name}')
    lp=highspy.HighsLp();lp.num_col_=matrix.shape[1];lp.num_row_=matrix.shape[0]
    lp.col_cost_=objective
    lp.col_lower_=np.array([-np.inf if a is None else a for a,b in bounds]);lp.col_upper_=np.array([np.inf if b is None else b for a,b in bounds])
    lp.row_lower_=rhs;lp.row_upper_=rhs
    lp.a_matrix_.format_=highspy.MatrixFormat.kRowwise;lp.a_matrix_.num_col_=matrix.shape[1];lp.a_matrix_.num_row_=matrix.shape[0]
    lp.a_matrix_.start_=matrix.indptr;lp.a_matrix_.index_=matrix.indices;lp.a_matrix_.value_=matrix.data
    loaded=h.passModel(lp)
    if loaded==highspy.HighsStatus.kError:raise ValueError('HiGHS rejected LP model')
    h.run();status=h.getModelStatus();solution=h.getSolution()
    return OptimizeResult(success=status==highspy.HighsModelStatus.kOptimal,message=h.modelStatusToString(status),
                          x=np.asarray(solution.col_value),fun=h.getObjectiveValue(),
                          eqlin=OptimizeResult(marginals=np.asarray(solution.row_dual)),
                          highs_version=h.version(),matrix_threshold=h.getOptionValue('small_matrix_value')[1])


def solve_weights(data,seeds,columns,backend='scipy',matrix_threshold=1e-9,max_weight=100):
    if backend not in ('scipy','native'):raise ValueError('Invalid LP backend')
    if not np.isfinite(max_weight) or max_weight<=0:raise ValueError('Positive finite square weight cap required')
    if backend=='scipy' and matrix_threshold!=1e-9:raise ValueError('Use native backend for matrix threshold control')
    n=data['projection'].shape[0];rank=len(data['basis']);nrays=len(seeds)
    matrix=sp.hstack([sp.csr_matrix(data['unit'][:,None]),sp.csr_matrix(data['a']),sp.csr_matrix(np.column_stack(columns)),sp.eye(n),-sp.eye(n)],format='csr')
    objective=np.r_[-1,np.zeros(rank+nrays),data['weights'],data['weights']]
    residual_bounds=[(0,0)]+[(0,None)]*(n-1)
    bounds=[(None,None)]*(1+rank)+[(0,max_weight)]*nrays+residual_bounds+residual_bounds
    started=time.monotonic()
    result=(native_lp(objective,matrix,data['rhs'],bounds,matrix_threshold) if backend=='native' else
            linprog(objective,A_eq=matrix,b_eq=data['rhs'],bounds=bounds,method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9}))
    status={'source':data['source'],'rows':n,'multiplier_rank':rank,'candidate_squares':nrays,
            'build_seconds':data['polish_build_seconds'],'lp_seconds':time.monotonic()-started,
            'status':result.message,'success':bool(result.success),'maximum_square_weight':max_weight,'backend':backend,
            'matrix_threshold':matrix_threshold,'highs_version':result.get('highs_version','bundled with SciPy')}
    if result.success:
        ray_weights=result.x[1+rank:1+rank+nrays]
        status['largest_square_weight']=float(np.max(ray_weights)) if nrays else 0.
        status['capped_square_weights']=int(np.count_nonzero(ray_weights>=max_weight-1e-6))
        status['positive_square_weights']=int(np.count_nonzero(ray_weights>1e-9))
        status['numerical_objective']=-float(result.fun)
        status['maximum_numerical_equation_residual']=float(np.max(abs(matrix@result.x-data['rhs'])))
        status['weighted_numerical_equation_residual']=float(data['weights']@abs(matrix@result.x-data['rhs']))
    return result,status


def export_weights(data,seeds,result):
    if not result.success:raise ValueError('Cannot export unsuccessful LP')
    modes=data['modes'];particles=data['particles'];original=data['original'];rank=len(data['basis']);nrays=len(seeds)
    def rational(c):return F(round(float(c)*10**12),10**12)
    selected={'orbit_squares':[],'direct_squares':[]};exchange=[(i+5)%modes for i in range(modes)]
    for (kind,p),weight in zip(seeds,result.x[1+rank:1+rank+nrays]):
        w=rational(max(0,weight));rounded={word:F(round(c*10**14),10**14) for word,c in p.items() if round(c*10**14)}
        if not w or not rounded:continue
        if kind=='exchange_squares':
            image=transform(rounded,exchange)
            if image==rounded or image==scale(rounded,-1):
                selected['direct_squares'].append({'polynomial':encode(rounded),'weight':str(w)})
            else:
                for poly in (rounded,image):
                    selected['direct_squares'].append({'polynomial':encode(poly),'weight':str(w/2)})
        else:selected[kind].append({'polynomial':encode(rounded),'weight':str(w)})
    multiplier=add(*(scale(p,rational(c)) for p,c in zip(data['basis'],result.x[1:1+rank])))
    multiplier=average_canonical_polynomial(multiplier,[list(range(modes)),exchange])
    certificate={'modes':modes,'particles':particles,'operator_degree':4,'hamiltonian':original['hamiltonian'],
                 'b':str(rational(result.x[0])),'number_multiplier':encode(multiplier),'permutations':original['permutations'],
                 'independent_upper':original['independent_upper'],**selected}
    certificate,_=normalize_scalar_residual(certificate)
    return certificate,replay(certificate)


def run(source,backend='native',matrix_threshold=1e-12):
    source=Path(source);data=prepare_polish(source);result,status=solve_weights(data,data['seeds'],data['columns'],backend,matrix_threshold)
    out=source.parent/'polished';out.mkdir(exist_ok=True)
    if not result.success:
        (out/'status.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status),flush=True);return status
    certificate,receipt=export_weights(data,data['seeds'],result);receipt.update(status)
    (out/'certificate.json').write_text(json.dumps(certificate)+'\n');(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('width_float','lower_float','orbit_squares','direct_squares','residual_l1_float','lp_seconds')}),flush=True)
    return receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('source',type=Path)
    parser.add_argument('--backend',choices=['scipy','native'],default='native');parser.add_argument('--matrix-threshold',type=float,default=1e-12)
    args=parser.parse_args();run(args.source,args.backend,args.matrix_threshold)
