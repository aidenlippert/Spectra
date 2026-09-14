"""Exact physical-row obstructions to two restricted positive metric families."""
from fractions import Fraction as F
import json
from pathlib import Path
from experiments.marginal_spin_reduction import SpinZeroOracle


def metric_class(oracle,state,family):
    m=oracle.modes//2
    if not oracle.valid_state(state):raise ValueError('Physical fixed-spin source required')
    charges=[((state>>(2*i))&3).bit_count()-1 for i in range(m)]
    d=charges.count(1)
    if not d:raise ValueError('Ionic source required')
    if family=='doublon_only':return d-1,m//2
    if family=='one_pair_constant_tail':
        if d>1:return 0,1+m*(m-1)
        i,j=charges.index(1),charges.index(-1)
        return 1+i*(m-1)+(j if j<i else j-1),1+m*(m-1)
    raise ValueError('Unsupported restricted metric family')


def row_coefficients(oracle,state,family,gamma):
    label,count=metric_class(oracle,state,family);row=[F(0)]*count
    action=oracle.action(state);row[label]=action.get(state,F(0))-gamma
    for target,value in action.items():
        if target==state or not any(((target>>(2*i))&3)==3 for i in range(oracle.modes//2)):continue
        group,_=metric_class(oracle,target,family);row[group]-=abs(value)
    return row


def replay(certificate):
    if certificate.get('kind')!='charge_metric_family_obstruction_v1':raise ValueError('Unsupported metric-family obstruction')
    oracle=SpinZeroOracle(certificate)
    if oracle.modes%4 or oracle.particles!=oracle.modes//2:raise ValueError('Even half-filled model required')
    if type(certificate.get('target_lower')) is not str:raise ValueError('Exact target required')
    family=certificate['metric_family'];gamma=F(certificate['target_lower']);raw=certificate['row_weights']
    if type(raw) is not list or not 1<=len(raw)<=512:raise ValueError('Bounded nonempty source-row weights required')
    total=None;seen=set();mass=F(0)
    for item in raw:
        state=item['state']
        if state in seen or type(item['weight']) is not str or F(item['weight'])<=0:raise ValueError('Distinct physical rows and positive exact weights required')
        coefficients=row_coefficients(oracle,state,family,gamma);weight=F(item['weight']);mass+=weight;seen.add(state)
        if total is None:total=[F(0)]*len(coefficients)
        total=[x+weight*y for x,y in zip(total,coefficients)]
    if any(x>=0 for x in total):raise ValueError('Strictly negative weighted sum in every metric class required')
    return {'target_lower':str(gamma),'metric_family':family,'metric_classes':len(total),'source_rows':len(raw),
            'row_weight_sum':str(mass),'class_coefficients':[str(x) for x in total],
            'least_negative_coefficient':str(max(total)),'determinant_actions':len(oracle.cache),
            'scope':'Exact impossibility of certifying this target by weighted Q-row dominance with any strictly positive weights in the stated metric family. A positive combination of actual Hamiltonian row inequalities has strictly negative coefficients in every class. Only witness rows are required. This is not an energy lower/upper bound or a no-go theorem for other metrics or positivity certificates.'}


def propose(source,output,family,target='-6.264',denominator=10**12):
    import numpy as np
    from scipy.optimize import linprog
    from experiments.marginal_spin_constructor import spin_states
    source,output=Path(source),Path(output)
    if output.exists():raise ValueError('Preserve previous metric obstruction')
    if type(denominator) is not int or denominator<1:raise ValueError('Positive rounding denominator required')
    data=json.loads(source.read_text());oracle=SpinZeroOracle(data)
    if oracle.sector_dimension>1000:raise ValueError('Discovery limited to1000 spin configurations')
    states=[s for s in spin_states(oracle) if any(((s>>(2*i))&3)==3 for i in range(oracle.modes//2))]
    a=np.array([[float(x) for x in row_coefficients(oracle,s,family,F(target))] for s in states])
    objective=np.r_[np.zeros(len(states)),1.0]
    result=linprog(objective,A_ub=np.column_stack([a.T,-np.ones(a.shape[1])]),b_ub=np.zeros(a.shape[1]),
        A_eq=[np.r_[np.ones(len(states)),0]],b_eq=[1],bounds=[(0,None)]*len(states)+[(None,None)],method='highs')
    if not result.success or result.x[-1]>=0:raise ValueError('No strictly negative numerical obstruction found')
    certificate={k:data[k] for k in ['modes','particles','hamiltonian']}
    certificate.update(kind='charge_metric_family_obstruction_v1',target_lower=str(F(target)),metric_family=family,
                       row_weights=[{'state':s,'weight':str(F(round(y*denominator),denominator))} for s,y in zip(states,result.x[:-1]) if round(y*denominator)>0])
    receipt=replay(certificate);output.mkdir(parents=True)
    proposal={'numerical_margin':float(result.x[-1]),'discovery_Q_states':len(states),'discovery_actions':len(oracle.cache),'scope':'Full finite-sector numerical proposal; independent exact replay binds selected witness rows to the physical Hamiltonian.'}
    for name,value in [('certificate',certificate),('receipt',receipt),('proposal',proposal)]: (output/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--verify');p.add_argument('--source');p.add_argument('--output');p.add_argument('--family');p.add_argument('--target',default='-6.264');args=p.parse_args()
    print(json.dumps(replay(json.loads(Path(args.verify).read_text())) if args.verify else propose(args.source,args.output,args.family,args.target),indent=2))
