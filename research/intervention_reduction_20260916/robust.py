"""Exact control-neighborhood certificate, with fresh nominal trajectory replay."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
from experiments.marginal_symbolic import decode, hermitian
from research.intervention_reduction_20260916.exact import integer_action,rational
from research.intervention_reduction_20260916.trajectory import check_trajectory


def control_norm_bound(encoded,modes):
    poly=decode(encoded,modes,4)
    if not hermitian(poly) or any(sum(2*c-1 for c,_ in w) for w in poly):
        raise ValueError('Hermitian number-conserving control required')
    support=sorted({p for w in poly for _,p in w})
    if len(support)>8:
        return sum(map(abs,poly.values()),F(0)),{'method':'CAR coefficient one-norm','support_modes':len(support)}
    order={p:i for i,p in enumerate(support)}
    local={tuple((c,order[p]) for c,p in w):v for w,v in poly.items()}
    # An even local CAR operator is unitarily equivalent, after permuting modes,
    # to this local Fock operator tensor identity on the complementary modes.
    den,columns=integer_action(local,range(1<<len(support)))
    rows=[0]*(1<<len(support))
    for col in columns:
        for i,value in col.items():
            rows[i]+=abs(value)
    return F(max(rows),den),{'method':'local Fock maximum absolute row sum',
        'support_modes':len(support),'local_dimension':len(rows)}


def extend_control_neighborhood(data,proposal,trajectory,nominal,l1_budgets):
    if len(l1_budgets)!=len(proposal['controls']) or any(x<0 for x in l1_budgets):
        raise ValueError('One nonnegative integrated control deviation per control')
    norms=[control_norm_bound(c['operator'],data['modes']) for c in proposal['controls']]
    allowance=sum(b*n[0] for b,n in zip(l1_budgets,norms))
    error=F(nominal['normalized_state_error_bound'])+allowance
    uniform=nominal.get('uniform_in_time_state_error_bound')
    uniform=None if uniform is None else F(uniform)+allowance
    expected=F(nominal['predicted_population']);initial=F(nominal['initial_population'])
    lo=max(F(0),expected-2*error);hi=min(F(2),expected+2*error)
    return {**nominal,'status':'accepted_control_neighborhood',
        'nominal_state_error_bound':nominal['normalized_state_error_bound'],
        'normalized_state_error_bound':str(error),'normalized_state_error_float':float(error),
        'target_met':error<=F(1,200),'control_deviation_state_allowance':str(allowance),
        'uniform_in_time_state_error_bound':None if uniform is None else str(uniform),
        'uniform_in_time_state_error_float':None if uniform is None else float(uniform),
        'uniform_in_time_target_met':uniform is not None and uniform<=F(1,200),
        'integrated_absolute_control_deviations_Ha_au':[str(x) for x in l1_budgets],
        'control_norms':[{'bound':str(n),**meta} for n,meta in norms],
        'population_interval':[str(lo),str(hi)],'population_interval_float':[float(lo),float(hi)],
        'population_change_interval':[str(lo-initial),str(hi-initial)],
        'population_change_interval_float':[float(lo-initial),float(hi-initial)],
        'control_scope':'all measurable perturbations satisfying the integrated absolute deviation budgets and the original amplitude box; same initial state',
        'proof':'unitary Duhamel: state deviation <= sum_l ||H_l|| integral |delta u_l(t)| dt'}


def check_robust(data,proposal,trajectory,l1_budgets):
    nominal=check_trajectory(data,proposal,trajectory)
    return extend_control_neighborhood(data,proposal,trajectory,nominal,l1_budgets)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture',type=Path);p.add_argument('proposal',type=Path)
    p.add_argument('trajectory',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--max-deviation-Ha',default='1/400000');a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    data,proposal,trajectory=[json.loads(x.read_text()) for x in (a.fixture,a.proposal,a.trajectory)]
    budget=rational(a.max_deviation_Ha)*rational(trajectory['horizon_atomic_time'])
    result=check_robust(data,proposal,trajectory,[budget]*len(proposal['controls']))
    a.output.write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps({k:result[k] for k in ('target_met','normalized_state_error_float','population_change_interval_float','control_norms')}),flush=True)
