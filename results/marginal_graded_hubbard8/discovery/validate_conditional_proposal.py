"""Finite validation of the last fresh conditional-search metric proposal.

This diagnostic enumerates charge patterns and is not part of discovery.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_joint_coefficient_constructor import prepare_features,vanishing_feature,input_digest
from experiments.marginal_polynomial_metric import JointPolynomial
root=Path('results/marginal_graded_hubbard8');data=json.loads((root/'hamiltonian.json').read_text())
assert data['hamiltonian']==build(8,4,1)['hamiltonian']
receipt=json.loads((root/'joint_reynolds_conditional/receipt.json').read_text())
assert receipt['source_sha256']==input_digest(data)
history=json.loads((root/'joint_reynolds_conditional/history.json').read_text())
proposal=next(r for r in reversed(history) if 'metric_coefficients' in r)
_,_,all_orbits,_,_,_=prepare_features(data,F(-14),feature_degree=2,quotient_degree=6)
orbits=[orbit for orbit in all_orbits if all(sum(p)%2==0 for p in orbit)]
assert len(orbits)==len(proposal['metric_coefficients'])
terms={}
for orbit,coefficient in zip(orbits,proposal['metric_coefficients']):
    integer=round(coefficient*10**9)
    for term in vanishing_feature(orbit,8)['terms']:
        powers=tuple(term['powers']);terms[powers]=terms.get(powers,0)+integer*term['coefficient']
candidate=dict(data,target_lower='-14',polynomial_metric={'denominator':2*10**9,'terms':[{'powers':list(p),'coefficient':c} for p,c in sorted(terms.items()) if c]})
ring = JointPolynomial(candidate); w,k,cost = ring.compile(F(-14)); minimum_w = minimum_k = None; count = 0; worst = None
for q in product((-1,0,1), repeat=8):
    if sum(q) or not any(q):continue
    symbols = [3 if x==1 else 0 for x in q]; runs=[]; start=0
    while start<8:
        if q[start]:start+=1;continue
        end=start
        while end<8 and q[end]==0:end+=1
        runs.append(list(range(start,end)));start=end
    odd = sum(len(run)%2 for run in runs); assert odd%2==0; assigned=0
    for run in runs:
        initial = 1
        if len(run)%2:
            initial = 1 if assigned<odd//2 else 2; assigned+=1
        for j,i in enumerate(run):symbols[i]=initial if j%2==0 else 3-initial
    state=sum(symbol << (2*i) for i,symbol in enumerate(symbols))
    assert all((state&spin).bit_count()==4 for spin in ring.spin_masks)
    evaluate=lambda p:sum(v for m,v in p.items() if state&m==m)
    weight=F(evaluate(w),cost['metric_scale']); numerator=F(evaluate(k),cost['numerator_scale'])
    minimum_w = weight if minimum_w is None else min(minimum_w,weight)
    if minimum_k is None or numerator<minimum_k:minimum_k=numerator;worst=state
    count+=1
assert count==1106
out=root/'conditional_metric_diagnostic';out.mkdir(exist_ok=True)
(out/'candidate.json').write_text(json.dumps(candidate,indent=2)+'\n')
result={'source_receipt':'joint_reynolds_conditional/receipt.json','source_round':proposal['round'],'source_sha256':input_digest(data),'gamma':'-14','metric_positive':minimum_w>0,'minimum_metric':str(minimum_w),'minimum_numerator':str(minimum_k),'minimum_numerator_float':float(minimum_k),'worst_representative_state':worst,'charge_patterns_checked':count,'physical_row_pass':minimum_w>0 and minimum_k>=0,'scope':'Exact finite charge-pattern diagnostic of a fresh algebraic-search proposal; enumerates1106 patterns after discovery. With positive charge metric, alternating balanced single-spin runs realize worst hopping rows. No atom certificate or scalable enumeration claim.'}
(out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
