"""Exact finite charge-pattern validation, not an implicit atom certificate."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_joint_coefficient_constructor import vanishing_feature
from experiments.marginal_polynomial_metric import JointPolynomial

root = Path('results/marginal_graded_hubbard8'); diagnostic = json.loads((root/'charge_metric_diagnostic.json').read_text())
data = json.loads((root/'hamiltonian.json').read_text()); assert data['hamiltonian'] == build(8, 4, 1)['hamiltonian']
proposal = next(r for r in diagnostic['results'] if r['gamma'] == -14)
terms = {}
for orbit, coefficient in zip(diagnostic['feature_orbits'], proposal['metric_coefficients']):
    integer = round(coefficient*10**9)
    for term in vanishing_feature(orbit, 8)['terms']:
        powers = tuple(term['powers']); terms[powers] = terms.get(powers, 0)+integer*term['coefficient']
candidate = dict(data, target_lower='-14', polynomial_metric={'denominator': 2*10**9,
    'terms': [{'powers': list(p), 'coefficient': c} for p,c in sorted(terms.items()) if c]})
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
assert count==1106 and minimum_w>0 and minimum_k>0
(root/'quadratic_metric_candidate.json').write_text(json.dumps(candidate,indent=2)+'\n')
receipt={'gamma':'-14','minimum_metric':str(minimum_w),'minimum_numerator':str(minimum_k),'minimum_numerator_float':float(minimum_k),'charge_patterns_checked':count,'worst_representative_state':worst,'actual_hamiltonian_compilation':cost,'scope':'Exact finite validation for the checked uniform Hubbard chain, enumerating1106 charge patterns and one worst-spin representative each. Alternating every single-occupation run maximizes positive target-metric hopping contributions; even total singles lets odd runs be balanced to Sz=0. All target charge metrics are positive on Q and zero on P. This establishes a physical weighted-row bound for this finite candidate, but no degree-six atom decomposition or implicit certificate has been constructed. No scalable enumeration claim.'}
(root/'quadratic_metric_finite_validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
