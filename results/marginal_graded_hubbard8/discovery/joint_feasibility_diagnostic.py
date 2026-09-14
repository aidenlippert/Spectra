"""Exact membership and normalization of the imported accepted metric.

This diagnoses the joint search; it is not fresh metric discovery.
"""
from pathlib import Path
from fractions import Fraction as F
import json,time
from math import comb
from experiments.marginal_polynomial_metric import JointPolynomial
from experiments.marginal_joint_coefficient_constructor import feature_orbits,vanishing_feature,CoefficientQuotient,sector_mean
from experiments.marginal_number_quotient import _add

root=Path('results/marginal_graded_hubbard8');candidate=json.loads((root/'quadratic_metric_candidate.json').read_text());bare={k:candidate[k] for k in ('modes','particles','hamiltonian')};ring=JointPolynomial(candidate);start=time.monotonic()
w,k,cost=ring.compile(F(candidate['target_lower']));target={m:F(v,cost['metric_scale']) for m,v in w.items()};q=CoefficientQuotient(ring,4);pivots={};orbits=[]
charges=[{0:-1,1<<(2*i):1,1<<(2*i+1):1} for i in range(ring.sites)]
for orbit in feature_orbits(ring.sites,2):
    if any(sum(p)%2 for p in orbit):continue
    trial=JointPolynomial(dict(bare,polynomial_metric=vanishing_feature(orbit,ring.sites)))
    polynomial=q.normal({m:F(v,2) for m,v in trial.metric_polynomial(charges).items()})
    j=len(orbits);expression={j:F(1)}
    while polynomial:
        pivot=min(polynomial,key=q.index.__getitem__);value=polynomial[pivot]
        if pivot not in pivots:
            pivots[pivot]=({m:v/value for m,v in polynomial.items()},{i:v/value for i,v in expression.items()});orbits.append(orbit);break
        previous,coordinates=pivots[pivot];_add(polynomial,previous,-value);_add(expression,coordinates,-value)
remainder=q.normal(target);coordinates={}
while remainder:
    pivot=min(remainder,key=q.index.__getitem__);value=remainder[pivot]
    if pivot not in pivots:raise ValueError('Metric lies outside even quadratic family')
    previous,expression=pivots[pivot];_add(remainder,previous,-value);_add(coordinates,expression,value)
# Independent synthesis of the metric and its compiled numerator in this basis.
w_sum={};k_sum={};q6=CoefficientQuotient(ring,6)
for j,value in coordinates.items():
    trial=JointPolynomial(dict(bare,polynomial_metric=vanishing_feature(orbits[j],ring.sites)))
    wi,ki,ci=trial.compile(F(candidate['target_lower']));_add(w_sum,{m:F(v,ci['metric_scale']) for m,v in wi.items()},value);_add(k_sum,{m:F(v,ci['numerator_scale']) for m,v in ki.items()},value)
if q.normal(w_sum)!=q.normal(target):raise ValueError('Metric synthesis failed')
if q6.normal(k_sum)!=q6.normal({m:F(v,cost['numerator_scale']) for m,v in k.items()}):raise ValueError('Numerator synthesis failed')
dimension=comb(ring.sites,ring.target);mean=F(dimension,dimension-1)*sector_mean(target,ring.sites,ring.target)
receipt=json.loads((root/'fixed_reynolds/receipt.json').read_text());lower={part:F(receipt[part+'_positivity']['lower'])/mean for part in ('weight','numerator')};floors={'weight':F(9,10000),'numerator':F(1,10000)}
if mean<=0 or any(lower[p]<=floors[p] for p in floors):raise ValueError('Normalized floors exceed accepted certificate guarantee')
proof=json.loads((root/'fixed_reynolds/certificate.json').read_text())
for part in ('weight','numerator'):
    item=proof[part+'_proof']
    if any(a['required'].bit_count()>6 for a in item['positive_indicators']):raise ValueError('Residual absorption degree unsupported')
    if any(a['required'].bit_count()>4 for a in item['charge_indicators']):raise ValueError('Charge degree unsupported')
    if any(a['mask'].bit_count()>5 for terms in item['number_multipliers'] for a in terms):raise ValueError('Number degree unsupported')
result={'gamma':candidate['target_lower'],'independent_even_features':len(orbits),'coordinates':[{'orbit':orbits[i],'coefficient':str(v)} for i,v in sorted(coordinates.items())],'mean_Q':str(mean),'normalized_weight_lower':str(lower['weight']),'normalized_numerator_lower':str(lower['numerator']),'required_weight_floor':str(floors['weight']),'required_numerator_floor':str(floors['numerator']),'residual_absorption':{'maximum_degree':6,'identity':'r+||r||_1=sum_(c>0) c*(1+n_M)+sum_(c<0) |c|*(1-n_M)','complement_identity':'1-product n_i=sum_j (1-n_j)*product_(i<j)n_i','admissibility':'Any nonzero degree<=6 indicator on two eight-mode four-particle slices has multiple completions; a singleton requires at least4 specifications in each spin, hence at least8 total. Zero-sector atoms are discarded. Group averaging preserves the exact invariant target.','conclusion':'The normalized complete degree-six joint LP has a feasible point; timeouts are search failures, not a cone obstruction.'},'metric_quotient_identity':True,'numerator_quotient_identity':True,'seconds':time.monotonic()-start,'scope':'Exact algebraic feasibility diagnostic for the normalized even quadratic joint family. Imports the already accepted metric and proof, so this does not constitute fresh joint discovery. Excess certified margins can be represented by the constant positive atom.'}
(root/'joint_feasibility_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
