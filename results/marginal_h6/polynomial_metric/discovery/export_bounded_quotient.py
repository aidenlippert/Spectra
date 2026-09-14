"""Rebuild number identities from coefficient algebra, preserving the supplied atoms."""
import argparse,json,time
from pathlib import Path
from fractions import Fraction as F
from unittest.mock import patch
from experiments.marginal_polynomial_metric import JointPolynomial,replay
from experiments.marginal_number_quotient import complete_bounded_number_ideals
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();c=json.loads(a.source.read_text());ring=JointPolynomial(c);start=time.monotonic()
with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No states')),patch('experiments.marginal_polynomial_metric.complete_number_ideals',side_effect=AssertionError('No full-population lift')),patch('experiments.marginal_spin_constructor.spin_states',side_effect=AssertionError('No sector generation')):
 weight,numerator,cost=ring.compile(F(c['target_lower']));construction={'source':str(a.source),'scope':'Inherited metric and positive atoms; number identities rebuilt in bounded-degree coefficient space, with no state enumeration or full-population lifting.','parts':{}}
 for name,poly,scale in [('weight',weight,cost['metric_scale']),('numerator',numerator,cost['numerator_scale'])]:
  proof=c[name+'_proof'];den=proof['denominator'];represented={0:F(proof['bound'],den)}
  for key in ['positive_indicators','charge_indicators']:
   for item in proof[key]:
    term={item['occupied']:F(item['weight'],den)}
    for i in range(ring.modes):
     if (item['required']^item['occupied'])&(1<<i):term=ring.multiply(term,{0:1,1<<i:-1})
    if key=='charge_indicators':term=ring.multiply(term,{0:-1,**{3<<(2*i):1 for i in range(ring.sites)}})
    represented=ring.add(represented,term)
  residual=ring.add({m:F(v,scale) for m,v in poly.items()},{m:-v for m,v in represented.items()});remainder,ideals,stats=complete_bounded_number_ideals(ring,residual);old=sum(map(len,proof['number_multipliers']));proof['number_multipliers']=[[{'mask':m,'coefficient':round(v*den)} for m,v in sorted(ideal.items()) if round(v*den)] for ideal in ideals];construction['parts'][name]={'quotient':stats,'remainder_terms':len(remainder),'remainder_l1':str(sum(map(abs,remainder.values()),F(0))),'old_number_multiplier_terms':old,'new_number_multiplier_terms':sum(map(len,proof['number_multipliers'])),'maximum_multiplier_degree':max((m.bit_count() for ideal in ideals for m in ideal),default=0)}
 r=replay(c)
construction['seconds']=time.monotonic()-start;a.out.mkdir(parents=True,exist_ok=True)
for name,value in [('certificate',c),('receipt',r),('construction',construction)]:(a.out/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
print(json.dumps({'weight':r['weight_positivity'],'numerator':r['numerator_positivity'],'construction':construction},indent=2))
