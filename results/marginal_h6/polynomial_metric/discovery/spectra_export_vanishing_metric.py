import json,time
from pathlib import Path
from fractions import Fraction as F
from math import comb
from experiments.marginal_polynomial_metric import JointPolynomial,complete_number_ideals,replay
root=Path('/Users/aidenlippert/Documents/Spectra/results/marginal_h6/polynomial_metric');out=root/'vanishing_metric_proof';out.mkdir(exist_ok=True);proposal=json.loads((root/'vanishing_joint_metric_cone/proposal_-6.264.json').read_text());c=json.loads((root/'candidate.json').read_text());metric_den=2*10**9;proof_den=10**12
metric_terms={}
for orbit,value in zip(proposal['orbits'],proposal['metric_coefficients']):
 integer=round(value*10**9)
 for powers in orbit:
  for site in range(6):
   expanded=list(powers);expanded[site]+=2;expanded=tuple(1 if x%2 else 2 if x else 0 for x in expanded);metric_terms[expanded]=metric_terms.get(expanded,0)+integer
c['polynomial_metric']={'denominator':metric_den,'terms':[{'powers':list(powers),'coefficient':value} for powers,value in metric_terms.items() if value]};c.update(kind='joint_polynomial_metric_gap_v1',target_lower=proposal['gamma']);ring=JointPolynomial(c);start=time.monotonic();weight,numerator,cost=ring.compile(F(c['target_lower']));construction={'compile_seconds':time.monotonic()-start,'compilation':cost,'discovery':'Joint metric/positivity LP on400 explicit spin-sector coordinates; exact number-ideal lift may have sector-sized cost.'}
for name,poly,scale in [('weight',weight,cost['metric_scale']),('numerator',numerator,cost['numerator_scale'])]:
 proof={'denominator':proof_den,'bound':round(proposal['metric_bound' if name=='weight' else 'numerator_bound']*proof_den),'positive_indicators':[],'charge_indicators':[],'number_multipliers':[[],[]]};represented={0:proof['bound']};localizer={0:-1,**{3<<(2*i):1 for i in range(ring.sites)}};negative=[];counts=[]
 for label,value in zip(proposal[name+'_labels'],proposal[name+'_values']):
  if value<0:negative.append(value)
  integer=max(0,round(value*proof_den))
  if not integer:continue
  family,required,occupied=label;atom={'required':required,'occupied':occupied,'weight':integer};proof['positive_indicators' if family=='positive' else 'charge_indicators'].append(atom);term={occupied:integer}
  for i in range(ring.modes):
   if (required^occupied)&(1<<i):term=ring.multiply(term,{0:1,1<<i:-1})
  if family=='charge':term=ring.multiply(term,localizer)
  represented=ring.add(represented,term)
  if family=='positive':
   count=1
   for spin in ring.spin_masks:
    available=ring.sites-(required&spin).bit_count();needed=ring.target-(occupied&spin).bit_count();count*=comb(available,needed) if 0<=needed<=available else 0
   if count<=1:raise ValueError('Singleton or empty positive atom')
   counts.append(count)
 residual=ring.add({m:F(v,scale) for m,v in poly.items()},{m:-F(v,proof_den) for m,v in represented.items()});remainder,ideals=complete_number_ideals(ring,residual)
 for i,ideal in enumerate(ideals):
  proof['number_multipliers'][i]=[{'mask':m,'coefficient':round(v*proof_den)} for m,v in ideal.items() if round(v*proof_den)]
 c[name+'_proof']=proof;construction[name]={'quotient_remainder_l1':str(sum(map(abs,remainder.values()),F(0))),'negative_proposal_count':len(negative),'most_negative_proposal':min(negative,default=0),'positive_atom_minimum_completions':min(counts),'positive_atoms':len(counts),'number_multiplier_terms':sum(map(len,proof['number_multipliers']))}
 print(name,construction[name],flush=True)
receipt=replay(c)
for name,value in [('certificate',c),('receipt',receipt),('construction',construction)]: (out/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
print(json.dumps(receipt,indent=2),flush=True)
