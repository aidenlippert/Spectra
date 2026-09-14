from pathlib import Path
from fractions import Fraction as F
from unittest.mock import patch
from math import comb
import json,sys,subprocess,time
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_polynomial_metric import replay
base=Path('results/marginal_hubbard_polynomial');base.mkdir(exist_ok=True);rows=[]
for sites in [4,6,8,10,16,24,32]:
 out=base/str(sites);out.mkdir(exist_ok=True);start=time.monotonic()
 with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No states permitted')),patch('experiments.marginal_polynomial_metric.complete_number_ideals',side_effect=AssertionError('No sector-sized lift permitted')),patch('experiments.marginal_spin_constructor.spin_states',side_effect=AssertionError('No enumeration permitted')):
  c=build(sites,4,F(1,3));built=time.monotonic();r=replay(c);end=time.monotonic()
 assert r['weight_positivity']['residual_l1']==r['numerator_positivity']['residual_l1']=='0'
 for name,obj in [('certificate',c),('receipt',r)]:(out/(name+'.json')).write_text(json.dumps(obj,indent=2)+'\n')
 with (out/'independent_replay.json').open('w') as f:subprocess.run([sys.executable,'-S','-m','experiments.marginal_polynomial_metric','--verify',str(out/'certificate.json')],stdout=f,check=True)
 assert r==json.loads((out/'independent_replay.json').read_text())
 row={'sites':sites,'spin_sector_dimension':comb(sites,sites//2)**2,'certificate_bytes':(out/'certificate.json').stat().st_size,'construction_seconds':built-start,'replay_seconds':end-built,'gamma':r['complement_lower'],'positive_atoms':len(c['numerator_proof']['positive_indicators']),'charge_atoms':len(c['numerator_proof']['charge_indicators']),'number_multiplier_terms':sum(map(len,c['numerator_proof']['number_multipliers'])),'compilation':r['compilation'],'verification_coefficient_products':r['verification_coefficient_products'],'maximum_occupation_degree':max(a['required'].bit_count() for a in c['numerator_proof']['positive_indicators']),'maximum_number_multiplier_degree':max(a['mask'].bit_count() for v in c['numerator_proof']['number_multipliers'] for a in v)};rows.append(row);print(json.dumps(row),flush=True)
 (base/'scaling.json').write_text(json.dumps({'scope':'Exact Q-gap certificates only for Hubbard chains U4 t1/3 slack1/1000; no energy interval, response compression or molecular transfer claimed. Constructor and replay refuse determinant actions and sector lifting during measurement.','rows':rows},indent=2)+'\n')
