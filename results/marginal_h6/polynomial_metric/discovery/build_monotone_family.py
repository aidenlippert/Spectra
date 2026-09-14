import json,subprocess,sys,time
from pathlib import Path
from fractions import Fraction as F
from unittest.mock import patch
from experiments.marginal_symbolic import decode,encode,add,scale
from experiments.marginal_polynomial_metric import JointPolynomial
from experiments.marginal_monotone_transfer import replay
root=Path('results/marginal_h6/polynomial_metric');base=root/'bounded_quotient';reference=json.loads((base/'proof/certificate.json').read_text());unit=json.loads((root/'vanishing_metric_transfer/candidate.json').read_text());h0=decode(reference['hamiltonian'],12,4);delta=scale(add(decode(unit['hamiltonian'],12,4),scale(h0,-1)),50);diagnostic=json.loads((base/'amplitude_transfer_diagnostic.json').read_text());limit=min(-2*F(g['old_range'][1]) for g in diagnostic['groups']);summary=[]
for name,strength in [('monotone_1x',F(1,50)),('monotone_3x',F(3,50)),('monotone_mixed_sign',F(3,20)),('monotone_limit',limit)]:
 c={'kind':'joint_polynomial_monotone_transfer_v1','modes':12,'particles':6,'hamiltonian':encode(add(h0,scale(delta,strength))),'target_lower':reference['target_lower'],'reference_certificate':reference};start=time.monotonic()
 with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No physical states')),patch('experiments.marginal_spin_constructor.spin_states',side_effect=AssertionError('No enumeration')):r=replay(c)
 direct_sign_gate=None
 if name=='monotone_mixed_sign':
  try:JointPolynomial(dict(c,polynomial_metric=reference['polynomial_metric'])).compile(F(c['target_lower']))
  except ValueError as e:
   if 'fixed-sign' not in str(e):raise
   direct_sign_gate=str(e)
  else:raise ValueError('Expected mixed-sign direct compiler refusal')
 out=base/name;out.mkdir(exist_ok=True)
 for filename,data in [('certificate',c),('receipt',r)]:(out/(filename+'.json')).write_text(json.dumps(data,indent=2)+'\n')
 with (out/'independent_replay.json').open('w') as f:subprocess.run([sys.executable,'-S','-m','experiments.marginal_monotone_transfer','--verify',str(out/'certificate.json')],stdout=f,check=True)
 assert r==json.loads((out/'independent_replay.json').read_text());summary.append({'name':name,'hopping_strength':str(strength),'hopping_strength_float':float(strength),'gap':r['complement_lower'],'changed_groups':len(r['changed_groups']),'certifies_entire_reference_target_segment':r['certifies_reference_target_segment'],'determinant_actions':r['determinant_actions'],'direct_compiler_refusal':direct_sign_gate,'elapsed_with_independent_replay':time.monotonic()-start});print(json.dumps(summary[-1]),flush=True);(base/'monotone_family.json').write_text(json.dumps({'scope':'Exact uniform Q-gap only. Existing metric/positive directions inherited; reference replay remains finite coefficient algebra. No new energy witness or response for the larger perturbations. Interval endpoint is sufficient, not asserted optimal.','rows':summary},indent=2)+'\n')
 if name=='monotone_1x':
  energy=json.loads((root/'vanishing_metric_transfer/proof/energy_certificate.json').read_text());red=energy['spin_symmetric_certificate'];assert red['hamiltonian']==c['hamiltonian'];red['complement_atoms']={'kind':c['kind'],'reference_certificate':reference};red['complement_lower']=c['target_lower'];(out/'energy_certificate.json').write_text(json.dumps(energy,indent=2)+'\n')
  with (out/'energy_independent_replay.json').open('w') as f:subprocess.run([sys.executable,'-S','-m','experiments.marginal_spin_temple','--verify',str(out/'energy_certificate.json')],stdout=f,check=True)
