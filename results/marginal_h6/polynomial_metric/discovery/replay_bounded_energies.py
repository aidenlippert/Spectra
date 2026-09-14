from pathlib import Path
import json,subprocess,sys
root=Path('results/marginal_h6/polynomial_metric');base=root/'bounded_quotient'
for folder,source in [('proof','vanishing_metric_proof'),('transfer_proof','vanishing_metric_transfer/proof')]:
 out=base/folder;gap=json.loads((out/'certificate.json').read_text());c=json.loads((root/source/'energy_certificate.json').read_text());r=c['spin_symmetric_certificate'];assert r['hamiltonian']==gap['hamiltonian'];r['complement_lower']=gap['target_lower'];r.pop('complement_reference',None);r['complement_atoms']={k:gap[k] for k in ['kind','polynomial_metric','weight_proof','numerator_proof']};(out/'energy_certificate.json').write_text(json.dumps(c,indent=2)+'\n')
 for module,name,result in [('experiments.marginal_polynomial_metric','certificate','independent_replay'),('experiments.marginal_spin_temple','energy_certificate','energy_independent_replay')]:
  with (out/(result+'.json')).open('w') as f:subprocess.run([sys.executable,'-S','-m',module,'--verify',str(out/(name+'.json'))],stdout=f,check=True)
 assert json.loads((out/'receipt.json').read_text())==json.loads((out/'independent_replay.json').read_text());print(folder,'independent gap and energy passed',flush=True)
