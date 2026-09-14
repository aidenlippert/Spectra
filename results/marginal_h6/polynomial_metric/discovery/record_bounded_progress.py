from pathlib import Path
from fractions import Fraction as F
import json,re,hashlib,shutil
ROOT=Path(__file__).resolve().parents[4];BASE=ROOT/'results/marginal_h6/polynomial_metric/bounded_quotient'
def read(p):return json.loads(p.read_text())
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
log=Path('/tmp/marginal_monotone_quotient_validation.log').read_text();match=re.search(r'Ran (\d+) tests in ([\d.]+)s\s+OK\s*$',log)
if not match or int(match.group(1))!=360:raise ValueError('Completed360-test validation required')
seconds=float(match.group(2));shutil.copyfile('/tmp/marginal_monotone_quotient_validation.log',BASE/'full_validation.log');checks=[];energies=[];proofs={}
for name in ['proof','transfer_proof','monotone_1x','monotone_3x','monotone_mixed_sign','monotone_limit']:
 folder=BASE/name;r=read(folder/'receipt.json');assert r==read(folder/'independent_replay.json');checks.append(str((folder/'independent_replay.json').relative_to(ROOT)));proofs[name]={'certificate':str((folder/'certificate.json').relative_to(ROOT)),'certificate_bytes':(folder/'certificate.json').stat().st_size,'receipt':r}
 if (folder/'energy_certificate.json').exists():
  e=read(folder/'energy_independent_replay.json');assert F(e['width'])>0;energies.append({'certificate':str((folder/'energy_certificate.json').relative_to(ROOT)),'independent':str((folder/'energy_independent_replay.json').relative_to(ROOT)),'width':e['width'],'width_float':float(F(e['width'])),'retained_dimension':e['retained_dimension'],'response_dimension':e['response_dimension'],'witness_support':e['witness_support'],'scope':'Existing witness and response reused and exactly rechecked.'})
assert len(checks)==6 and len(energies)==3
source=read(BASE/'proof/certificate.json');transfer=read(BASE/'transfer_proof/certificate.json');assert source['polynomial_metric']==transfer['polynomial_metric'];assert source['weight_proof']==transfer['weight_proof']
family=read(BASE/'monotone_family.json');assert len(family['rows'])==4;assert family['rows'][-1]['hopping_strength']=='14275927873493/50000000000000';assert family['rows'][-2]['direct_compiler_refusal']
progress={'status':'Exact bounded-degree identity reconstruction, fresh coefficient-space positivity transfer and uniform monotone hopping-family gap accepted; overall goal remains active.','reports':['research/marginal_bounded_quotient.md','research/marginal_monotone_transfer.md'],'previous_goal_turn':'progress','tests_passed':360,'elapsed_seconds':seconds,'gap_replays':checks,'energy_replays':energies,'new_stdlib_replay_count':9,'proofs':proofs,'monotone_family':family,'coefficient_space_lp':read(BASE/'transfer_lp/receipt.json'),'threefold_lp_proposal':read(BASE/'transfer_3x/lp/receipt.json'),'threefold_scope':'The3x numerical numerator proposal was not exported as its own certificate; an independently accepted monotone transfer proves this target gap.','remaining':['Metric and positive seed directions still inherit explicit finite-sector discovery.','H6 degree reaches the full slice quotient:400 coordinates and1364 multiplier terms remain.','Useful molecular accuracy at a controlled positivity degree and construction cost remains unproved.','Retained-response and upper-witness construction still use explicit configurations.','Larger-system molecular transfer, general representability, thermodynamics, kinetics and synthesis remain open.']}
save(BASE/'progress.json',progress)
vpath=ROOT/'results/marginal_final_validation.json';v=read(vpath);v.update(tests_passed=360,elapsed_seconds=seconds,latest_validation_scope='360 full marginal tests; nine new standard-library replays: two bounded-quotient gap proofs, four monotone transfer gaps, and three energies using existing witnesses/responses. Uniform H6 hopping-family Q-gap certified; no new larger-perturbation energy intervals or general representability claim.',bounded_quotient_progress=str((BASE/'progress.json').relative_to(ROOT)),bounded_quotient_reports=progress['reports'],bounded_quotient_stdlib_replays=checks+[e['independent'] for e in energies])
for e in energies:
 if e['certificate'] not in v['stdlib_replays_passed']:v['stdlib_replays_passed'].append(e['certificate'])
assert len(v['stdlib_replays_passed'])==145;save(vpath,v)
jpath=ROOT/'results/marginal_h6/polynomial_metric/joint_progress.json';j=read(jpath);j['bounded_quotient_followup']=str((BASE/'progress.json').relative_to(ROOT));save(jpath,j)
intro='The latest [bounded-degree coefficient quotient](marginal_bounded_quotient.md) rebuilds the accurate H6 proof and discovers fresh perturbed-H positivity coefficients without physical-configuration evaluation or full-population lifting. An [exact monotone-transfer theorem](marginal_monotone_transfer.md) preserves the −6.264-Ha complement bound throughout the specified hopping interval0≤δ≤0.28551855746986, including mixed-sign targets refused by direct compilation. **360 full tests and nine new independent replays pass.** The inherited metric/direction discovery,400-coordinate H6 quotient, and explicit energy-response work remain uncompressed.\n\n'
for name in ['research/marginal_blocker_progress.md','research/marginal_joint_metric.md']:
 p=ROOT/name;s=p.read_text();lines=s.splitlines(True)
 if 'The latest [bounded-degree coefficient quotient]' not in s:p.write_text(lines[0]+'\n'+intro+''.join(lines[2:]))
for name in progress['reports']:
 p=ROOT/name;s=p.read_text()
 if 'Final validation:' not in s:s+='\nFinal validation: **360 full marginal tests pass**. The current progress record accounts for all nine new independent standard-library replays.\n'
 p.write_text(s)
paths={'experiments/marginal_number_quotient.py','experiments/marginal_monotone_transfer.py','experiments/marginal_spin_reduction.py','tests/test_marginal_number_quotient.py','tests/test_marginal_monotone_transfer.py','research/marginal_bounded_quotient.md','research/marginal_monotone_transfer.md','research/marginal_blocker_progress.md','research/marginal_joint_metric.md','results/marginal_final_validation.json','results/marginal_h6/polynomial_metric/joint_progress.json'}
parent=BASE.parent;paths.update(str(p.relative_to(ROOT)) for p in parent.rglob('*') if p.is_file());mpath=ROOT/'results/marginal_progress_manifest.json';manifest=read(mpath)
for path in sorted(paths):manifest['sha256'][path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
manifest.update(file_count=len(manifest['sha256']),last_refreshed_paths=sorted(paths),last_refresh_scope='Bounded-degree number quotient and exact witnesses, coefficient-space transfer discovery, monotone hopping-family gap and energy integration, validated source/tests/reports, and labeled numerical proposals.');save(mpath,manifest)
for path in paths:assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==manifest['sha256'][path]
print(json.dumps({'tests':360,'seconds':seconds,'new_stdlib_replays':9,'energy_ledger':145,'manifest_files':manifest['file_count'],'refreshed_paths':len(paths),'goal_complete':False},indent=2))
