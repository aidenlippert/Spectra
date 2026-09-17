"""Audit input bytes, classify artifacts and summarize the measured run ledger."""
import hashlib,json,platform,zipfile
from datetime import datetime,timezone
from pathlib import Path


def save(path,value):path.write_text(json.dumps(value,indent=2)+'\n')


def run():
    repo=Path(__file__).resolve().parents[2];root=repo/'results/direct_control_20260916'
    imported=json.loads((root/'import.json').read_text());archive=Path(imported['source'])
    digest=hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest!=imported['archive_sha256']:raise AssertionError('Source archive changed')
    checked=0
    with zipfile.ZipFile(archive) as z:
        for info in z.infolist():
            if info.is_dir():continue
            relative=Path(info.filename)
            if relative.is_absolute() or '..' in relative.parts:raise ValueError('Archive path')
            path=root/'imported'/relative
            if path.read_bytes()!=z.read(info):raise AssertionError('Imported bytes changed: '+info.filename)
            checked+=1
    save(root/'source_integrity.json',{'source_archive_sha256':digest,'imported_files_compared_byte_for_byte':checked,'all_unchanged':True})
    exclusions={
       'tensor_operator/summary.json':'Invalid early CAR encoding; false eight-word result. Replaced by exact_h8_mpo.json and independent action tests.',
       'canonical_schmidt.json':'Invalid square-root scale argument. Replaced by canonical_schmidt_fixed.json.',
       'collective_control.json':'Superseded prototype metadata. Use collective_control_bound64.json and fresh exact replay.',
       'collective_obstruction.json':'Earlier fixed-angle-only calculation. Use the final checked all-duration obstruction.',
       'collective_obstruction_all_durations.json':'Earlier prototype without standalone initial-moment rebinding. Use collective_obstruction_final.json or standalone_replay/obstruction.json.',
       'local_basis/':'Failed first local-basis attempt. Valid output is local_basis_retry/.',
       'Spectra_direct_control_component.zip':'Superseded preliminary package; not the final independently extracted bundle.',
       'early subagent MPS and commutator claims':'Rejected CAR/tensor-scale/denominator mistakes; removed or replaced implementations. No accepted result depends on them.'}
    save(root/'artifact_disposition.json',{'excluded_or_superseded':exclusions,
      'exact_positive_results':['collective_control_bound64.json','collective_transfer/h6_asymmetric.json','collective_transfer/water_asymmetric.json','standalone_replay/complete.json'],
      'exact_operator_component':'exact_h8_mpo.json',
      'exact_scoped_obstruction':'standalone_replay/obstruction.json',
      'conditional_numeric_diagnostics':['canonical_schmidt_fixed.json','local_schmidt_fixed.json','local_basis_retry/record.json'],
      'numeric_reference_only':'collective_oracle.json',
      'original_reference_only':'reference_replay/complete.json'})
    runs=[json.loads(p.read_text()) for p in sorted((root/'runs').glob('*.json'))]
    accounting={'as_of_UTC':datetime.now(timezone.utc).isoformat(),'host':platform.platform(),
       'recorded_processes':len(runs),'sum_recorded_wall_seconds':sum(r.get('wall_seconds',0) for r in runs),
       'sum_recorded_user_seconds':sum(r.get('child_user_seconds',0) for r in runs),
       'sum_recorded_system_seconds':sum(r.get('child_system_seconds',0) for r in runs),
       'maximum_recorded_process_RSS_bytes':max(r.get('peak_child_RSS_bytes',0) for r in runs),
       'runs':runs,'paid_compute_used':False,'new_libraries_installed':False,
       'excluded_costs':['Inherited integrals and rational-Hamiltonian preparation','Inherited initial-MPS discovery','Prior one-particle orbital proposal','Reasoning, editing, inspection and packaging','Early unmetered subagent diagnostics and tests'],
       'interpretation':'Sum of measured sequential child-process clocks, including failures and repeats. Not complete cold-discovery cost and not a matched-task speedup.'}
    save(root/'ACCOUNTING.json',accounting)
    standalone=json.loads((root/'standalone_replay/complete.json').read_text())
    save(root/'RESULT.json',{'original_weak_control_enumeration_free_goal_solved':False,
       'general_many_body_problem_solved':False,'strong_control_component_certified':True,
       'strong_control_transfer_cases':[r['case'] for r in standalone['cases']],
       'H8_control_amplitude_Ha':'64','original_H8_amplitude_cap_Ha':'1/2',
       'all_duration_constant_W_norm_envelope_obstruction':True,
       'standalone_replay_seconds':standalone['seconds'],'standalone_standard_library_tests':standalone['tests_run'],
       'full_focused_suite_tests':17,'source_archive_unchanged':True,
       'status':'Conditional component success; main original objective still open'})
    manifest={}
    paths=list((repo/'research/direct_control_20260916').glob('*.py'))+list((repo/'research/direct_control_20260916').glob('*.md'))
    for pattern in ('collective_control_bound64.json','collective_transfer/*.json','collective_replay_final.json','collective_obstruction_final.json','exact_h8_mpo.json','canonical_schmidt_fixed.json','local_schmidt_fixed.json','local_basis_retry/record.json','source_integrity.json','artifact_disposition.json','ACCOUNTING.json','RESULT.json','standalone_replay/*.json','Spectra_direct_control_verified.zip'):
        paths.extend(root.glob(pattern))
    for path in sorted(set(paths)):
        manifest[str(path.relative_to(repo))]=hashlib.sha256(path.read_bytes()).hexdigest()
    save(root/'MANIFEST.json',{'kind':'reviewed_work_artifact_hashes_v1','files_sha256':manifest,
       'note':'Integrity index, not a substitute for fresh proof replay. Numeric diagnostics retain their stated arithmetic assumptions.'})
    print(json.dumps({'source_files_unchanged':checked,'recorded_processes':len(runs),'sum_recorded_wall_seconds':accounting['sum_recorded_wall_seconds'],'peak_RSS_bytes':accounting['maximum_recorded_process_RSS_bytes'],'standalone_seconds':standalone['seconds']}))


if __name__=='__main__':run()
