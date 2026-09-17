"""Receipt/dependency audit and complete instrumented campaign accounting.

This does not replace mathematical replay. It checks bindings, exact preservation
of inherited columns, cost records, and the claimed finite control objective.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import sys
from research.intervention_reduction_20260916.exact import digest

ROOT=Path('results/intervention_reduction_20260916')
SOURCES=Path('research/intervention_reduction_20260916')
LEGACY=Path('results/transfer_solver_20260915')


def read(path):
    return json.loads(path.read_text())


def audit():
    inputs={};models={}
    for case in ('h6_asymmetric','h8_cold','water_asymmetric'):
        base=LEGACY/'cases'/case;data=read(base/'fixture.json');state=read(base/'mps/state.json')
        if state['fixture_sha256']!=digest(data):
            raise AssertionError('Inherited state fixture binding')
        generate=read(LEGACY/'runs'/f'{case}_generate.json');discovery=read(LEGACY/'runs'/f'{case}_state.json')
        inputs[case]={'fixture_sha256':digest(data),'inherited_MPS_sha256':digest(state),
            'integral_generation_historical_process_seconds':generate['wall_seconds'],
            'state_stage_historical_process_seconds':discovery['wall_seconds'],
            'interpretation':'archived input-acquisition costs; not a fresh complete timing for this campaign'}
        models[digest(data)]=case
    proposals={}
    for path in ROOT.rglob('proposal_r*.json'):
        if '.construction.' in path.name:
            continue
        p=read(path)
        if p.get('kind')!='integer_control_subspace_v1':
            continue
        if p['fixture_sha256'] not in models:
            raise AssertionError('Unknown source model')
        proposals[digest(p)]=(path,p)
    for path in ROOT.rglob('discovery.json'):
        d=read(path)
        if 'inherited_MPS_sha256' in d:
            expected=inputs[models[d['fixture_sha256']]]['inherited_MPS_sha256']
            if expected!=d['inherited_MPS_sha256']:
                raise AssertionError('Inherited MPS changed after discovery')
    preservation=[]
    for path in ROOT.rglob('*.construction.json'):
        d=read(path)
        if 'parent_sha256' not in d:
            continue
        parent_path,parent=proposals[d['parent_sha256']];child_path,child=proposals[d['proposal_sha256']]
        rank=len(parent['vectors'][0]);old=dict(zip(parent['configurations'],parent['vectors']))
        new=dict(zip(child['configurations'],child['vectors']))
        okay=(parent['denominator']==child['denominator'] and parent['controls']==child['controls']
              and all(s in new for s in old)
              and all(row[:rank]==old.get(s,[0]*rank) for s,row in new.items()))
        if not okay:
            raise AssertionError('Enrichment changed an inherited column or control')
        preservation.append({'parent':str(parent_path),'child':str(child_path),'all_old_columns_preserved_exactly':True})
    receipts=[];oracles=[]
    for path in ROOT.rglob('*.json'):
        if path.parent.name=='runs':
            continue
        d=read(path)
        if not isinstance(d,dict):
            continue
        if d.get('status') in ('accepted_molecular_trajectory','accepted_control_neighborhood','accepted_uniform_control_reduction'):
            if d['fixture_sha256'] not in models or d['proposal_sha256'] not in proposals:
                raise AssertionError('Unbound accepting receipt')
            receipts.append({'path':str(path),'status':d['status'],'target_met':d['target_met'],
                'rank':d['reduced_dimension'],'state_bound':d.get('normalized_state_error_float',d.get('state_vector_error_float')),
                'uniform_in_time_target_met':d.get('uniform_in_time_target_met')})
        if 'inside_state_bound' in d:
            if not d['inside_state_bound'] or not d['inside_population_interval']:
                raise AssertionError('Failed independent numerical diagnostic')
            oracles.append({'path':str(path),'reference_dimension':d['reference_dimension'],
                'observed_state_error':d['observed_state_error'],'certified_state_error':d['certified_state_error']})
    comparison=read(ROOT/'h6_matched_comparison.json')
    if not comparison['all_mathematical_fields_identical']:
        raise AssertionError('Compiled replay mismatch')
    for directory in ('h6_bundle','water40_bundle'):
        for path in sorted((ROOT/directory).glob('query_*.json')):
            if not read(path)['uniform_in_time_target_met']:
                raise AssertionError('Claimed final response target not met')
    design=read(ROOT/'h6_design/design.json');robust=read(ROOT/'h6_functional_robustness/query_1.json')
    requested=F(design['requested_population_increase']);lower=F(robust['population_change_interval'][0])
    if lower<requested:
        raise AssertionError('Functional inverse-control target not met')
    policy={'kind':'certified_finite_model_policy_v1','fixture_sha256':robust['fixture_sha256'],
        'proposal_sha256':robust['proposal_sha256'],'initial_state':robust['initial_state'],
        'observed_spatial_orbital':robust['observed_spatial_orbital'],
        'nominal_schedule':design['top_candidates'][0]['schedule'],
        'control_operators_source':'results/intervention_reduction_20260916/h6_snapshot/proposal_r24.json',
        'requested_population_increase':str(requested),'certified_population_increase_interval':robust['population_change_interval'],
        'certified_population_increase_interval_float':robust['population_change_interval_float'],
        'functional_target_met':True,'strict_0p005_state_target_met':robust['target_met'],
        'allowed_uniform_control_deviation_Ha':'1/2000',
        'amplitude_box_per_control_Ha':['-1/100','1/100'],
        'perturbation_scope':'all measurable deviations of at most 1/2000 Ha per channel, while remaining in the original amplitude box',
        'proof_receipt':str(ROOT/'h6_functional_robustness/query_1.json'),
        'global_control_optimality_proved':False,'laboratory_control_mapping_provided':False,
        'physical_model_error_certified':False}
    (ROOT/'CONTROL_POLICY.json').write_text(json.dumps(policy,indent=2)+'\n')
    runs=[read(p) for p in sorted((ROOT/'runs').glob('*.json'))]
    if any(r['status'] not in ('passed','failed','timeout') for r in runs):
        raise AssertionError('Unfinished instrumented process')
    report={'status':'receipt_and_dependency_audit_passed','input_dependencies':inputs,
        'enrichment_preservation':preservation,'accepting_receipts':receipts,'independent_numerical_diagnostics':oracles,
        'instrumented_process_count':len(runs),'process_wall_seconds':sum(r['wall_seconds'] for r in runs),
        'process_CPU_seconds':sum(r.get('child_user_seconds',0)+r.get('child_system_seconds',0) for r in runs),
        'peak_single_child_RSS_bytes':max(r.get('peak_child_RSS_bytes',0) for r in runs),
        'process_failures':[{'name':r['name'],'status':r['status'],'wall_seconds':r['wall_seconds']} for r in runs if r['status']!='passed'],
        'valid_certificates_missing_numerical_target':sum(not r['target_met'] for r in receipts),
        'unexecuted_candidate_proposals_counted_as_success':False,
        'accounting_scope':'all instrumented current-campaign processes, including failed targets, repeated validations and references; not time for one fresh solution',
        'excluded_costs':'human/agent reasoning, editing, browsing, filesystem inspection, documentation, this audit, unmetered focused tests, and earlier project R&D; inherited input acquisition is listed separately',
        'host':{'platform':platform.platform(),'machine':platform.machine(),'Python':sys.version},
        'runs':runs}
    (ROOT/'ACCOUNTING.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','instrumented_process_count','process_wall_seconds','process_CPU_seconds','peak_single_child_RSS_bytes','valid_certificates_missing_numerical_target')}))


if __name__=='__main__':
    audit()
