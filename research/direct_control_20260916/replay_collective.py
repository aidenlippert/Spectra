"""Fresh exact replay, including refusal and certificate-mutation gates."""
import json,time,argparse,copy,sys
from fractions import Fraction as F
from pathlib import Path
from .collective_control import replay,construct
from .collective_obstruction import calculate


def run(root,out):
    start=time.monotonic();inputs=root/'imported/Spectra_control_reduction/inputs'
    data=json.loads((inputs/'fixture.json').read_text());state=json.loads((inputs/'state.json').read_text())
    artifact=json.loads((root/'collective_control_bound64.json').read_text())
    fresh=replay(data,state,artifact,F(64));accepted_seconds=time.monotonic()-start
    refusal=construct(data,state,max_amplitude=F(1,2))
    obstruction=calculate(data,state,artifact)
    if refusal['accepted'] or not obstruction['family_cannot_certify'] or not obstruction['all_durations_cannot_certify']:raise AssertionError('Required refusal did not reproduce')
    mutations=[]
    for name,mutate in [
      ('changed bound',lambda a:a['attempts'][-1].update(D_interval=['-2','-3/5'])),
      ('changed control',lambda a:a['attempts'][-1].update(v_Ha='1/2')),
      ('changed initialization allowance',lambda a:a['uncertainties'].update(initial_trace_distance='0'))]:
        bad=copy.deepcopy(artifact);mutate(bad)
        try:replay(data,state,bad,F(64))
        except ValueError:mutations.append(name)
        else:raise AssertionError('Corrupt certificate accepted: '+name)
    bad=copy.deepcopy(artifact);bad['initial_Y']='1'
    try:calculate(data,state,bad)
    except ValueError:mutations.append('obstruction nonzero initial Y')
    else:raise AssertionError('Obstruction accepted an incompatible initial Y')
    if any(x in sys.modules for x in ('numpy','scipy','quimb','pyscf')):raise AssertionError('Numerical library on exact acceptance path')
    result={'status':'accepted_strong_control_and_expected_low_amplitude_refusal',
       'accepted_replay_seconds':accepted_seconds,'total_with_mutation_tests_seconds':time.monotonic()-start,
       'selected_protocol':fresh['attempts'][-1],'accepted':fresh['accepted'],'low_amplitude_target_proved':refusal['accepted'],
       'commutator_sum_family_obstruction':obstruction['family_cannot_certify'],'mutations_rejected':mutations,
       'all_durations_commutator_sum_obstruction':obstruction['all_durations_cannot_certify'],
       'exact_arithmetic':True,'enumerated_configurations':fresh['stats']['enumerated_determinants'],
       'many_body_matrix_entries':fresh['stats']['many_body_matrix_entries'],
       'frozen_four_phase_problem_solved':False,
       'scope':'Strong-control variant. The original amplitude-.5 four-phase problem remains unresolved without enumeration.'}
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='selected_protocol'}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('results/direct_control_20260916'));p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.root,a.output)
