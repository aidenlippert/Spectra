"""Standalone exact accepting replay, copied to the root of the release bundle."""
import argparse,copy,hashlib,json,sys,time,unittest
from fractions import Fraction as F
from pathlib import Path


def run(output):
    start=time.monotonic();base=Path(__file__).resolve().parent;output=output.resolve()
    if output.exists() or output==base or base in output.parents:
        raise ValueError('Output must be new and outside the bundle')
    manifest=json.loads((base/'MANIFEST.json').read_text())
    for name,expected in manifest['files_sha256'].items():
        path=(base/name).resolve()
        if base not in path.parents or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('Bundle integrity failure: '+name)
    # Import only after file integrity is checked. This is corruption detection,
    # not authentication against an adversary who can replace the manifest too.
    from research.direct_control_20260916.collective_control import replay,construct
    from research.direct_control_20260916.collective_obstruction import calculate
    from research.direct_control_20260916.transfer_collective import choose_pair
    from research.direct_control_20260916.tensor_operator import build,expand_products,product_terms
    output.mkdir(parents=True);rows=[]
    for case in ('h8_original','h6_asymmetric','water_asymmetric'):
        directory=base/'cases'/case
        data=json.loads((directory/'fixture.json').read_text());state=json.loads((directory/'state.json').read_text())
        artifact=json.loads((directory/'certificate.json').read_text());st=time.monotonic()
        if case=='h8_original':pair=(6,10);cap=F(64)
        else:pair,_,_=choose_pair(data,state);cap=F(1024)
        fresh=replay(data,state,artifact,cap,expected_pair=pair)
        if not fresh['accepted']:raise AssertionError('Expected strong-control case refused')
        (output/(case+'.json')).write_text(json.dumps(fresh,indent=2)+'\n')
        row={'case':case,'pair_even_spin_indices':list(pair),'seconds':time.monotonic()-st,
             'protocol':fresh['attempts'][-1],'enumerated_configurations':fresh['stats']['enumerated_determinants'],
             'many_body_matrix_entries':fresh['stats']['many_body_matrix_entries']}
        rows.append(row)
        if case=='h8_original':
            h8=(data,state,artifact)
    data,state,artifact=h8
    refused=construct(data,state,max_amplitude=F(1,2))
    if refused['accepted']:raise AssertionError('Low-amplitude refusal failed')
    obstruction=calculate(data,state,artifact)
    if not obstruction['all_durations_cannot_certify']:raise AssertionError('Envelope obstruction failed')
    (output/'obstruction.json').write_text(json.dumps(obstruction,indent=2)+'\n')
    rejected=[]
    for label,change in (
      ('endpoint',lambda a:a['attempts'][-1].update(D_interval=['-2','-3/5'])),
      ('amplitude',lambda a:a['attempts'][-1].update(v_Ha='1/2')),
      ('uncertainty',lambda a:a['uncertainties'].update(initial_trace_distance='0'))):
        bad=copy.deepcopy(artifact);change(bad)
        try:replay(data,state,bad,F(64))
        except ValueError:rejected.append(label)
        else:raise AssertionError('Corrupt certificate accepted')
    bad=copy.deepcopy(artifact);bad['initial_Y']='1'
    try:calculate(data,state,bad)
    except ValueError:rejected.append('obstruction initial Y')
    else:raise AssertionError('Invalid obstruction premise accepted')
    mpo=build(data)
    if expand_products(mpo)!=product_terms(data):raise AssertionError('MPO reconstruction failed')
    tests=unittest.defaultTestLoader.loadTestsFromNames([
       'research.direct_control_20260916.test_tensor_operator',
       'research.direct_control_20260916.test_collective_control'])
    with (output/'tests.log').open('w') as stream:
        result=unittest.TextTestRunner(stream=stream,verbosity=2).run(tests)
    if not result.wasSuccessful():raise AssertionError('Focused test failure')
    if any(x in sys.modules for x in ('numpy','scipy','quimb','pyscf')):
        raise AssertionError('Numerical library on accepting path')
    receipt={'status':'accepted_strong_control_component_only','cases':rows,
       'all_durations_weak_control_envelope_obstruction':True,'corruptions_rejected':rejected,
       'tests_run':result.testsRun,'mpo_maximum_bond':max(mpo['widths']),
       'seconds':time.monotonic()-start,'arithmetic':'Python integers and Fractions',
       'original_weak_control_dynamics_solved':False,'general_many_body_problem_solved':False,
       'dependencies':'Supplied rational Hamiltonians and initial MPS states. Their discovery is inherited. Stronger control amplitudes are part of these new task specifications.',
       'scope':'No state trajectories or sector matrices in this accepting replay. The single-determinant norm witness is used only for the proof-family obstruction.'}
    (output/'complete.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='cases'}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);args=p.parse_args();run(args.out)
