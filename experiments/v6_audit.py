"""Post-run audits; no model tuning or autonomous method acquisition."""
import ast
from pathlib import Path
from hashlib import sha256
import json
import re
import numpy as np
from experiments.v6_headroom import run,measured_traces,physical_system,HORIZON

ROOT=Path(__file__).resolve().parents[1]


def excitation(train,planned):
    """Check whether planned input contrasts are identified by training inputs.

    This necessary condition is not a sufficient identifiability certificate.
    It uses input values only and never evaluates target outputs.
    """
    U=np.concatenate([t.u for t in train]);mean=U.mean(axis=0)
    _,values,vh=np.linalg.svd(U-mean,full_matrices=False)
    rank=int(np.sum(values>max(float(values[0]),1.)*1e-10))
    basis=vh[:rank]
    flags=[]
    for trace in planned:
        for start in range(8,len(trace.u)-HORIZON+1,HORIZON):
            contrast=trace.u[start:start+HORIZON]-mean
            remainder=contrast-(contrast@basis.T)@basis
            flags.append(bool(np.max(np.abs(remainder))>1e-8))
    return {'training_input_contrast_rank':rank,'input_dimension':U.shape[1],
            'planned_blocks':len(flags),'blocks_with_unidentified_input_contrast':sum(flags),
            'uses_future_output':False}


def local_dependencies(roots):
    pending=list(roots);seen=set();edges={}
    while pending:
        name=pending.pop()
        if name in seen:continue
        seen.add(name);path=ROOT/(name.replace('.','/')+'.py')
        if not path.exists():raise AssertionError('missing local module')
        tree=ast.parse(path.read_text());deps=set()
        for node in ast.walk(tree):
            if isinstance(node,ast.ImportFrom) and node.module and node.module.startswith('experiments.'):
                deps.add(node.module)
            if isinstance(node,ast.Import):
                deps.update(a.name for a in node.names if a.name.startswith('experiments.'))
        edges[name]=sorted(deps);pending.extend(deps-seen)
        if len(seen)>30:raise AssertionError('unexpected dependency expansion')
    return edges


def strip_timing(value):
    if isinstance(value,dict):
        return {k:strip_timing(v) for k,v in value.items() if not k.endswith('_seconds')}
    if isinstance(value,list):return [strip_timing(v) for v in value]
    return value


def audit():
    saved=json.loads((ROOT/'results/v6/headroom_results.json').read_text())
    replay=run()
    assert strip_timing(saved)==strip_timing(replay),'saved result did not reproduce'
    assert saved['headroom_gate']=='not_established'
    assert not saved['autonomous_acquisition_performed']
    previous=json.loads((ROOT/'results/v5/verification_receipt.json').read_text())
    preserved=0
    for group in ('source_hashes','result_hashes'):
        for path,digest in previous[group].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==digest,'V5 regression fixture changed'
            preserved+=1
    dependencies=local_dependencies(('experiments.v5_run','experiments.v5_complementarity','experiments.v5_account'))
    assert all(not module.startswith('experiments.v4') for module in dependencies)
    train,selection,cal,test=measured_traces()
    physical=excitation(train,test)
    synthetic=[]
    for seed in range(62001,62033):
        train,_,_,test,_=physical_system(seed)
        synthetic.append(excitation(train,test))
    assert physical['training_input_contrast_rank']==0
    assert all(s['training_input_contrast_rank']==2 for s in synthetic)
    log=(ROOT/'results/v6/test_log.txt').read_text()
    assert log.rstrip().endswith('OK')
    tests=int(re.search(r'Ran (\d+) tests',log).group(1))
    sources=sorted(list(ROOT.glob('experiments/v6*.py'))+list(ROOT.glob('tests/test_v6*.py')))
    data=sorted((ROOT/'data/v6').glob('*'))
    result={'stage':'post-run verification and validity diagnosis','prospective_gate_revised':False,
            'tests_passed':tests,'saved_full_protocol_replayed':True,
            'prior_v5_source_and_result_hashes_preserved':preserved,
            'v5_transitive_local_imports':dependencies,'v4_calibration_required_by_v5':False,
            'v4_576000_episode_cost_category':'separate historical experiment, not a V5 or V6 runtime dependency',
            'measured_input_excitation':physical,'synthetic_input_excitation':synthetic,
            'physical_diagnosis':'constant training heater inputs cannot identify response to new input contrasts without additional priors',
            'autonomous_method_acquisition':'not_run_headroom_gate_failed',
            'm1_causes_cheaper_m2_acquisition':'not_tested','two_frozen_method_transfer_gains':'not_tested',
            'method_acquisition_net_benefit':'not_tested','two_generation_pass_succeeded':False,
            'physical_measurements_new':0,
            'source_hashes':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources},
            'data_hashes':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in data},
            'headroom_result_sha256':sha256((ROOT/'results/v6/headroom_results.json').read_bytes()).hexdigest()}
    (ROOT/'results/v6/verification_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    print(json.dumps({k:v for k,v in audit().items() if k not in ('synthetic_input_excitation','source_hashes','data_hashes')},indent=2))
