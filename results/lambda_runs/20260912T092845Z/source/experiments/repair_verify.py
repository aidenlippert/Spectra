"""Replay the representation repair and verify historical artifacts unchanged."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json,re
from .repair_realization import unpack_model,verify_against_records,learn,pack,mv
from .repair_run import equation_check,ports,alt_ports,predict,feasible
from .mission_policy_cover import Sample,certify

ROOT=Path(__file__).resolve().parents[1]


def read(path):return json.loads((ROOT/path).read_text())
def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def require(condition,message):
    if not condition:raise AssertionError(message)


def run():
    result=read('results/repair/result.json')
    initial=read('results/repair/initial_evidence.json')
    data=read('results/repair/acquisition_evidence.json')
    first=read('results/repair/initial_response.json')
    second=read('results/repair/learned_response.json')
    model=unpack_model(second['model'])
    evidence=verify_against_records(model,data)
    require(pack(model)==result['model'],'saved model differs from checked model')
    require(pack(evidence)==result['evidence_certificate'],'evidence receipt differs')
    require(pack(equation_check(model))==result['equation_certificate'],'equation receipt differs')
    require(first==result['first_response']==pack(learn(initial)),'initial refusal differs')
    require(first['status']=='needs_evidence','initial evidence unexpectedly accepted')

    old=read('results/mission/mathematical_examples.json')
    old_thermal=next(x for x in old['rows'] if x['name']=='two_node_thermal')
    require(len(initial['records'])==len(old_thermal['result']['samples'])==9,'initial record count')
    for sample,record in zip(old_thermal['result']['samples'],initial['records']):
        u=[F(x) for x in sample['point']]
        intervals=[[F(x) for x in pair] for pair in sample['intervals']]
        require(all(lo==hi for lo,hi in intervals),'old data not exact')
        y=[intervals[2][0]+F(3,10),F(1,10)-intervals[1][0]]
        require(u==[F(x) for x in record['controls']],'initial controls changed')
        require([[x,x] for x in y]==[[F(x) for x in pair] for pair in record['output_intervals']],'initial outputs changed')
        require(ports(u)==alt_ports(u)==y,'compatible-model witness does not match old evidence')
    require(data['records'][:9]==initial['records'],'initial data omitted from acquisition')
    experiment=[F(x) for x in first['experiment']['controls']]
    require(len(data['records'][9:])==len(experiment)==8,'experiment record count')
    for length,record in enumerate(data['records'][9:],1):
        u=[F(x) for x in record['controls']]
        require(u==experiment[:length],'acquisition departed from requested experiment')
        require([[v,v] for v in ports(u)]==[[F(x) for x in pair] for pair in record['output_intervals']],
                'saved mathematical experiment changed')

    policy=result['policy'];prep=[F(x) for x in policy['preparation']]
    operating=[F(x) for x in policy['operating']];full=prep+operating
    other=[F(x) for x in policy['initial_hot_only_alias']['other_preparation']]
    y=predict(model,full);y_alt=alt_ports(full);y_other=ports(other+operating)
    require(all(full!=[F(x) for x in record['controls']] for record in data['records']),
            'policy answer was supplied as an acquisition record')
    require(y==ports(full)==[F(x) for x in policy['final_ports']],'policy endpoint mismatch')
    require(y_alt==[F(x) for x in policy['old_evidence_compatible_alternative_final']],'alternate endpoint mismatch')
    require(ports(prep)[0]==ports(other)[0],'histories do not alias under old representation')
    require(feasible(y) and not feasible(y_alt) and not feasible(y_other),'policy status change failed')
    discrepancy=abs(ports(prep+[F(0)])[0]-ports(other+[F(0)])[0])
    require(discrepancy==F(policy['initial_hot_only_alias']['one_step_output_difference']),
            'future-output separation mismatch')
    now=mv(model['output'],model['input'])
    earlier=mv(model['output'],mv(model['transition'],model['input']))
    lipschitz=[F(2),abs(now[1])+abs(earlier[1]),abs(now[0])+abs(earlier[0])]
    values=[sum(operating),F(1,10)-y[1],y[0]-F(3,10)]
    sample=Sample(tuple(operating),tuple((v,v) for v in values))
    bridge=certify([''],[sample],lipschitz,2,F(0))
    require(pack(bridge)==result['bridge']['receipt'],'bridge replay differs')
    require(bridge['status']=='unresolved' and all(hi<=0 for lo,hi in bridge['policy_bounds'][1:]),
            'point feasibility/global comparison distinction failed')

    for p,h in result['source_sha256'].items():require(digest(p)==h,'repair source drift: '+p)
    for p,h in result['input_sha256'].items():require(digest('results/repair/'+p)==h,'repair input drift: '+p)
    require(digest('results/mission/mathematical_examples.json')==result['old_evidence_sha256'],'old evidence drift')
    source_checks={};archive_checks={}
    for version in ('v8','v9','v10','v11','v12','mission'):
        receipt=read(f'results/{version}/verification_receipt.json')
        for p,h in receipt['source_sha256'].items():source_checks[p]=digest(p)==h
    for version in ('v7','v8','v9','v10','v11','v12'):
        hashes=read(f'results/{version}/verification_receipt.json')['archive_sha256']
        if isinstance(hashes,str):hashes={'certificate_archive.json':hashes}
        for p,h in hashes.items():archive_checks[f'results/{version}/{p}']=digest(f'results/{version}/{p}')==h
    v7=read('results/v7/verification_receipt.json');v6=read('results/v6/verification_receipt.json')
    v5=read('results/v5/verification_receipt.json')
    old_hashes={**v7['source_hashes'],**v6['source_hashes'],**v6['data_hashes'],
        'results/v6/headroom_results.json':v6['headroom_result_sha256'],**v5['source_hashes'],**v5['result_hashes']}
    historical_checks={p:digest(p)==h for p,h in old_hashes.items()}
    require(all(source_checks.values()) and all(archive_checks.values()) and all(historical_checks.values()),
            'historical artifact drift')
    log=(ROOT/'results/repair/tests.log').read_text();counts=re.findall(r'Ran (\d+) tests',log)
    require(bool(counts) and log.rstrip().endswith('OK'),'test suite did not pass')
    sources=list(result['source_sha256'])+['experiments/repair_verify.py','tests/test_repair_realization.py']
    artifact_names=['initial_evidence.json','acquisition_evidence.json','initial_response.json',
                    'learned_response.json','result.json','tests.log']
    receipt=dict(status='verified',evidence_replay=evidence,policy_endpoint=pack(y),
        compatible_alternative_endpoint=pack(y_alt),hot_only_minimax_error_lower_bound=str(discrepancy/2),
        bridge_replayed=True,feasible_policy=True,global_comparison='unresolved',
        historical_hash_checks=historical_checks,historical_source_checks=source_checks,
        historical_archive_checks=archive_checks,tests_passed=int(counts[-1]),
        source_sha256={p:digest(p) for p in sources},
        artifact_sha256={p:digest('results/repair/'+p) for p in artifact_names},
        physical_validation=False,novel_algorithm=False,compounding_demonstrated=False)
    (ROOT/'results/repair/verification_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(status='verified',feasible_policy=True,global_comparison='unresolved',
                         tests_passed=int(counts[-1]),historical_hashes=len(historical_checks),
                         historical_sources=len(source_checks),historical_archives=len(archive_checks)),indent=2))


if __name__=='__main__':run()
