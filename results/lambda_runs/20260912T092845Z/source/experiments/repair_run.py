"""Observable-data acquisition, independent equation check and policy repair."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
from hashlib import sha256
from time import perf_counter
import subprocess,sys,json
from .mission_examples import thermal
from .mission_policy_cover import Sample,certify
from .repair_realization import pack,unpack_model,mv,mm,rank,state_after,step,verify,verify_against_records
ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'results/repair'


def ports(controls):
    # Use existing task-output interface; discard returned internal history.
    values,_=thermal(controls)
    return [values[2]+F(3,10),F(1,10)-values[1]]


def row(controls,outputs):
    return dict(controls=pack(controls),output_intervals=[[str(v),str(v)] for v in outputs],reset_confirmed=True)


def predict(model,controls):return mv(model['output'],state_after(model,list(controls),reset_confirmed=True))
def feasible(y):return y[1]>=F(1,10) and y[0]<=F(3,10)


def equation_check(model):
    # These matrices belong ONLY to the independent model-based checker.
    A=[[F(1,2),F(1,4)],[F(1,4),F(1,2)]];B=[F(1,2),F(0)]
    identity=[[F(1),F(0)],[F(0),F(1)]];lift=[]
    for port,lag in model['selected_tests']:
        v=identity[port][:]
        for _ in range(lag):v=mm([v],A)[0]
        lift.append(v)
    if mm(lift,A)!=mm(model['transition'],lift):raise ValueError('transition intertwining failed')
    if mv(lift,B)!=model['input']:raise ValueError('input intertwining failed')
    if mm(model['output'],lift)!=identity:raise ValueError('port preservation failed')
    if mv(lift,[F(0),F(0)])!=model['reset']:raise ValueError('reset failed')
    # Linear specialization of constrained invariant-row-space construction.
    basis=[[F(1),F(0)]]
    while True:
        old=len(basis)
        for r in list(basis):
            v=mm([r],A)[0]
            if rank(basis+[v])>len(basis):basis.append(v)
        if len(basis)==old:break
    return dict(status='exact_model_intertwining',lift=lift,known_equation_linear_lumping_rank=len(basis),
        scope='all finite input histories and causal observed-port input feedback under the unchanged supplied dynamics')


def alt_ports(controls):
    # An explicit passive, symmetric compatible model; never sent to learner.
    z=[F(0),F(0)]
    for u in controls:z=[z[0]/2+z[1]/4+u/2,z[0]/4+z[1]/4]
    return z


def invoke(input_name,output_name):
    start=perf_counter()
    subprocess.run([sys.executable,'-m','experiments.repair_realization',str(DEST/input_name),str(DEST/output_name)],
                   cwd=ROOT,check=True,capture_output=True,text=True,timeout=60)
    return json.loads((DEST/output_name).read_text()),perf_counter()-start


def run():
    DEST.mkdir(exist_ok=True)
    old_path=ROOT/'results/mission/mathematical_examples.json'
    old=json.loads(old_path.read_text())
    thermal_row=next(r for r in old['rows'] if r['name']=='two_node_thermal')
    records=[]
    for s in thermal_row['result']['samples']:
        us=[F(x) for x in s['point']];iv=[[F(x) for x in pair] for pair in s['intervals']]
        if any(lo!=hi for lo,hi in iv):raise AssertionError('old evidence is not exact')
        y=[iv[2][0]+F(3,10),F(1,10)-iv[1][0]]
        if ports(us)!=y or alt_ports(us)!=y:raise AssertionError('compatible model witness failed on old data')
        records.append(row(us,y))
    data=dict(prior='zero_reset_strictly_proper_LTI',maximum_dimension=4,output_ports=2,records=records)
    (DEST/'initial_evidence.json').write_text(json.dumps(data,indent=2)+'\n')
    first,t_first=invoke('initial_evidence.json','initial_response.json')
    if first['status']!='needs_evidence':raise AssertionError('initial evidence unexpectedly complete')
    experiment=first['experiment'];controls=[F(x) for x in experiment['controls']]
    observed=[];probe_start=perf_counter()
    for i in range(1,len(controls)+1):
        us=controls[:i];y=ports(us);observed.append(row(us,y))
    probe_seconds=perf_counter()-probe_start
    data['records']=records+observed
    (DEST/'acquisition_evidence.json').write_text(json.dumps(data,indent=2)+'\n')
    second,t_second=invoke('acquisition_evidence.json','learned_response.json')
    if second['status']!='constructed':raise AssertionError(second)
    checked=perf_counter()
    model=unpack_model(second['model']);moment_receipt=verify(model)
    evidence_receipt=verify_against_records(model,data)
    equation_receipt=equation_check(model)
    words=0
    for length in range(7):
        for word in product((F(0),F(1,2),F(1)),repeat=length):
            if predict(model,word)!=ports(word):raise AssertionError('fixed input-word check failed')
            words+=1
    # New bounded port-feedback inputs, generated from each arm's own outputs.
    learned_z=model['reset'];history=[]
    for _ in range(16):
        y=mv(model['output'],learned_z);actual_y=ports(history)
        u=max(F(0),min(F(1),F(1,2)-y[0]+y[1]/2))
        actual_u=max(F(0),min(F(1),F(1,2)-actual_y[0]+actual_y[1]/2))
        if u!=actual_u:raise AssertionError('feedback divergence')
        learned_z=step(model,learned_z,u);history.append(u)
        if mv(model['output'],learned_z)!=ports(history):raise AssertionError('feedback output mismatch')
    prep_a=[F(1),F(1,4)];prep_b=[F(0),F(3,4)];operating=None;attempts=[]
    candidates=sorted(product((F(0),F(1,4),F(1,2),F(3,4),F(1)),repeat=2),key=lambda u:(sum(u),u))
    for candidate in candidates:
        y=predict(model,prep_a+list(candidate));ok=feasible(y)
        attempts.append(dict(controls=candidate,predicted_ports=y,feasible=ok))
        if ok:operating=list(candidate);break
    if operating is None:raise AssertionError('no feasible operating policy')
    full=prep_a+operating;other=prep_b+operating
    if any([F(x) for x in r['controls']]==full for r in data['records']):raise AssertionError('policy evaluation leaked into acquisition')
    before_a=ports(prep_a);before_b=ports(prep_b)
    predicted=predict(model,full)
    actual=ports(full);alternative=alt_ports(full);other_y=ports(other)
    if predicted!=actual:raise AssertionError('independent policy evaluation disagrees with learned prediction')
    if before_a[0]!=before_b[0] or not feasible(actual) or feasible(alternative) or feasible(other_y):
        raise AssertionError('consequential ambiguity/repair failed')
    # Supply the previously assumed response oracle to the unchanged bridge.
    # Its global optimum remains unresolved; point feasibility is a separate fact.
    now=mv(model['output'],model['input'])
    earlier=mv(model['output'],mv(model['transition'],model['input']))
    L=[F(2),abs(now[1])+abs(earlier[1]),abs(now[0])+abs(earlier[0])]
    values=[sum(operating),F(1,10)-predicted[1],predicted[0]-F(3,10)]
    supplied=Sample(tuple(operating),tuple((v,v) for v in values))
    bridge=certify([''],[supplied],L,2,F(0))
    if bridge.get('policy')!=tuple(operating) or any(hi>0 for _,hi in bridge['policy_bounds'][1:]):
        raise AssertionError('bridge did not retain the feasible policy')
    proof_time=perf_counter()-checked
    files=['research/repair/PROTOCOL.md','experiments/repair_realization.py','experiments/repair_run.py',
           'experiments/mission_examples.py','experiments/mission_policy_cover.py','research/MISSION.md',
           'research/mission/CONSTRUCTIVE_BRIDGE.md']
    report=dict(status='constructed_and_checked',initial_records=len(records),first_response=first,
        experiment=dict(controls=controls,planned_streaming_resets=1,planned_streaming_steps=len(controls),scalar_observations=len(observed)*2,
                        exact_prefix_evaluator_calls=len(observed),prefix_step_evaluations=sum(range(1,len(controls)+1)),
                        evidence='supplied equation evaluations, not laboratory measurements'),
        model=model,moment_certificate=moment_receipt,equation_certificate=equation_receipt,
        evidence_certificate=evidence_receipt,
        fixed_words_checked=words,feedback_steps_checked=len(history),
        bridge=dict(receipt=bridge,learned_response_lipschitz=L,
                    response_source='learned model; independent equation evaluation used only to check it',
                    claim='point feasibility supported; complete domain optimization unresolved'),
        policy=dict(preparation=prep_a,operating=operating,total_control=sum(full),attempts=attempts,
            final_ports=actual,cold_margin=actual[1]-F(1,10),hot_margin=F(3,10)-actual[0],
            old_evidence_compatible_alternative_final=alternative,
            initial_hot_only_alias=dict(other_preparation=prep_b,current_port=before_a[0],
                other_final_ports=other_y,one_step_output_difference=ports(prep_a+[F(0)])[0]-ports(prep_b+[F(0)])[0]),
            before='no guarantee: compatible models disagree',after='feasible within the exact admitted model class',
            global_preparation_optimality='not_claimed'),
        cost_seconds=dict(initial_inference=t_first,additional_equation_calls=probe_seconds,
                          acquisition_subprocess=t_second,checking_and_evaluation=proof_time),
        input_sha256={name:sha256((DEST/name).read_bytes()).hexdigest() for name in ('initial_evidence.json','acquisition_evidence.json')},
        old_evidence_sha256=sha256(old_path.read_bytes()).hexdigest(),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        physical_validation=False,novel_algorithm=False,compounding_demonstrated=False,
        uncertainty=dict(representation='zero under checked equations',numerical='exact rational',
            state_estimation='zero only with verified reset and full input history',
            physical_model='not bounded by this experiment',measurements='exact mathematical data only',
            reaction_thermal_coupling='not identified or authorized by isolated records'))
    (DEST/'result.json').write_text(json.dumps(pack(report),indent=2)+'\n')
    print(json.dumps(pack({k:report[k] for k in ('status','initial_records','experiment','model','fixed_words_checked','policy','uncertainty')}),indent=2))


if __name__=='__main__':run()
