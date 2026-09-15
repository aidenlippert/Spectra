"""Privileged-method headroom diagnostic. No autonomous acquisition arm.

Hidden physical matrices are confined to generation. Model fitting receives
only input/output traces. Future measured outputs are not read inside rollout.
"""
from dataclasses import dataclass, field
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
ORDERS=(1,2,3,4,6,8)
RIDGES=(1e-8,1e-5,1e-2,1e-1)
HORIZON=16
METHODS=('privileged_history_state','stateful_arx','history_retrieval','sparse_history','current_state')


@dataclass
class Trace:
    u: np.ndarray
    y: np.ndarray

    def check(self):
        if self.u.ndim!=2 or self.y.ndim!=1 or len(self.y)!=len(self.u)+1:
            raise ValueError('input transitions and observed endpoints required')
        if len(self.u)<32 or not np.all(np.isfinite(self.u)) or not np.all(np.isfinite(self.y)):
            raise ValueError('finite sufficiently long traces required')


@dataclass
class Scaling:
    ym: float
    ys: float
    um: np.ndarray
    us: np.ndarray


def scaling(traces):
    for t in traces:t.check()
    y=np.concatenate([t.y for t in traces]);u=np.concatenate([t.u for t in traces])
    return Scaling(float(y.mean()),max(float(y.std()),1e-6),u.mean(axis=0),
                   np.where(u.std(axis=0)>1e-8,u.std(axis=0),1.))


def features(trace,p,norm):
    trace.check()
    if p not in ORDERS:raise ValueError('history order outside declared search')
    y=(trace.y-norm.ym)/norm.ys;u=(trace.u-norm.um)/norm.us
    times=np.arange(p-1,len(u))
    X=np.column_stack([y[times-i] for i in range(p)]+
                      [u[times-i,j] for i in range(p) for j in range(u.shape[1])]+[np.ones(len(times))])
    return X,y[times+1]


@dataclass
class Model:
    p: int
    norm: Scaling
    coef: np.ndarray | None = None
    bank_x: np.ndarray | None = None
    bank_y: np.ndarray | None = None
    neighbors: int = 1
    matrix: np.ndarray | None = None
    input_matrix: np.ndarray | None = None
    work: dict = field(default_factory=dict)


def compile_state(model):
    if model.coef is None:raise ValueError('linear law required')
    p=model.p;k=len(model.norm.um);d=p+(p-1)*k+1
    A=np.zeros((d,d));B=np.zeros((d,k))
    A[0,:p]=model.coef[:p]
    A[0,p:-1]=model.coef[p+k:-1]
    A[0,-1]=model.coef[-1];B[0]=model.coef[p:p+k]
    for i in range(1,p):A[i,i-1]=1
    if p>1:
        B[p:p+k]=np.eye(k)
        for i in range(1,p-1):A[p+i*k:p+(i+1)*k,p+(i-1)*k:p+i*k]=np.eye(k)
    A[-1,-1]=1
    model.matrix=A;model.input_matrix=B
    return model


def forecast(model,observed_y,previous_u,future_u,compiled=False):
    """Only the observed prefix and requested future inputs cross this boundary."""
    p=model.p;k=len(model.norm.um)
    if len(observed_y)<p or len(previous_u)<p-1 or future_u.ndim!=2 or future_u.shape[1]!=k:
        raise ValueError('insufficient state history or wrong input shape')
    yh=((np.asarray(observed_y)[-p:]-model.norm.ym)/model.norm.ys)[::-1].copy()
    uh=((np.asarray(previous_u)[-(p-1):]-model.norm.um)/model.norm.us)[::-1].copy() if p>1 else np.zeros((0,k))
    uf=(future_u-model.norm.um)/model.norm.us
    output=[]
    if compiled:
        if model.matrix is None:raise ValueError('uncompiled model')
        state=np.r_[yh,uh.ravel(),1.]
        for u in uf:
            state=model.matrix@state+model.input_matrix@u
            output.append(state[0])
    else:
        for u in uf:
            x=np.r_[yh,u,uh.ravel(),1.]
            if model.bank_x is None:pred=float(x@model.coef)
            else:
                distances=np.sum((model.bank_x-x)**2,axis=1)
                chosen=np.argpartition(distances,model.neighbors-1)[:model.neighbors]
                pred=float(model.bank_y[chosen].mean())
            output.append(pred)
            yh=np.r_[pred,yh[:-1]]
            if p>1:uh=np.vstack((u,uh[:-1]))
    return np.asarray(output)*model.norm.ys+model.norm.ym


def blocks(model,traces,compiled=False):
    errors=[]
    for trace in traces:
        for start in range(8,len(trace.u)-HORIZON+1,HORIZON):
            pred=forecast(model,trace.y[:start+1],trace.u[:start],trace.u[start:start+HORIZON],compiled)
            errors.append((pred-trace.y[start+1:start+HORIZON+1])/model.norm.ys)
    if not errors:raise ValueError('no complete forecast blocks')
    return np.asarray(errors)


def fit(traces,validation,method):
    if method not in METHODS:raise ValueError('unknown diagnostic method')
    began=perf_counter();norm=scaling(traces)
    work={'candidate_fits':0,'refits':0,'training_rows_processed':0,'normal_equation_work_proxy':0,
          'selection_forecast_blocks':0,'retrieval_distance_coordinates':0,'failed_candidates':0}
    candidates=[]
    orders=(1,) if method=='current_state' else ((8,) if method in ('history_retrieval','sparse_history') else ORDERS)
    for p in orders:
        pairs=[features(t,p,norm) for t in traces]
        X=np.concatenate([a for a,b in pairs]);y=np.concatenate([b for a,b in pairs])
        if method=='history_retrieval':
            for neighbors in (1,4,16):
                m=Model(p,norm,bank_x=X,bank_y=y,neighbors=neighbors)
                candidates.append(m)
            continue
        configurations=(0.,1e-3,1e-2,1e-1) if method=='sparse_history' else RIDGES
        gram=X.T@X/len(y);rhs=X.T@y/len(y)
        work['normal_equation_work_proxy']+=2*len(y)*X.shape[1]**2
        work['training_rows_processed']+=len(y)
        for value in configurations:
            lam=1e-5 if method=='sparse_history' else value
            penalty=np.eye(X.shape[1])*lam;penalty[-1,-1]=0
            coef=np.linalg.solve(gram+penalty,rhs)
            work['candidate_fits']+=1;work['normal_equation_work_proxy']+=X.shape[1]**3
            if method=='sparse_history':
                for _ in range(3):
                    keep=np.abs(coef)>=value;keep[-1]=True
                    small=np.linalg.solve((gram+penalty)[np.ix_(keep,keep)],rhs[keep])
                    coef=np.zeros(X.shape[1]);coef[keep]=small
                    work['refits']+=1;work['normal_equation_work_proxy']+=int(keep.sum())**3
            candidates.append(Model(p,norm,coef=coef))
    scored=[]
    for model in candidates:
        with np.errstate(over='ignore',invalid='ignore'):
            err=blocks(model,validation)
        work['selection_forecast_blocks']+=len(err)
        if model.bank_x is not None:
            work['retrieval_distance_coordinates']+=err.size*model.bank_x.size
        score=float(np.mean(err**2))
        if not np.isfinite(score):work['failed_candidates']+=1;score=float('inf')
        scored.append(score)
    best=int(np.argmin(scored))
    if not np.isfinite(scored[best]):raise ValueError('all diagnostic candidates failed')
    selected=candidates[best]
    if method=='privileged_history_state':compile_state(selected)
    work['selection_nrmse']=float(np.sqrt(scored[best]));work['fit_elapsed_seconds']=perf_counter()-began
    work['validation_refusal']=work['selection_nrmse']>.15
    selected.work=work
    return selected


def evaluate(model,calibration,test,method):
    began=perf_counter();compiled=method=='privileged_history_state'
    cal=blocks(model,calibration,compiled)
    maxima=np.max(np.abs(cal),axis=1)
    rank=int(np.ceil((len(maxima)+1)*.9))
    radius=float(np.sort(maxima)[rank-1]) if rank<=len(maxima) else float('inf')
    errors=blocks(model,test,compiled)
    nrmse=float(np.sqrt(np.mean(errors**2)))
    coverage=float(np.mean(np.max(np.abs(errors),axis=1)<=radius))
    p=model.p;k=len(model.norm.um);state=p+(p-1)*k+1
    coefficients=p*(1+k)+1
    if compiled:arithmetic=2*(state*state+state*k);storage=state*state+state*k+state
    elif model.bank_x is not None:arithmetic=3*model.bank_x.size;storage=model.bank_x.size+len(model.bank_y)+state
    else:arithmetic=2*coefficients;storage=coefficients+state
    return {'order':p,'nrmse':nrmse,'block_coverage':coverage,'calibration_radius_normalized':radius,
            'test_blocks':len(errors),'quality_pass':nrmse<=.15 and coverage>=.9 and not model.work['validation_refusal'],
            'forecast_arithmetic_proxy':int(arithmetic*errors.size),'retained_numeric_scalars':int(storage),
            'evaluation_elapsed_seconds':perf_counter()-began,'work':model.work,
            'signed_normalized_errors':errors.tolist(),
            'coefficients':None if model.coef is None else model.coef.tolist()}


def physical_system(seed):
    rng=np.random.default_rng(seed);n=int(rng.integers(3,7))
    weights=rng.uniform(.03,.18,n-1);ground=rng.uniform(.005,.025,n)
    K=np.diag(ground)
    for i,w in enumerate(weights):
        K[i,i]+=w;K[i+1,i+1]+=w;K[i,i+1]-=w;K[i+1,i]-=w
    # Exact discretization of a stable passive linear thermal network.
    values,vectors=np.linalg.eigh(K)
    A=(vectors*np.exp(-values))@vectors.T
    heaters=np.zeros((n,2));heaters[0,0]=.08;heaters[-1,1]=.08
    B=(vectors*((1-np.exp(-values))/values))@vectors.T@heaters
    regime='colored' if seed%2 else 'white'
    traces=[]
    for index in range(10):
        length=192;t=np.arange(length)
        if index<8:
            u=np.repeat(rng.uniform(-1,1,(24,2)),8,axis=0)
        elif index==8:
            u=np.column_stack((np.where((t//24)%2,1.,-1.),np.where((t//32)%2,.7,-.7)))
        else:
            u=np.column_stack((np.sin(2*np.pi*t/43),.8*np.sin(2*np.pi*t/67+.7)))
        x=rng.normal(0,.1,n);noise=0.;y=[float(x[0])]
        for drive in u:
            x=A@x+B@drive
            noise=(.75*noise+float(rng.normal(0,.002*np.sqrt(1-.75**2)))) if regime=='colored' else float(rng.normal(0,.002))
            y.append(float(x[0])+noise)
        traces.append(Trace(u,np.asarray(y)))
    fixture={'seed':seed,'dimension':n,'regime':regime,'generator_eigenvalues':values.tolist(),
             'trajectory_count':10,'simulator_transitions':1920,'new_physical_measurements':0}
    return traces[:4],traces[4:6],traces[6:8],traces[8:],fixture


def measured_traces():
    directory=ROOT/'data/v6'
    step=np.genfromtxt(directory/'tclab_step_test.csv',delimiter=',',names=True)
    sine=np.genfromtxt(directory/'tclab_sine_test_5min_period.csv',delimiter=',',names=True)
    def trace(table,lo,hi):
        return Trace(np.column_stack((table['Q1'][lo:hi],table['Q2'][lo:hi])),table['T1'][lo:hi+1])
    return [trace(step,0,540)],[trace(step,540,720)],[trace(step,720,900)],[trace(sine,0,900)]


def assess(train,selection,calibration,test,fixture):
    data_hash=sha256(b''.join(t.u.tobytes()+t.y.tobytes() for t in train+selection+calibration+test)).hexdigest()
    row={'fixture':fixture,'data_sha256':data_hash,'methods':{},
         'common_observation_records':sum(len(t.y) for t in train+selection+calibration+test),
         'common_input_transition_records':sum(len(t.u) for t in train+selection+calibration+test)}
    for method in METHODS:
        model=fit(train,selection,method)
        row['methods'][method]=evaluate(model,calibration,test,method)
    assert data_hash==sha256(b''.join(t.u.tobytes()+t.y.tobytes() for t in train+selection+calibration+test)).hexdigest(), 'an arm changed shared observations'
    a=row['methods']['privileged_history_state'];b=row['methods']['stateful_arx']
    row['privileged_vs_arx_max_forecast_difference_normalized']=float(np.max(np.abs(np.asarray(a['signed_normalized_errors'])-np.asarray(b['signed_normalized_errors']))))
    return row


def run():
    protocol=ROOT/'research/v6/PROTOCOL.md'
    frozen=json.loads((ROOT/'research/v6/protocol_freeze.json').read_text())
    assert sha256(protocol.read_bytes()).hexdigest()==frozen['sha256'],'protocol changed after freeze'
    rows=[]
    for seed in range(62001,62033):
        train,selection,calibration,test,fixture=physical_system(seed)
        rows.append(assess(train,selection,calibration,test,fixture))
    measured=assess(*measured_traces(),{'kind':'independently_measured_TCLab','independent_apparatus_count':1,
                      'new_physical_measurements':0,'active_intervention_savings_claim':False})
    summary={}
    for method in METHODS:
        records=[row['methods'][method] for row in rows]
        summary[method]={'mean_nrmse':float(np.mean([r['nrmse'] for r in records])),
                         'mean_block_coverage':float(np.mean([r['block_coverage'] for r in records])),
                         'quality_pass_systems':sum(r['quality_pass'] for r in records),
                         'validation_refusals':sum(r['work']['validation_refusal'] for r in records),
                         'candidate_fits':sum(r['work']['candidate_fits'] for r in records),
                         'refits':sum(r['work']['refits'] for r in records),
                         'normal_equation_work_proxy':sum(r['work']['normal_equation_work_proxy'] for r in records),
                         'forecast_arithmetic_proxy':sum(r['forecast_arithmetic_proxy'] for r in records),
                         'fit_seconds':sum(r['work']['fit_elapsed_seconds'] for r in records)}
    equivalence=max(r['privileged_vs_arx_max_forecast_difference_normalized'] for r in rows+[measured])
    if equivalence>1e-10:raise AssertionError('unexpected companion/recurrence discrepancy; investigate before interpretation')
    charged_headroom=all(r['methods']['privileged_history_state']['forecast_arithmetic_proxy']
                        <=.9*r['methods']['stateful_arx']['forecast_arithmetic_proxy']
                        and r['methods']['privileged_history_state']['work']['normal_equation_work_proxy']
                        <=r['methods']['stateful_arx']['work']['normal_equation_work_proxy'] for r in rows)
    return {'stage':'privileged positive-control headroom only','autonomous_acquisition_performed':False,
            'protocol_sha256':frozen['sha256'],'systems':rows,'measured_anchor':measured,'summary':summary,
            'maximum_companion_arx_difference':equivalence,
            'headroom_gate':'established' if charged_headroom else 'not_established',
            'reason':'equal predictive laws; assess repeatable arithmetic separately from implementation timing',
            'physical_measurements_new':0,'simulator_trajectories':320,'simulator_transitions':61440,
            'independent_simulated_systems':32,'independent_measured_apparatuses':1}


if __name__=='__main__':
    result=run();directory=ROOT/'results/v6';directory.mkdir(parents=True,exist_ok=True)
    (directory/'headroom_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('systems','measured_anchor')},indent=2))
