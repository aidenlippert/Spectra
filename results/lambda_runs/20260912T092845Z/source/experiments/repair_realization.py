"""Data-only exact linear predictive realization; established baseline method.

No physical evaluator, hidden state, system matrix or intended rank is imported.
Model-class assumptions are explicit and are not inferred from the data.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
from .v7_identifiability import _rref,_solve
from .v7_certificate import proof_fraction


def q(x):
    if isinstance(x,float) or isinstance(x,bool):raise ValueError('exact data required')
    v=proof_fraction(x)
    if max(abs(v.numerator).bit_length(),v.denominator.bit_length())>256:
        raise ValueError('input rational budget')
    return v


def mv(a,x):return [sum((v*w for v,w in zip(row,x)),F(0)) for row in a]
def mm(a,b):return [[sum((u*v for u,v in zip(row,col)),F(0)) for col in zip(*b)] for row in a]
def transpose(a):return [list(c) for c in zip(*a)]
def rank(a):return len(_rref(a)[1]) if a else 0


def inverse(a):
    n=len(a);columns=[]
    if rank(a)!=n:raise ValueError('singular selected basis')
    for j in range(n):columns.append(_solve(a,[F(i==j) for i in range(n)]))
    return transpose(columns)


def parse_evidence(data):
    if not isinstance(data,dict) or not {'maximum_dimension','output_ports','records'}<=data.keys():
        raise ValueError('evidence mapping requires dimension, output ports and records')
    D=data['maximum_dimension'];p=data['output_ports']
    if type(D) is not int or not 1<=D<=8 or type(p) is not int or not 1<=p<=4:
        raise ValueError('dimension/port budget')
    if data.get('prior')!='zero_reset_strictly_proper_LTI':raise ValueError('unsupported model prior')
    records=data['records']
    if not isinstance(records,list) or not 1<=len(records)<=128:raise ValueError('record budget')
    X=[];Y=[]
    for record in records:
        if not isinstance(record,dict) or not {'controls','output_intervals'}<=record.keys():
            raise ValueError('record mapping requires controls and output intervals')
        if record.get('reset_confirmed') is not True:raise ValueError('unknown initialization')
        us=record['controls'];intervals=record['output_intervals']
        if not isinstance(us,list) or not 1<=len(us)<=2*D:raise ValueError('history budget')
        us=[q(u) for u in us]
        if any(not 0<=u<=1 for u in us):raise ValueError('control outside interface')
        if not isinstance(intervals,list) or len(intervals)!=p:raise ValueError('port shape')
        values=[]
        for pair in intervals:
            if not isinstance(pair,list) or len(pair)!=2:raise ValueError('interval shape')
            lo,hi=map(q,pair)
            if lo!=hi:raise ValueError('uncertain measurement: exact realization refused')
            values.append(lo)
        X.append(list(reversed(us))+[F(0)]*(2*D-len(us)));Y.append(values)
    return D,p,X,Y


def learn(data):
    D,p,X,Y=parse_evidence(data);n=2*D
    rr,piv=_rref(X);solutions=[]
    for j in range(p):
        v=_solve(X,[row[j] for row in Y])
        if v is None:raise ValueError('data inconsistent with declared linear input-output law')
        solutions.append(v)
    identified=[]
    for k in range(n):
        weights=_solve(transpose(X),[F(i==k) for i in range(n)])
        if weights is not None:identified.append(k)
    partial={str(k):[str(s[k]) for s in solutions] for k in identified}
    if len(identified)<n:
        return dict(status='needs_evidence',identified_markov_vectors=partial,
            experiment=dict(reset='declared_zero_reset',controls=['1']+['0']*(n-1),
                            observe_each_step=True,reason='complete 2D moment bound, not inferred rank saturation'),
            class_assumptions=data['prior'],maximum_dimension=D)
    g=transpose(solutions)
    rows=[(j,k) for j in range(p) for k in range(D)]
    H=[[g[k+l][j] for l in range(D)] for j,k in rows]
    selected=[];basis=[]
    for i,row in enumerate(H):
        if rank(basis+[row])>len(basis):selected.append(i);basis.append(row)
    r=len(basis)
    if r:
        columns=_rref(basis)[1]
        K=[[H[i][j] for j in columns] for i in selected];inv=inverse(K)
        shift=[[g[rows[i][1]+l+1][rows[i][0]] for l in columns] for i in selected]
        M=mm(shift,inv)
        C=mm([[g[l][j] for l in columns] for j in range(p)],inv)
        B=[g[rows[i][1]][rows[i][0]] for i in selected]
    else:columns=[];M=[];C=[[] for _ in range(p)];B=[]
    model=dict(schema='observable_realization_1',maximum_dimension=D,output_ports=p,
        rank=r,selected_tests=[list(rows[i]) for i in selected],selected_columns=columns,
        transition=M,input=B,output=C,reset=[F(0)]*r,
        observed_moments=g,assumptions='exact zero-reset LTI, order bound and observed ports',
        acquisition=dict(records=len(X),scalar_outputs=len(X)*p,data_rank=len(piv),
            hankel_rows=len(H),hankel_columns=D,rank_tests=len(H)))
    verify(model)
    return dict(status='constructed',model=model)


def verify(model):
    D=model['maximum_dimension'];p=model['output_ports'];r=model['rank']
    if type(D) is not int or not 1<=D<=8 or type(p) is not int or not 1<=p<=4 or type(r) is not int or not 0<=r<=D:
        raise ValueError('model dimensions')
    M=model['transition'];B=model['input'];C=model['output'];g=model['observed_moments']
    if len(M)!=r or any(len(row)!=r for row in M) or len(B)!=r or len(C)!=p or any(len(row)!=r for row in C):raise ValueError('model shape')
    if len(g)!=2*D or any(len(row)!=p for row in g):raise ValueError('insufficient moment certificate')
    v=list(B)
    for k in range(2*D):
        if mv(C,v)!=g[k]:raise ValueError('moment reproduction failed')
        v=mv(M,v)
    if model['reset']!=[F(0)]*r:raise ValueError('invalid reset')
    return dict(status='exact_finite_order_certificate',moment_vectors=2*D,rank=r,
                universal_scope='conditional on true order<=D and exact LTI reset/port assumptions')


def verify_against_records(model,data):
    """Bind candidate moments to independently supplied observable records."""
    D,p,X,Y=parse_evidence(data)
    if D!=model['maximum_dimension'] or p!=model['output_ports']:raise ValueError('prior mismatch')
    verify(model)
    if rank(X)!=2*D:raise ValueError('evidence does not identify the required moments')
    if mm(X,model['observed_moments'])!=Y:raise ValueError('candidate does not reproduce observable evidence')
    g=model['observed_moments']
    H=[[g[k+l][j] for l in range(D)] for j in range(p) for k in range(D)]
    if rank(H)!=model['rank']:raise ValueError('nonminimal or incorrect rank claim')
    return dict(status='verified_against_observable_records',records=len(X),scalar_outputs=len(X)*p,
                identified_moment_vectors=2*D,rank=model['rank'])


def state_after(model,history,*,reset_confirmed=False):
    if reset_confirmed is not True or history is None:raise ValueError('recorded history from known reset required')
    if not isinstance(history,(tuple,list)) or len(history)>256:raise ValueError('initialization history budget')
    z=list(model['reset'])
    for u in history:z=step(model,z,q(u))
    return z


def step(model,z,u):
    if len(z)!=model['rank']:raise ValueError('state dimension')
    u=q(u)
    if not 0<=u<=1:raise ValueError('input outside declared interface')
    return [a+b*u for a,b in zip(mv(model['transition'],z),model['input'])]


def pack(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,dict):return {k:pack(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [pack(x) for x in v]
    return v


def unpack_model(model):
    result=dict(model)
    for key in ('transition','output','observed_moments'):result[key]=[[q(x) for x in row] for row in model[key]]
    for key in ('input','reset'):result[key]=[q(x) for x in model[key]]
    verify(result);return result


if __name__=='__main__':
    try:
        evidence=json.loads(Path(sys.argv[1]).read_text())
        answer=pack(learn(evidence))
    except (ValueError,TypeError) as exc:answer=dict(status='refused',reason=str(exc))
    Path(sys.argv[2]).write_text(json.dumps(answer,indent=2)+'\n')
