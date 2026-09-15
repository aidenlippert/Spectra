"""Local-tensor certificate for a known supersymmetric open t=V chain.

The qbar construction is the spin reversal of Hagendorf--Lienardy,
arXiv:1612.02951, equations 4--12 and 16. This is a narrow structural
recognizer, not a new general theorem or the main fixed-M CAR/SOS format.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import time

from research.side_routes_20260913.positive_chain import rational


def transpose(a): return [list(row) for row in zip(*a)]


def multiply(a,b):
    if len(a[0]) != len(b): raise ValueError('Incompatible tensor shapes')
    return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*b)] for row in a]


def kron(a,b):
    return [[x*y for x in row_a for y in row_b] for row_a in a for row_b in b]


def linear(*terms):
    rows,cols=len(terms[0][1]),len(terms[0][1][0])
    if any(len(a)!=rows or any(len(row)!=cols for row in a) for _,a in terms):
        raise ValueError('Incompatible tensor sum')
    return [[sum((c*a[i][j] for c,a in terms),F(0)) for j in range(cols)] for i in range(rows)]


def exact_matrix(value,rows,cols):
    if not isinstance(value,list) or len(value)!=rows or any(not isinstance(row,list) or len(row)!=cols for row in value):
        raise ValueError('Invalid local tensor shape')
    return [[rational(x) for x in row] for row in value]


def parameters(model):
    allowed={'modes','particles','hopping','interaction','fields','offset','energy_unit'}
    if set(model)-allowed: raise ValueError('Unsupported Hamiltonian fields')
    m,n=model['modes'],model['particles']
    if type(m) is not int or not 2<=m<=1000000 or type(n) is not int or not 0<=n<=m:
        raise ValueError('Invalid chain sector or recognition budget')
    data=[]
    for key,length in (('hopping',m-1),('interaction',m-1),('fields',m)):
        raw=model[key]
        if not isinstance(raw,list) or len(raw)!=length: raise ValueError('Invalid Hamiltonian array')
        data.append([rational(x) for x in raw])
    t,v,h=data; scale,mu=t[0],h[0]
    if scale<=0 or any(x!=scale for x in t+v) or any(x!=mu for x in h):
        raise ValueError('Hamiltonian does not match the uniform positive t=V structural class')
    return m,n,scale,mu,rational(model.get('offset',0))


def discover(model):
    m,n,scale,mu,offset=parameters(model)
    lower=offset+mu*n-scale*(m-n)
    return {'q':[[0,0],[0,0],[0,0],[1,0]],'cycle_pair':[[0],[1],[1],[0]],
            'dual_pair':[0,1],'lower':str(lower),
            'upper':str(lower) if m%2==0 and n==m//2 else None,
            'scope':'Known length-changing supersymmetry identity; exact saturation certified only for even half-filled chains.'}


def verify(model,witness):
    start=time.monotonic()
    m,n,scale,mu,offset=parameters(model)
    q=exact_matrix(witness['q'],4,2); qd=transpose(q); identity=[[F(1),F(0)],[F(0),F(1)]]
    split_left,split_right=kron(q,identity),kron(identity,q)
    coassoc=multiply(linear((1,split_left),(-1,split_right)),q)
    if any(x for row in coassoc for x in row): raise ValueError('Local nilpotency identity failed')
    number_change=[i.bit_count()-j for i in range(4) for j in range(2) if q[i][j]]
    if number_change!=[2]: raise ValueError('Wrong particle grading of local supercharge')
    qdq=multiply(qd,q)
    density=linear((-1,multiply(kron(identity,qd),split_left)),
                   (-1,multiply(kron(qd,identity),split_right)),
                   (1,multiply(q,qd)),(F(1,2),kron(qdq,identity)),(F(1,2),kron(identity,qdq)))
    expected=[[F(1),0,0,0],[0,F(1,2),-1,0],[0,-1,F(1,2),0],[0,0,0,F(1)]]
    boundary=linear((F(1,2),qdq))
    if density!=expected or boundary!=[[F(1,2),0],[0,0]]:
        raise ValueError('Local tensor does not reproduce the declared bulk and boundaries')
    lower=offset+mu*n-scale*(m-n)
    if rational(witness['lower'])!=lower: raise ValueError('False structural lower claim')
    upper=None
    if witness.get('upper') is not None:
        if m%2 or n!=m//2: raise ValueError('Saturation witness requires even half filling')
        chi=exact_matrix(witness['cycle_pair'],4,1)
        if any(chi[i][0] and i.bit_count()!=1 for i in range(4)):
            raise ValueError('Cycle has the wrong fixed-particle grading')
        if any(x[0] for x in multiply(linear((-1,split_left),(1,split_right)),chi)):
            raise ValueError('Two-site cycle is not closed')
        pair=witness['dual_pair']
        if pair!=[0,1]: raise ValueError('Unsupported alternating dual pattern')
        # Both neighboring patterns occur in the repeated 01 product state.
        for code in (1,2):
            basis=[[F(int(i==code))] for i in range(4)]
            if any(x[0] for x in multiply(qd,basis)):
                raise ValueError('Dual product state is not co-closed')
        if chi[1][0]!=1: raise ValueError('Cycle/dual overlap is not the certified nonzero value')
        if rational(witness['upper'])!=lower: raise ValueError('False saturation claim')
        upper=lower
    return {'method':'known_dynamic_supersymmetry_local_tensor_v1','modes':m,'particles':n,
            'lower':str(lower),'upper':str(upper) if upper is not None else None,
            'width':'0' if upper is not None else None,'model_sha256':hashlib.sha256(json.dumps(model,sort_keys=True).encode()).hexdigest(),
            'local_tensor_shape':[4,2],'largest_local_identity_shape':[8,4],
            'hamiltonian_coefficients_checked':3*m-2,'many_body_states_enumerated':0,
            'witness_json_bytes':len(json.dumps(witness).encode()),'wall_seconds':time.monotonic()-start,
            'scope':'Exact recognized-class certificate using a length-changing complex. No generic chemistry or fixed-M CAR/SOS claim.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    cases=[]
    for m in (4,8,12,16,32,64,256,1024,10000):
        model={'modes':m,'particles':m//2,'hopping':['1']*(m-1),'interaction':['1']*(m-1),'fields':['0']*m,
               'energy_unit':'t_reference=1; abstract lattice model'}
        start=time.monotonic();witness=discover(model);result=verify(model,witness)
        cases.append({'model':model,'witness':witness,'claim':result,'discovery_and_verification_seconds':time.monotonic()-start})
    args.out.write_text(json.dumps({'cases':cases,'known_source':'https://arxiv.org/html/1612.02951'},indent=2)+'\n')
    print(json.dumps([{'modes':a['model']['modes'],'energy':a['claim']['lower'],'width':a['claim']['width'],'seconds':a['discovery_and_verification_seconds']} for a in cases]))
