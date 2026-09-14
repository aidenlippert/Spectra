"""Verify the interaction-contracted rule's exact equation-of-motion identity."""
from fractions import Fraction as F
from pathlib import Path
import json
import time

from experiments.marginal_symbolic import add,scale,product,mono
from research.molecular_collective_20260913.core import extract,retained_polynomial,factor_operators,digest


def check(data,tail):
    start=time.monotonic();p=extract(data,tail['center_number']);H=retained_polynomial(p,tail)
    factors=factor_operators(p,tail);checked=[]
    for mode in range(p['modes']):
        a=mono(((0,mode),));actual=add(product(H,a),scale(product(a,H),-1))
        expected=add(product(p['one'],a),scale(product(a,p['one']),-1))
        for weight,Q in factors:
            T={((0,j),):Q.get(((1,mode),(0,j)),F(0)) for j in range(p['modes']) if Q.get(((1,mode),(0,j)),0)}
            linear={}
            for j in range(p['modes']):
                c=sum(Q.get(((1,mode),(0,k)),F(0))*Q.get(((1,k),(0,j)),F(0)) for k in range(p['modes']))
                if c:linear[((0,j),)]=c
            expected=add(expected,scale(product(Q,T),-weight),scale(linear,-weight/2))
        residual=add(actual,scale(expected,-1))
        if residual:raise AssertionError('Equation-of-motion contraction identity failed')
        checked.append({'mode':mode,'commutator_terms':len(actual),'maximum_degree':max(map(len,actual),default=0)})
    return {'fixture_sha256':digest(data),'tail_sha256':digest(tail),'checked':checked,
        'residual_terms':0,'wall_seconds':time.monotonic()-start,
        'identity':'[Hret,a_i]=[h_tilde,a_i]-sum_k lambda_k Q_k (L_k a)_i - (1/2)sum_k lambda_k (L_k^2 a)_i',
        'scope':'Exact algebra for the alpha=1 spin-summed member before direction rounding. The alpha=1/2 variant and individual spin components are candidate extensions; this is not an energy or completeness theorem.'}


if __name__=='__main__':
    root=Path(__file__).resolve().parents[2];p=root/'results/molecular_collective_20260913/campaign/h6'
    result=check(json.loads((p/'fixture.json').read_text()),json.loads((p/'rank_10/tail.json').read_text()))
    (root/'results/mechanism_transfer_20260913/identity_check.json').write_text(json.dumps(result,indent=2)+'\n')
