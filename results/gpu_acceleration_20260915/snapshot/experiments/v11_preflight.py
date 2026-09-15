"""Exact-polynomial, complete-cost paired matrix arithmetic preflight."""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from random import Random
from hashlib import sha256
import json
from .v7_headroom import model,TOL
from .v7_certificate import Generator,check_certificate,evaluate,clean
from .v8_integer_taylor import fraction_free_taylor
from .v11_matrix_taylor import matrix_taylor
from .v7_verify import packed
ROOT=Path(__file__).resolve().parents[1]


def calculation(arm,n=3,family='xxz',gamma=F(2),T=F(1,2)):
    start=perf_counter();h,o=model(n,family);g=Generator(h,gamma,n,512)
    try:
        method=matrix_taylor if arm=='matrix' else fraction_free_taylor
        p,w,c=method(g,o,T,TOL);constructed=perf_counter()
        receipt=check_certificate(Generator(h,gamma,n,512),o,[p],w,TOL,expected_time=T)
        if receipt['status']!='certified':raise ValueError('checker: '+str(receipt))
        checked=perf_counter();out=clean(evaluate(p.coefficients,T),n,512);done=perf_counter()
        row=dict(arm=arm,status='certified',bound=w['claimed_bound'],degree=len(p.coefficients)-1,
                 construction_seconds=constructed-start,checking_seconds=checked-constructed,
                 endpoint_seconds=done-checked,complete_seconds=done-start,cost=c,checker=receipt)
        artifact=packed(n,h,gamma,o,T,p,w,dict(arm=arm,family=family))
        artifact['output']={p:str(v) for p,v in out.items()}
        return row,artifact,p,w
    except ValueError as exc:
        return dict(arm=arm,status='refused',reason=str(exc),complete_seconds=perf_counter()-start,
                    cost=getattr(exc,'construction_cost',dict(g.cost))),None,None,None


def run():
    warm=[];by_arm={};archive=[];pairs=[]
    for arm in ('integer','matrix'):
        row,artifact,p,w=calculation(arm);warm.append(row)
        if row['status']!='certified':raise AssertionError(row)
        by_arm[arm]=(p,w)
        archive.append(artifact)
    if by_arm['integer'][0]!=by_arm['matrix'][0]:raise AssertionError('polynomial mismatch')
    if by_arm['integer'][1]['witnesses']!=by_arm['matrix'][1]['witnesses']:raise AssertionError('norm witness mismatch')
    rng=Random(11)
    for repetition in range(7):
        order=['integer','matrix'];rng.shuffle(order)
        for arm in order:
            row,_,p,w=calculation(arm)
            if row['status']!='certified' or p!=by_arm[arm][0]:raise AssertionError('paired output drift')
            row['repetition']=repetition;pairs.append(row)
            print(repetition,arm,round(row['complete_seconds'],6),flush=True)
    files=['research/v11/PREFLIGHT_PROTOCOL.md','experiments/v11_matrix_taylor.py','experiments/v11_preflight.py',
           'experiments/v8_integer_taylor.py','experiments/v7_certificate.py']
    result=dict(schema='v11-paired-1',scope='one existing development case; no learned method',
                warmup=warm,rows=pairs,exact_polynomial_equality=True,exact_norm_witness_equality=True,
                source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    (ROOT/'results/v11/preflight.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'results/v11/certificate_archive.json').write_text(json.dumps(archive,separators=(',',':'))+'\n')
    return result


if __name__=='__main__':run()
