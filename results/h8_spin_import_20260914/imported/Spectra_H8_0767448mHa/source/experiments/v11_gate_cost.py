"""Fixed-workload paired study of conventional lower-norm rejection guards."""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from random import Random
from hashlib import sha256
import json
from .v7_headroom import model,TOL
from .v7_certificate import Generator,check_certificate,evaluate,clean
from .v8_integer_taylor import fraction_free_taylor
from .v11_guarded_taylor import guarded_taylor
from .v7_verify import packed
ROOT=Path(__file__).resolve().parents[1]


def attempt(n,family,gamma,T,arm):
    start=perf_counter();h,o=model(n,family);g=Generator(h,gamma,n,512)
    try:
        p,w,c=(fraction_free_taylor(g,o,T,TOL) if arm=='integer' else guarded_taylor(g,o,T,TOL,guard=arm))
        constructed=perf_counter()
        receipt=check_certificate(Generator(h,gamma,n,512),o,[p],w,TOL,expected_time=T)
        if receipt['status']!='certified':raise ValueError('checker refusal: '+str(receipt))
        checked=perf_counter();out=clean(evaluate(p.coefficients,T),n,512);done=perf_counter()
        row=dict(status='certified',bound=w['claimed_bound'],degree=len(p.coefficients)-1,
                 construction_seconds=constructed-start,checking_seconds=checked-constructed,
                 endpoint_seconds=done-checked,complete_seconds=done-start,cost=c,checker=receipt)
        artifact=packed(n,h,gamma,o,T,p,w,dict(arm=arm,family=family));artifact['output']={p:str(v) for p,v in out.items()}
        signature=(p,w['witnesses'],w['claimed_bound'],out)
        return row,artifact,signature
    except ValueError as exc:
        return dict(status='refused',reason=str(exc),complete_seconds=perf_counter()-start,
                    cost=getattr(exc,'construction_cost',dict(g.cost))),None,None


def run():
    rows=[];archive=[];signatures={};statuses={};rng=Random(111)
    for repetition in range(5):
        for n in (3,4,6):
            for family in ('xxz','mixed'):
                for gamma in (F(0),F(1,5),F(2)):
                    for T in (F(1,5),F(1,2)):
                        key=(n,family,str(gamma),str(T));arms=['integer','max','frobenius','cascade'];rng.shuffle(arms)
                        for arm in arms:
                            row,artifact,signature=attempt(n,family,gamma,T,arm)
                            if key in statuses and row['status']!=statuses[key]:raise AssertionError('gate changed success/refusal outcome')
                            statuses[key]=row['status']
                            if signature is not None:
                                if key in signatures and signature!=signatures[key]:raise AssertionError('gate changed accepted calculation')
                                signatures[key]=signature
                                if repetition==0:archive.append(artifact)
                            row.update(repetition=repetition,n=n,family=family,gamma=str(gamma),time=str(T),arm=arm);rows.append(row)
                        print(repetition,n,family,gamma,T,row['status'],flush=True)
    files=['research/v11/REJECTION_GATE_PROTOCOL.md','experiments/v11_gate_cost.py','experiments/v11_guarded_taylor.py',
           'experiments/v8_integer_taylor.py','experiments/v7_certificate.py']
    result=dict(schema='v11-rejection-guard-1',rows=rows,source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                matched_physical_inputs=len(statuses),exact_accepted_output_groups=len(signatures),
                exact_polynomial_norm_and_endpoint_equality=True,scope='supplied conventional baseline; no acquired methods or heldout evaluation')
    (ROOT/'results/v11/gate_cost.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'results/v11/gate_certificate_archive.json').write_text(json.dumps(archive,separators=(',',':'))+'\n')
    print('complete',len(rows),'archived',len(archive),flush=True)


if __name__=='__main__':run()
