"""Complete-cost test of a small exact cover versus the strongest baseline."""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from random import Random
from hashlib import sha256
from statistics import median
import json
from .v7_headroom import model,TOL
from .v7_certificate import Generator,check_certificate,evaluate,clean
from .v8_fractional_cover import check_evolution_cover
from .v11_guarded_taylor import guarded_taylor
from .v12_small_cover_taylor import small_cover_taylor
from .v7_verify import packed
ROOT=Path(__file__).resolve().parents[1]


def calculation(arm):
    start=perf_counter();h,o=model(3,'xxz');gamma=F(2);T=F(1,2)
    g=Generator(h,gamma,3,512)
    try:
        if arm=='baseline':p,w,c=guarded_taylor(g,o,T,TOL,guard='frobenius')
        elif arm in ('small_cover','disabled'):
            p,w,c=small_cover_taylor(g,o,T,TOL,enabled=arm=='small_cover')
        else:raise ValueError('unknown arm')
        constructed=perf_counter()
        checker=check_evolution_cover if w['schema']=='v8-evolution-cover-1' else check_certificate
        receipt=checker(Generator(h,gamma,3,512),o,[p],w,TOL,expected_time=T)
        if receipt['status']!='certified':raise AssertionError(receipt)
        checked=perf_counter();out=clean(evaluate(p.coefficients,T),3,512);done=perf_counter()
        row=dict(arm=arm,status='certified',bound=w['claimed_bound'],degree=len(p.coefficients)-1,
                 construction_seconds=constructed-start,checking_seconds=checked-constructed,
                 endpoint_seconds=done-checked,complete_seconds=done-start,cost=c,checker=receipt)
        artifact=packed(3,h,gamma,o,T,p,w,dict(arm=arm,family='xxz'))
        artifact['output']={p:str(v) for p,v in out.items()}
        return row,artifact,p,w,out
    except ValueError as exc:
        return dict(arm=arm,status='refused',reason=str(exc),complete_seconds=perf_counter()-start,
                    cost=getattr(exc,'construction_cost',dict(g.cost))),None,None,None,None


def run():
    warm=[];references={};archive=[];rows=[];arms=('baseline','small_cover','disabled')
    for arm in arms:
        row,artifact,p,w,out=calculation(arm);warm.append(row)
        if row['status']!='certified':raise AssertionError(row)
        references[arm]=(p,w,out);archive.append(artifact)
    if references['baseline'][0]!=references['disabled'][0] or references['baseline'][2]!=references['disabled'][2]:
        raise AssertionError('disabled baseline mismatch')
    if references['baseline'][1]['witnesses']!=references['disabled'][1]['witnesses']:
        raise AssertionError('disabled witness mismatch')
    rng=Random(1201)
    for repetition in range(21):
        order=list(arms);rng.shuffle(order)
        for arm in order:
            row,_,p,w,out=calculation(arm)
            ref=references[arm]
            if row['status']!='certified' or p!=ref[0] or out!=ref[2] or w['witnesses']!=ref[1]['witnesses']:
                raise AssertionError('paired output drift')
            row['repetition']=repetition;rows.append(row)
    base={r['repetition']:r for r in rows if r['arm']=='baseline'};stats={}
    for arm in arms:
        selected=[r for r in rows if r['arm']==arm]
        ratios=[r['complete_seconds']/base[r['repetition']]['complete_seconds'] for r in selected]
        stats[arm]=dict(total_seconds=sum(r['complete_seconds'] for r in selected),
            median_paired_ratio=median(ratios),wins=sum(x<1 for x in ratios),
            mean_construction_seconds=sum(r['construction_seconds'] for r in selected)/len(selected),
            mean_checking_seconds=sum(r['checking_seconds'] for r in selected)/len(selected),
            mean_endpoint_seconds=sum(r['endpoint_seconds'] for r in selected)/len(selected))
    s=stats['small_cover'];b=stats['baseline']
    passed=(s['total_seconds']<=.95*b['total_seconds'] and s['median_paired_ratio']<=.95 and s['wins']>=18)
    paths=['research/v12/PREFLIGHT_PROTOCOL.md','experiments/v12_preflight.py',
           'experiments/v12_small_cover.py','experiments/v12_small_cover_taylor.py',
           'experiments/v11_guarded_taylor.py','experiments/v7_certificate.py','experiments/v8_fractional_cover.py']
    result=dict(schema='v12-paired-1',scope='one existing development case; supplied conventional methods',
        warmup=warm,rows=rows,summary=stats,headroom_gate='passed' if passed else 'failed',
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    (ROOT/'results/v12/preflight.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'results/v12/certificate_archive.json').write_text(json.dumps(archive,separators=(',',':'))+'\n')
    print(json.dumps(dict(summary=stats,headroom_gate=result['headroom_gate']),indent=2),flush=True)
    return result


if __name__=='__main__':run()
