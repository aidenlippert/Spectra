"""Fixed-physics complete-cost reference expansion diagnostic."""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from random import Random
import hashlib,json
from .v7_headroom import model,TOL
from .v7_certificate import Generator,check_certificate,evaluate,clean
from .v8_integer_taylor import fraction_free_taylor
from .v9_galerkin import adaptive_galerkin
from .v10_reference import construct,check_expansion,evaluate_expansion
from .v7_verify import packed,replay
ROOT=Path(__file__).resolve().parents[1]
TRUNCATION=F(999,1000000)
READOUT=F(1,1000000)


def pack_reference(n,h,gamma,o,T,layers,output,bound):
    return dict(schema='v10-reference-1',n=n,h={p:str(c) for p,c in h.items()},
                gamma=str(gamma),initial={p:str(c) for p,c in o.items()},time=str(T),
                tolerance=str(TOL),truncation_tolerance=str(TRUNCATION),readout_tolerance=str(READOUT),
                layers=[[[p,str(a),str(b),k,str(c[0]),str(c[1])] for (p,a,b,k),c in d.items()] for d in layers],
                output={p:str(c) for p,c in output.items()},bound=str(bound))


def replay_reference(record):
    h={p:F(c) for p,c in record['h'].items()};o={p:F(c) for p,c in record['initial'].items()}
    gamma,T,n=F(record['gamma']),F(record['time']),record['n']
    layers=[]
    for rows in record['layers']:
        layer={(p,F(a),F(b),k):(F(cr),F(ci)) for p,a,b,k,cr,ci in rows}
        if len(layer)!=len(rows):raise ValueError('duplicate serialized mode')
        layers.append(layer)
    ans=check_expansion(h,gamma,n,o,T,F(record['truncation_tolerance']),layers)
    if ans['status']!='certified':raise ValueError('reference replay tail refused')
    output,error,cost=evaluate_expansion(layers,n,T,F(record['readout_tolerance']))
    total=F(ans['bound'])+error
    if total>F(record['tolerance']) or total!=F(record['bound']):raise ValueError('reference total bound mismatch')
    if output!={p:F(c) for p,c in record['output'].items()}:raise ValueError('reference output mismatch')
    return dict(status='certified',bound=str(total),order=ans['order'],entries=ans['total_entries'])


def attempt(n,family,gamma,T,arm):
    start=perf_counter();g=None;cost={};stage='initialization'
    try:
        h,o=model(n,family)
        stage='construction'
        if arm=='reference':
            layers,cost=construct(h,gamma,n,o,T,TRUNCATION)
            constructed=perf_counter()
            stage='checking';receipt=check_expansion(h,gamma,n,o,T,TRUNCATION,layers)
            if receipt['status']!='certified':raise ValueError('reference checker refused')
            checked=perf_counter()
            stage='readout';output,error,ec=evaluate_expansion(layers,n,T,READOUT)
            bound=F(receipt['bound'])+error
            if bound>TOL:raise ValueError('total tolerance exceeded')
            done=perf_counter()
            row=dict(status='certified',bound=str(bound),truncation_bound=receipt['bound'],readout_bound=str(error),
                     construction_seconds=constructed-start,checking_seconds=checked-constructed,
                     readout_seconds=done-checked,complete_seconds=done-start,
                     order=len(layers)-1,entries=sum(map(len,layers)),output_terms=len(output),
                     construction=cost,checker=receipt,readout=ec)
            artifact=pack_reference(n,h,gamma,o,T,layers,output,bound)
        else:
            g=Generator(h,gamma,n,512)
            if arm=='integer':
                p,w,cost=fraction_free_taylor(g,o,T,TOL)
                constructed=perf_counter();stage='checking'
                receipt=check_certificate(Generator(h,gamma,n,512),o,[p],w,TOL,expected_time=T)
            elif arm=='bfs':
                cost=adaptive_galerkin(g,o,T,TOL,growth='bfs')
                if cost['status']!='certified':raise ValueError(cost['reason'])
                p,w,receipt=cost.pop('piece'),cost.pop('certificate'),cost['checker']
                constructed=None  # Adaptive implementation includes its final check.
            else:raise ValueError('unknown arm')
            if receipt['status']!='certified':raise ValueError('polynomial checker refused')
            checked=perf_counter();stage='readout'
            output=clean(evaluate(p.coefficients,T),n,512)
            done=perf_counter()
            row=dict(status='certified',bound=w['claimed_bound'],readout_bound='0',
                     construction_and_checking_seconds=checked-start,readout_seconds=done-checked,
                     complete_seconds=done-start,order=len(p.coefficients)-1,
                     entries=sum(map(len,p.coefficients)),output_terms=len(output),cost=cost,checker=receipt)
            artifact=packed(n,h,gamma,o,T,p,w,dict(arm=arm,family=family))
            artifact['output']={p:str(c) for p,c in output.items()}
        return row,artifact
    except ValueError as exc:
        return dict(status='refused',reason=str(exc),stage=stage,complete_seconds=perf_counter()-start,
                    cost=getattr(exc,'construction_cost',cost),generator=dict(g.cost) if g else {}),None


def run():
    rows=[];archive=[];rng=Random(10)
    for n in (3,4,6):
        for family in ('xxz','mixed'):
            for gamma in (F(0),F(1,5),F(2)):
                for T in (F(1,5),F(1,2)):
                    arms=['integer','bfs','reference'];rng.shuffle(arms)
                    for arm in arms:
                        row,artifact=attempt(n,family,gamma,T,arm)
                        row.update(n=n,family=family,gamma=str(gamma),time=str(T),arm=arm);rows.append(row)
                        if artifact:archive.append(artifact)
                        print(n,family,gamma,T,arm,row['status'],round(row['complete_seconds'],3),row.get('reason',''),flush=True)
                        (ROOT/'results/v10/progress.json').write_text(json.dumps(rows,indent=2)+'\n')
    files=['research/v10/HEADROOM_PROTOCOL.md','experiments/v10_headroom.py','experiments/v10_reference.py',
           'experiments/v10_exp_readout.py','experiments/v9_galerkin.py','experiments/v8_integer_taylor.py','experiments/v7_certificate.py']
    output=dict(schema='v10-development-1',rows=rows,source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                scope='one development diagnostic; no autonomous acquisition; reserved evaluation unopened')
    (ROOT/'results/v10/headroom.json').write_text(json.dumps(output,indent=2)+'\n')
    (ROOT/'results/v10/certificate_archive.json').write_text(json.dumps(archive,separators=(',',':'))+'\n')
    print('complete',len(rows),'accepted',len(archive),flush=True)


if __name__=='__main__':run()
