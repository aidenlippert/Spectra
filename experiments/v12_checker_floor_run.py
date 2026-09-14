from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from hashlib import sha256
import json
from experiments.v7_headroom import model, TOL
from experiments.v7_certificate import Generator, Piece, residual_records, check_certificate, evaluate, clean
from experiments.v11_guarded_taylor import guarded_taylor
from experiments.v8_fractional_cover import check_evolution_cover, check_fractional_cover

ROOT=Path(__file__).resolve().parents[1]

def neg_cover(w):
    return dict(w, groups=[dict(g, coefficients={p:str(-F(c)) for p,c in g['coefficients'].items()}) for g in w['groups']])

def main():
    h,o=model(3,'xxz'); T=F(1,2); gamma=F(2)
    t0=perf_counter(); gb=Generator(h,gamma,3,512)
    p12,w12,work=guarded_taylor(gb,o,T,TOL,max_order=12,guard='frobenius')
    construction_total=perf_counter()-t0
    p11=Piece(T,p12.coefficients[:12])
    # Independently derive -G c11 and verify against the known row-10 residual.
    dg=Generator(h,gamma,3,512); actual=residual_records(dg,o,[p11],'power')[-1][1]
    old=json.loads((ROOT/'results/v8/fractional_probe.json').read_text())['rows'][10]
    oldres={p:F(c) for p,c in old['residual'].items()}
    assert actual and oldres=={p:-c for p,c in actual.items()}
    cover=neg_cover(old['cover'])
    empty={'schema':'v8-cover-2','groups':[],'claimed_bound':'0'}
    witnesses={'jump:0':empty}
    for k in range(11): witnesses[f'residual:0:{k}']=empty
    witnesses['residual:0:11']=cover
    evolution={'schema':'v8-evolution-cover-1','integration_basis':'power','witnesses':witnesses,
               'claimed_bound':old['integrated_bound']}
    # Independent exact evolution check once before timing pairs.
    verify=check_evolution_cover(Generator(h,gamma,3,512),o,[p11],evolution,TOL,expected_time=T)
    assert verify['status']=='certified',verify
    # Baseline checker and exact endpoint are freshly instantiated each repetition.
    baseline=[]; overlap=[]
    for i in range(7):
        g=Generator(h,gamma,3,512); t=perf_counter(); r=check_certificate(g,o,[Piece(T,p12.coefficients)],w12,TOL,expected_time=T); c=perf_counter()-t
        assert r['status']=='certified',r
        t=perf_counter(); endpoint=clean(evaluate(p12.coefficients,T),3,512); e=perf_counter()-t
        baseline.append(dict(repetition=i,checker_seconds=c,endpoint_seconds=e,checker_cost=r))
        g=Generator(h,gamma,3,512); t=perf_counter(); r=check_evolution_cover(g,o,[p11],evolution,TOL,expected_time=T); c=perf_counter()-t
        assert r['status']=='certified',r
        t=perf_counter(); endpoint=clean(evaluate(p11.coefficients,T),3,512); e=perf_counter()-t
        overlap.append(dict(repetition=i,checker_seconds=c,endpoint_seconds=e,checker_cost=r))
    def avg(rows,k): return sum(x[k] for x in rows)/len(rows)
    files=['experiments/v11_guarded_taylor.py','experiments/v8_fractional_cover.py','experiments/v7_certificate.py','experiments/v7_headroom.py','experiments/v12_checker_floor_run.py','results/v8/fractional_probe.json']
    out={'schema':'v12-checker-floor-1','status':'certified','scope':'already-public n3 XXZ gamma2 T=1/2 eps=.001; seven paired fresh checker plus exact endpoint operations','baseline':{'method':'experiments.v11_guarded_taylor.guarded_taylor','guard':'frobenius','max_order':12,'degree':12,'bound':w12['claimed_bound'],'construction_seconds_total':construction_total,'construction_work':work,'coefficient_term_counts':[len(c) for c in p12.coefficients]},'accepted_v8_overlap':{'degree':11,'residual_definition':'actual -G c11 via Generator.apply / residual_records','row10_sign_relation':'row10 residual = - actual residual; cover coefficients negated','integrated_bound':old['integrated_bound'],'evolution_verification':verify},'paired':{'baseline':baseline,'overlap':overlap,'means':{'baseline_checker_seconds':avg(baseline,'checker_seconds'),'overlap_checker_seconds':avg(overlap,'checker_seconds'),'checker_saved_seconds':avg(baseline,'checker_seconds')-avg(overlap,'checker_seconds'),'baseline_endpoint_seconds':avg(baseline,'endpoint_seconds'),'overlap_endpoint_seconds':avg(overlap,'endpoint_seconds'),'endpoint_saved_seconds':avg(baseline,'endpoint_seconds')-avg(overlap,'endpoint_seconds'),'combined_saved_seconds':avg(baseline,'checker_seconds')+avg(baseline,'endpoint_seconds')-avg(overlap,'checker_seconds')-avg(overlap,'endpoint_seconds')}},'no_claims':['checking floor is not headroom','candidate construction cost for overlap proposer is excluded from paired floor timing','no held-out data'],'source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files}}
    (ROOT/'results/v12').mkdir(exist_ok=True)
    (ROOT/'research/v12').mkdir(exist_ok=True)
    (ROOT/'results/v12/checker_floor.json').write_text(json.dumps(out,indent=2)+'\n')
    return out
if __name__=='__main__': main()
