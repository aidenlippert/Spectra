from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from hashlib import sha256
import json
from experiments.v7_headroom import model,TOL
from experiments.v7_certificate import Generator,Piece,check_certificate,evaluate,clean
from experiments.v11_guarded_taylor import guarded_taylor
from experiments.v12_small_cover_taylor import small_cover_taylor
from experiments.v8_fractional_cover import check_evolution_cover

ROOT=Path(__file__).resolve().parents[1]

def main():
    h,o=model(3,'xxz'); T=F(1,2); gamma=F(2)
    p12,w12,base_work=guarded_taylor(Generator(h,gamma,3,512),o,T,TOL,max_order=12,guard='frobenius')
    p11,sw,small_work=small_cover_taylor(Generator(h,gamma,3,512),o,T,TOL,max_order=12,guard='frobenius')
    assert len(p12.coefficients)-1==12 and len(p11.coefficients)-1==11
    assert sw['schema']=='v8-evolution-cover-1'
    groups=len(sw['witnesses']['residual:0:11']['groups']); assert groups==14
    # Verify both actual certificates before timing.
    bcheck=check_certificate(Generator(h,gamma,3,512),o,[p12],w12,TOL,expected_time=T)
    scheck=check_evolution_cover(Generator(h,gamma,3,512),o,[p11],sw,TOL,expected_time=T)
    assert bcheck['status']=='certified',bcheck; assert scheck['status']=='certified',scheck
    baseline=[]; small=[]
    for i in range(7):
        g=Generator(h,gamma,3,512); t=perf_counter(); r=check_certificate(g,o,[p12],w12,TOL,expected_time=T); c=perf_counter()-t
        assert r['status']=='certified',r
        t=perf_counter(); clean(evaluate(p12.coefficients,T),3,512); e=perf_counter()-t
        baseline.append(dict(repetition=i,checker_seconds=c,endpoint_seconds=e,checker_cost=r))
        g=Generator(h,gamma,3,512); t=perf_counter(); r=check_evolution_cover(g,o,[p11],sw,TOL,expected_time=T); c=perf_counter()-t
        assert r['status']=='certified',r
        t=perf_counter(); clean(evaluate(p11.coefficients,T),3,512); e=perf_counter()-t
        small.append(dict(repetition=i,checker_seconds=c,endpoint_seconds=e,checker_cost=r))
    def avg(rows,k): return sum(x[k] for x in rows)/len(rows)
    files=['experiments/v11_guarded_taylor.py','experiments/v12_small_cover.py','experiments/v12_small_cover_taylor.py','experiments/v8_fractional_cover.py','experiments/v7_certificate.py','experiments/v7_headroom.py','experiments/v12_small_cover_floor_run.py']
    out={'schema':'v12-small-cover-floor-1','status':'certified','scope':'n3 XXZ gamma2 T=1/2 eps=.001; seven paired fresh checker plus exact endpoint operations','baseline_v11_guarded':{'method':'guarded_taylor','guard':'frobenius','degree':12,'bound':w12['claimed_bound'],'work':base_work},'small_cover_actual':{'method':'small_cover_taylor','guard':'frobenius','degree':11,'witness_schema':sw['schema'],'groups':groups,'incidences':small_work['cover_calls'][0]['work']['incidences'],'bound':sw['claimed_bound'],'work':small_work},'certification':{'baseline':bcheck,'small_cover':scheck},'old_61_group_witness':{'label':'oldwitness','used_for_ranking':False,'note':'retained as historical diagnostic only; no old cost reused or reranked'},'paired':{'baseline':baseline,'small_cover':small,'means':{'baseline_checker_seconds':avg(baseline,'checker_seconds'),'small_cover_checker_seconds':avg(small,'checker_seconds'),'checker_saved_seconds':avg(baseline,'checker_seconds')-avg(small,'checker_seconds'),'baseline_endpoint_seconds':avg(baseline,'endpoint_seconds'),'small_cover_endpoint_seconds':avg(small,'endpoint_seconds'),'endpoint_saved_seconds':avg(baseline,'endpoint_seconds')-avg(small,'endpoint_seconds'),'combined_saved_seconds':avg(baseline,'checker_seconds')+avg(baseline,'endpoint_seconds')-avg(small,'checker_seconds')-avg(small,'endpoint_seconds')}},'no_claims':['checker timing is not headroom','small-cover construction is excluded from paired floor timing','no held-out data'],'source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files}}
    (ROOT/'results/v12').mkdir(exist_ok=True); (ROOT/'research/v12').mkdir(exist_ok=True)
    (ROOT/'results/v12/small_cover_floor.json').write_text(json.dumps(out,indent=2)+'\n')
    return out
if __name__=='__main__':main()
