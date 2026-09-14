"""Bounded exact Clifford-frame moment comparison and size probes."""
from pathlib import Path
import json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_clifford_moments import CliffordMomentOracle
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'clifford_moments';out.mkdir(exist_ok=True);rows=[]
    for sites in [8,12,16,20,24]:
        o=CliffordMomentOracle(build(sites,4,1))
        for order in [4,8,12,18,24]:
            start=time.monotonic()
            try:M,r=o.dimer_moments(order)
            except ValueError as error:
                row={'sites':sites,'order':order,'accepted':False,'reason':str(error),'source_frame_states':len(o.cache),'referenced_frame_states':len(o.referenced),'seconds':time.monotonic()-start}
            else:row={'sites':sites,'order':order,'accepted':True,'last_moment':str(M[-1][0][0]),'work':r,'seconds':time.monotonic()-start}
            rows.append(row);print(row,flush=True)
            (out/'scaling.json').write_text(json.dumps({'results':rows,'scope':'Single analytic dimer boundary; different orders are not equal-accuracy energy comparisons. Frame-source cap4096 and vector-support cap65536 remain fixed.'},indent=2)+'\n')
            if not row['accepted']:break
    e=json.loads((ROOT/'singlet_valence_embedding/certificate.json').read_text());V=[{int(s):a for s,a in v.items()} for v in e['basis']]
    h=json.loads((ROOT/'hamiltonian.json').read_text());o=CliffordMomentOracle(h);reference=json.loads((ROOT/'symmetry_moments24/moments.json').read_text());checks=[]
    for order in [4,8,12,14,18,24]:
        start=time.monotonic()
        try:
            M,r=o.moments(V,order)
            if M!=reference[:order+1]:raise ValueError('Exact H8 moment mismatch')
        except ValueError as error:row={'order':order,'accepted':False,'reason':str(error),'source_frame_states':len(o.cache),'seconds':time.monotonic()-start}
        else:row={'order':order,'accepted':True,'all_entries_match_signed_orbit_moments':True,'work':r,'seconds':time.monotonic()-start}
        checks.append(row);print('H8 all boundaries',row,flush=True)
        (out/'h8_comparison.json').write_text(json.dumps(checks,indent=2)+'\n')
        if not row['accepted']:break
if __name__=='__main__':main()
