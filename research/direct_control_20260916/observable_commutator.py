"""Exact CAR coefficient bounds; does not enumerate any state sector."""
from fractions import Fraction as F
import json
from pathlib import Path
from experiments.marginal_symbolic import decode
from experiments.marginal_hunt_car import add, mul, scale


def commutator(a,b):
    return add(mul(a,b),scale(mul(b,a),-1))


def calculate(data):
    h=decode(data['hamiltonian'],data['modes'],4)
    d={((1,p),(0,p)):F(c) for p,c in [(6,1),(7,1),(10,-1),(11,-1)]}
    w={((1,p),(0,q)):F(1) for p,q in [(6,10),(10,6),(7,11),(11,7)]}
    hd=commutator(h,d); hw=commutator(h,w)
    cd=sum(map(abs,hd.values()),F()); cw=sum(map(abs,hw.values()),F())
    phases=[(F(-88687,1000000),F(1,2)),(F(1,2),F(1,2)),(F(1,2),F(1,2)),(F(311663,1000000),F(1,2))]
    # ||(H-E) Uc(t) psi|| <= sigma_initial + integral_0^t ||[H,C(s)]|| ds.
    # Hence integrated defect <= T*sigma_initial + integral_0^T (T-s)c(s) ds.
    tau=F(1,2); T=4*tau
    cb=F()
    for i,(u,v) in enumerate(phases):
        cb+=(abs(u)*cd+abs(v)*cw)*(T*tau-F(2*i+1,2)*tau*tau)
    return {'status':'coefficient_norm_bound_only','HD_terms':len(hd),'HW_terms':len(hw),
            'HD_operator_norm_upper_Ha':str(cd),'HW_operator_norm_upper_Ha':str(cw),
            'control_commutator_integrated_defect_allowance':str(cb),
            'allowance_float':float(cb),'initial_variance_additional_term':'2 * sigma_initial',
            'enumerated_configurations':0,'conclusion':'This sufficient bound is inconclusive, not a family impossibility proof.'}


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('--output',required=True);a=p.parse_args()
    result=calculate(json.loads(Path(a.fixture).read_text()))
    out=Path(a.output)
    if out.exists(): raise FileExistsError(out)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
