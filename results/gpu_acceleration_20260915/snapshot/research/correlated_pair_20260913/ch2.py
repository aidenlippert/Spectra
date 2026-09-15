"""Exact target-spin upper inferred from a charge MPS and bounded leakage."""
from fractions import Fraction as F
import json
from pathlib import Path
import time
from experiments.marginal_symbolic import decode,product,add,scale
from research.ch2_validation_20260913.certify import spin_square,check as old_check
from research.correlated_pair_20260913.mps_exact import State,check


def projected_upper(data,cert,spin):
    if type(spin) is not int or spin not in (0,1) or cert['spin_counts']!=[3+spin,3-spin]:raise ValueError('Declared CH2 target and M_S required')
    basic=check(data,cert);state=State(data,cert);h=decode(data['hamiltonian'],12,4);s2=spin_square(12)
    if add(product(h,s2),scale(product(s2,h),-1)):raise ValueError('H does not exactly commute with S squared')
    U=F(basic['upper_Ha']);svalue=state.expectation(s2);excess=svalue-spin*(spin+1);gap=2*(spin+1)
    if excess<0:raise ValueError('Impossible S squared expectation in exact M_S sector')
    leakage=excess/gap
    if leakage>=1:raise ValueError('No positive target-spin weight certified')
    ell=h.get((),F(0))-sum(abs(c) for w,c in h.items() if w)
    if U<ell:raise ValueError('Inconsistent global bound')
    # Let w be target-spin weight. U >= w E_target+(1-w)ell.
    # w >= 1-leakage > 0, so E_target <= ell+(U-ell)/(1-leakage).
    bound=ell+(U-ell)/(1-leakage)
    return {'raw_MPS':basic,'S_squared_expectation':str(svalue),'target_spin':spin,'target_spin_weight_lower':str(1-leakage),
            'target_spin_projected_upper_Ha':str(bound),'upper_float_Ha':float(bound),'projection_allowance_mHa':float((bound-U)*1000),
            'global_coefficient_L1_lower_Ha':str(ell),'exact_H_commutes_S2':True,
            'spin_purity_claim':'The implicit nonzero target-spin projection has this upper; the rounded unprojected MPS is not claimed spin pure.'}


def run():
    start=time.monotonic();base=Path('results/correlated_pair_20260913/ch2');data=json.loads(Path('results/ch2_validation_20260913/fixture.json').read_text());rows=[]
    for spin,name in ((0,'singlet'),(1,'triplet')):
        cert=json.loads((base/name/'state.json').read_text());new=projected_upper(data,cert,spin)
        old=old_check(data,json.loads(Path(f'results/ch2_validation_20260913/spin_{spin}_certificate.json').read_text()))
        L=F(old['lower_Ha']);U=F(new['target_spin_projected_upper_Ha'])
        if U<L:raise ValueError('Spin endpoints inconsistent')
        rows.append({'spin':spin,'lower_Ha':str(L),'upper_Ha':str(U),'width_mHa':float((U-L)*1000),'new_upper':new,'old_enumerated_lower_control':old})
    S,T=rows;lo=F(S['lower_Ha'])-F(T['upper_Ha']);hi=F(S['upper_Ha'])-F(T['lower_Ha'])
    result={'states':rows,'gap_convention':'E_S-E_T','gap_lower_Ha':str(lo),'gap_upper_Ha':str(hi),
            'gap_lower_mHa':float(lo*1000),'gap_upper_mHa':float(hi*1000),'gap_width_mHa':float((hi-lo)*1000),
            'triplet_below_singlet_certified_for_model':lo>0,'wall_seconds':time.monotonic()-start,
            'scope':'Frozen-core finite-basis fixed-geometry rational CH2 Hamiltonian. Inherited lower uses enumerated spin-projection matrices. No experimental or model-error claim.'}
    (base/'interval.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='states'},indent=2))


if __name__=='__main__':run()
