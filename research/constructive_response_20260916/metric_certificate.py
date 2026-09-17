"""Exact energy-weighted response certificate and a four-site model checker.

The checker explicitly enumerates 36 balanced four-site configurations. It is
a correctness experiment, not a determinant-free many-body constructor.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
from math import isqrt

from research.constructive_response_20260916.composition_exact import (
    add,sub,sc,mm,tr,eye,block,is_psd,
)
from research.constructive_response_20260916.bond_coercivity import verify_bond_all_sectors


def square_model():
    """Independently generate exact t=1,U=8 four-site square Hamiltonian."""
    labels=sorted(sum(1<<(2*i) for i in up)+sum(1<<(2*i+1) for i in down)
                  for up in combinations(range(4),2) for down in combinations(range(4),2))
    index={x:i for i,x in enumerate(labels)}
    H=[[F(0) for _ in labels] for _ in labels]
    edges=[(0,1),(2,3),(0,2),(1,3)]
    for j,label in enumerate(labels):
        H[j][j]=F(8*sum((label>>(2*i))&3==3 for i in range(4)))
        for a,b in edges:
            for spin in (0,1):
                for dst,src in ((2*a+spin,2*b+spin),(2*b+spin,2*a+spin)):
                    if label&(1<<src) and not label&(1<<dst):
                        lo,hi=sorted((dst,src))
                        between=((1<<hi)-1)^((1<<(lo+1))-1)
                        sign=-1 if (label&between).bit_count()%2 else 1
                        out=label^(1<<src)^(1<<dst)
                        H[index[out]][j]-=sign
    if H != tr(H):
        raise AssertionError('Independent model is not symmetric')
    P=[i for i,x in enumerate(labels) if all((x>>(2*j))&3!=3 for j in range(4))]
    Q=[i for i in range(len(labels)) if i not in P]
    return labels,H,P,Q


def energy_bound(A,B,D,X,delta,rho,verify_gap=True):
    """Return eta with block H-e >= -eta I under checked exact premises."""
    delta=F(delta);rho=F(rho)
    if delta <= 0 or rho < 0 or 2*rho >= delta:
        raise ValueError('Need delta>2rho>=0')
    if verify_gap and not is_psd(sub(D,sc(delta,eye(len(D))))):
        raise ValueError('Eliminated gap premise failed')
    if D != tr(D) or A != tr(A):
        raise ValueError('Hermitian blocks required')
    R=sub(B,mm(D,X))
    M=add(eye(len(A)),mm(tr(X),X))
    K=add(sub(sub(A,mm(tr(B),X)),mm(tr(X),B)),mm(mm(tr(X),D),X))
    if not is_psd(K):
        raise ValueError('Retained positivity failed')
    if not is_psd(sub(sc(rho*rho,M),mm(tr(R),R))):
        raise ValueError('Relative residual inequality failed')
    eta=rho*rho/(delta-2*rho)
    return eta,K,M,R


def sharpness_counterexample(delta,rho,claimed_eta):
    """Exact counterexample to a stronger bound using only delta,rho,K>=0.

If delta<=2rho the same construction defeats any finite nonnegative eta.
This is a limit of the abstract information, not of a particular Hamiltonian.
"""
    delta=F(delta);rho=F(rho);claimed_eta=F(claimed_eta)
    if delta<=0 or rho<=0 or claimed_eta<0:
        raise ValueError('Positive delta,rho and nonnegative claimed eta required')
    a=delta-2*rho
    negative_coefficient=rho*rho-a*claimed_eta
    if negative_coefficient<=0:
        raise ValueError('Claimed bound is not stronger than the sharp bound')
    threshold=claimed_eta*(delta+claimed_eta)/negative_coefficient
    k=isqrt(threshold.numerator//threshold.denominator)+1
    X=[[F(k)]];D=[[delta]];B=[[(delta-rho)*k]];A=[[a*k*k]]
    determinant=(A[0][0]+claimed_eta)*(delta+claimed_eta)-B[0][0]**2
    if determinant>=0:
        raise AssertionError('Counterexample determinant did not become negative')
    K=A[0][0]-2*B[0][0]*k+delta*k*k
    residual=B[0][0]-delta*k
    if K!=0 or residual**2>rho**2*(1+k*k):
        raise AssertionError('Counterexample does not satisfy original hypotheses')
    return {'delta':str(delta),'rho':str(rho),'proposed_eta':str(claimed_eta),
            'X':k,'H_minus_e':[[str(A[0][0]),str(B[0][0])],
                                  [str(B[0][0]),str(delta)]],
            'determinant_at_proposed_lower':str(determinant),
            'original_abstract_hypotheses_satisfied':True}


def check_square(payload):
    if payload.get('kind')!='square_hubbard_metric_response_v1':
        raise ValueError('Wrong exact model/certificate kind')
    labels,H,P,Q=square_model()
    e=F(payload['reference_energy_over_t']);rho=F(payload['rho_over_t'])
    X=[[F(x) for x in row] for row in payload['response']]
    v=[[F(x)] for x in payload['retained_trial']]
    if len(X)!=len(Q) or any(len(row)!=len(P) for row in X) or len(v)!=len(P):
        raise ValueError('Response/trial dimensions do not match fixed sector')
    if e>=0:
        raise ValueError('Local gap proof needs a negative reference energy')
    # Each of four physical bonds has T_ij >= -2 I, independently verified on
    # all 16 local occupations. On Q, total doublon number >=1; hence QHQ>=0.
    passed,b=verify_bond_all_sectors(F(0),F(1))
    if not passed or b!=2 or 8-4*b!=0:
        raise AssertionError('Local eliminated-sector lower bound failed')
    A=sub(block(H,P,P),sc(e,eye(len(P))))
    B=block(H,Q,P);D=sub(block(H,Q,Q),sc(e,eye(len(Q))))
    delta=-e
    # The determinant-sized PSD test is deliberately NOT used for the gap;
    # its premise is the independently reconstructed local bond proof above.
    eta,K,M,R=energy_bound(A,B,D,X,delta,rho,verify_gap=False)
    norm=mm(mm(tr(v),M),v)[0][0]
    if norm<=0:
        raise ValueError('Nonzero trial state required')
    upper=e+mm(mm(tr(v),K),v)[0][0]/norm
    lower=e-eta
    absolute_rho=F(payload['absolute_rho_over_t'])
    if absolute_rho<0 or not is_psd(sub(sc(absolute_rho**2,eye(len(P))),mm(tr(R),R))):
        raise ValueError('Absolute residual comparison was not verified')
    absolute_eta=absolute_rho**2/delta
    if lower>upper:
        raise AssertionError('Inconsistent certified endpoints')
    return {'status':'accepted_exact_rational_finite_model',
        'energy_lower_over_t':str(lower),'energy_upper_over_t':str(upper),
        'width_over_t':str(upper-lower),'relative_allowance_over_t':str(eta),
        'absolute_allowance_over_t':str(absolute_eta),
        'global_determinants_enumerated':len(labels),'retained_dimension':len(P),
        'eliminated_dimension':len(Q),'local_gap_check_dimension':16,
        'eliminated_gap_over_t':str(delta),'retained_psd_check_dimension':len(P),
        'global_ground_state_teacher_used':False,
        'not_a_general_compact_solver':True}


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('certificate',type=Path)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():
        raise FileExistsError(a.out)
    result=check_square(json.loads(a.certificate.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
