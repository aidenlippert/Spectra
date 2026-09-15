"""Exact finite separators for occupation-polynomial square cones."""
from fractions import Fraction as F
from itertools import combinations
from math import lcm

from experiments.marginal_polynomial_metric import JointPolynomial, check_nonsingleton_separator


def integer_psd(matrix,initial_divisor=1):
    """Symmetric fraction-free elimination, including exact null pivots.

    A positive initial divisor can remove a known common factor from first
    Schur numerators. Every division is checked exactly. Each accepted update
    is a positive multiple of a Schur complement, including this first step.
    """
    n=len(matrix)
    if not 1<=n<=300 or any(len(row)!=n for row in matrix):
        raise ValueError('Bounded square moment matrix required')
    if any(type(x) is not int for row in matrix for x in row):
        raise ValueError('Integer moment matrix required')
    if any(matrix[i][j]!=matrix[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Exact symmetric moment matrix required')
    if type(initial_divisor) is not int or initial_divisor<=0:
        raise ValueError('Positive integer initial divisor required')
    a=[list(row) for row in matrix];previous=initial_divisor;rank=0
    for k in range(n):
        pivot=a[k][k]
        if pivot<0:raise ValueError('Negative moment pivot')
        if not pivot:
            if any(a[k][j] for j in range(k+1,n)):
                raise ValueError('Nonzero coupling from a null moment pivot')
            continue
        for i in range(k+1,n):
            for j in range(i,n):
                numerator=pivot*a[i][j]-a[i][k]*a[k][j]
                if numerator%previous:raise ValueError('Nonexact fraction-free moment elimination')
                a[i][j]=a[j][i]=numerator//previous
        previous=pivot;rank+=1
    return {'dimension':n,'rank':rank,'nullity':n-rank,'positive_semidefinite':True}


def polynomial_moments(functional,sites=None):
    """Caller first validates the signed physical-slice functional."""
    fractions=[F(x) for x in functional['weights']]
    denominator=lcm(*(x.denominator for x in fractions));moments={}
    for state,value in zip(functional['states'],fractions):
        weight=int(value*denominator);mask=state
        if sites is not None:weight*=sum(((state>>(2*i))&3)==3 for i in range(sites))-1
        while True:
            moments[mask]=moments.get(mask,0)+weight
            if not mask:break
            mask=(mask-1)&state
    return moments,denominator


def replay_separator(certificate):
    if certificate.get('kind')!='joint_polynomial_sos_separator_v1':
        raise ValueError('Unsupported polynomial-square separator')
    degree=certificate['square_degree']
    if type(degree) is not int or not 1<=degree<=3:
        raise ValueError('Square-polynomial degree must be one, two or three')
    charge_degree=certificate.get('charge_square_degree',0)
    if type(charge_degree) is not int or not 0<=charge_degree<=2:
        raise ValueError('Localized square degree must be zero (absent), one or two')
    if type(certificate.get('target_lower')) is not str:
        raise ValueError('Exact gap threshold required')
    ring=JointPolynomial(certificate)
    if ring.modes>12:raise ValueError('Finite square separator supports at most twelve modes')
    _,poly,cost=ring.compile(F(certificate['target_lower']))
    result=check_nonsingleton_separator(ring,poly,cost['numerator_scale'],certificate['functional'])
    moments,_=polynomial_moments(certificate['functional'])
    # Use the COMPLETE monomial basis. No numerical rank or selected span is trusted.
    basis=[sum(1<<i for i in indices) for k in range(degree+1)
           for indices in combinations(range(ring.modes),k)]
    matrix=[[moments.get(i|j,0) for j in basis] for i in basis]
    psd=integer_psd(matrix)
    result.update(square_degree=degree,moment_psd=psd,compilation=cost)
    if charge_degree:
        moments,_=polynomial_moments(certificate['functional'],ring.sites)
        basis=[sum(1<<i for i in indices) for k in range(charge_degree+1)
               for indices in combinations(range(ring.modes),k)]
        matrix=[[moments.get(i|j,0) for j in basis] for i in basis]
        result.update(charge_square_degree=charge_degree,charge_moment_psd=integer_psd(matrix))
    return dict(result,
        scope='Exact finite separator for the fixed metric and threshold: nonsingleton occupation indicators of degree<=6, (D-1) indicators of degree<=4, arbitrary unlocalized polynomial squares through square_degree, optional (D-1)-localized polynomial squares through charge_square_degree, and all fixed-spin number identities. Includes coefficient-L1 residual allowances. No restriction on the number of admitted squares. Does not exclude other localizers, higher-degree squares, other metrics or general representability.')


if __name__=='__main__':
    import argparse,json
    from pathlib import Path
    parser=argparse.ArgumentParser();parser.add_argument('--verify',required=True);args=parser.parse_args()
    print(json.dumps(replay_separator(json.loads(Path(args.verify).read_text())),indent=2))
