"""Exact six-dimensional Schur bounds for centered unequal pair hopping.

For the specified M=10,N=5 model each conserved pair-charge sector has
H0=15/4+Jz^2-(2/5)Jx. The maximal-spin j=5/2 subspace P occurs only when
all five pairs are singly occupied. Every spin block in Q=I-P has j=3/2
or j=1/2; the d=1,2 Jacobi matrices bound them, including multiplicities.

V=-sum_i (i-2) sigma_x(i) preserves the pair charges. Symmetric averaging
gives PVP=0, PV^2P=(25 I-4 Jx^2)/2, and ||V||<=6. If QH0Q>=c Q and
q=c-6|epsilon|-b>0, the Schur condition
    PH0P-bP-epsilon^2 PV^2P/q > 0
proves H0+epsilon V >= bI. This is a specialized spin-sector certificate;
it neither supplies a general SOS certificate nor locates a marginal cone.
"""
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import time

from experiments.marginal_sector_reference import bracket
from experiments.marginal_asymmetric_reference import model
from experiments.marginal_symbolic import decode, encode
from experiments.marginal_transfer_verify import apply_word

ROOT=Path(__file__).resolve().parents[1]


def matmul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def ldl_pivots(matrix):
    n=len(matrix)
    if not n or any(len(row)!=n for row in matrix):raise ValueError('Square matrix required')
    a=[[F(x) for x in row] for row in matrix]
    if any(a[i][j]!=a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Exact symmetric matrix required')
    for i in range(n):
        pivot=a[i][i]-sum(a[i][k]**2*a[k][k] for k in range(i))
        if pivot<=0:return None
        for j in range(i+1,n):
            a[j][i]=(a[j][i]-sum(a[j][k]*a[k][k]*a[i][k] for k in range(i)))/pivot
        a[i][i]=pivot
    return [a[i][i] for i in range(n)]


def collective_matrix():
    h=[[F(0) for _ in range(6)] for _ in range(6)]
    for k in range(6):
        h[k][k]=F(15,4)+(k-F(5,2))**2
        if k<5:
            h[k][k+1]=-F(5-k,5)
            h[k+1][k]=-F(k+1,5)
    return h


def square_B():
    x=[[F(0) for _ in range(6)] for _ in range(6)]
    for k in range(5):
        x[k][k+1]=F(5-k,2);x[k+1][k]=F(k+1,2)
    x2=matmul(x,x)
    return [[(F(25,2) if i==j else F(0))-2*x2[i][j]
             for j in range(6)] for i in range(6)]


def complement_lower():
    return min(F(bracket(10,F(1,5),d)['lower']) for d in (1,2))


def parity_matrix(matrix,sign=1):
    n=len(matrix);half=n//2
    return [[matrix[i][j]+sign*matrix[i][n-1-j] for j in range(half)] for i in range(half)]


def parity_sector_checks():
    # In five spins, physical F=product(sigma_x) acts on spin j as
    # (-1)^(5/2-j) times Dicke reflection: each singlet contributes -1.
    # Thus the F=+ complement consists of j=3/2 reflection- and j=1/2
    # reflection+. The latter has energy exactly 19/5.
    h3=[[F(0) for _ in range(4)] for _ in range(4)]
    for k in range(4):
        h3[k][k]=F(15,4)+(k-F(3,2))**2
        if k<3:h3[k][k+1]=-F(3-k,5);h3[k+1][k]=-F(k+1,5)
    odd3=parity_matrix(h3,-1)
    tested=[[2*comb(3,i)*(odd3[i][j]-(F(19,5) if i==j else 0)) for j in range(2)] for i in range(2)]
    if ldl_pivots(tested) is None:raise ValueError('Complement parity bound failed')
    # All F=- states and all other pair-charge sectors lie above 7/2.
    # For the maximal spin F=- block check directly; the lower spins and
    # other charges are bounded by the d=1,2 sector brackets.
    odd5=parity_matrix(collective_matrix(),-1)
    tested=[[2*comb(5,i)*(odd5[i][j]-(F(7,2) if i==j else 0)) for j in range(3)] for i in range(3)]
    if ldl_pivots(tested) is None or complement_lower()<=F(7,2):
        raise ValueError('Excluded-sector bound failed')
    return F(19,5),F(7,2)


def coupling_support_check():
    # Exact auxiliary five-spin check: Z contains symmetric functions f(k).
    # J^2=-5/4 I+sum_(i<j) Swap_ij. Hence sum Swap(VZ)=5 VZ proves
    # VZ lies entirely in j=3/2, since j(j+1)=15/4 there.
    image=[[sum(-(i-2) for i in range(5) if (s^(1<<i)).bit_count()==k)
            for k in range(6)] for s in range(32)]
    for s in range(32):
        for k in range(6):
            total=0
            for i in range(5):
                for j in range(i+1,5):
                    swapped=s^((1<<i)|(1<<j)) if ((s>>i)^(s>>j))&1 else s
                    total+=image[swapped][k]
            if total!=5*image[s][k]:raise ValueError('Spin selection rule failed')
    for k in range(6):
        for l in range(6):
            if sum(image[s][l] for s in range(32) if s.bit_count()==k):
                raise ValueError('Centered perturbation has a symmetric component')
    return True


def resolved_complement_lower():
    coupling_support_check()
    # j=3/2 reflection- block in its weighted two-dimensional basis.
    c=F(216411,50000)
    block=[[F(6),-F(3,5)],[-F(1,5),F(22,5)]]
    weighted=[[(2,6)[i]*(block[i][j]-(c if i==j else 0)) for j in range(2)] for i in range(2)]
    if ldl_pivots(weighted) is None:raise ValueError('Resolved spin-gap bound failed')
    return c


def schur_pivots(epsilon,b,qsector=None,parity=False,resolved=False):
    epsilon=F(epsilon);b=F(b)
    parity=parity or resolved
    if parity:
        qsector,excluded=parity_sector_checks()
        if excluded-6*abs(epsilon)<=b:return None
    else:qsector=complement_lower() if qsector is None else F(qsector)
    q=qsector-6*abs(epsilon)-b
    if q<=0:return None
    if resolved:
        q=resolved_complement_lower()-6*abs(epsilon)-b-36*epsilon**2/q
        if q<=0:return None
    h=collective_matrix();square=square_B()
    if parity:h=parity_matrix(h);square=parity_matrix(square)
    size=len(h);metric=[(2 if parity else 1)*comb(5,i) for i in range(size)]
    # f(k) is the common amplitude of each weight-k bitstring. Its metric
    # is diag(binomial(5,k)); this congruence must precede rational LDL.
    matrix=[[metric[i]*(h[i][j]-(b if i==j else 0)-epsilon**2*square[i][j]/q)
             for j in range(size)] for i in range(size)]
    return ldl_pivots(matrix)


def lower(epsilon=F(1,1000),tolerance=F(1,10**10),parity=False,resolved=False):
    epsilon=F(epsilon);tolerance=F(tolerance)
    parity=parity or resolved
    if abs(epsilon)>F(1,100):raise ValueError('This search supports |epsilon|<=1/100')
    if not 0<tolerance<F(1,10):raise ValueError('Invalid positive search tolerance')
    qsector=complement_lower();lo=F(3)-6*abs(epsilon);hi=F(33,10)
    if schur_pivots(epsilon,lo,qsector,parity,resolved) is None:raise ValueError('Initial lower bound failed')
    if schur_pivots(epsilon,hi,qsector,parity,resolved) is not None:raise ValueError('Initial upper search endpoint passed')
    while hi-lo>tolerance:
        mid=(lo+hi)/2
        if schur_pivots(epsilon,mid,qsector,parity,resolved) is not None:lo=mid
        else:hi=mid
    return {'lower':str(lo),'failed_search_endpoint':str(hi),'qsector':str(resolved_complement_lower() if resolved else (F(19,5) if parity else qsector)),'parity':parity,'resolved':resolved,
            'epsilon':str(epsilon),'tolerance':str(tolerance),'operator_norm_bound':'6'}


def replay(certificate):
    if certificate.get('kind') not in ('matched_spin_schur_v1','matched_spin_schur_parity_v1','matched_spin_schur_resolved_v1') or certificate.get('modes')!=10 or certificate.get('particles')!=5:
        raise ValueError('Unsupported spin-sector certificate')
    epsilon=F(certificate['epsilon']);h=decode(certificate['hamiltonian'],10,4)
    if h!=model(epsilon):raise ValueError('Hamiltonian does not match the spin-sector model')
    resolved=certificate['kind']=='matched_spin_schur_resolved_v1'
    parity=certificate['kind']=='matched_spin_schur_parity_v1' or resolved
    b=F(certificate['lower']);pivots=schur_pivots(epsilon,b,parity=parity,resolved=resolved)
    if pivots is None:raise ValueError('Schur lower bound failed exact positivity')
    states=[s for s in range(1<<10) if s.bit_count()==5];index={s:i for i,s in enumerate(states)}
    witness=certificate.get('independent_upper',{});amps=witness.get('amplitudes')
    if type(amps) is not list or len(amps)!=len(states) or any(type(a) is not int for a in amps) or not any(amps):
        raise ValueError('Expected nonzero integer fixed-sector upper witness')
    norm=sum(a*a for a in amps);energy=F(0)
    for word,c in h.items():
        for state,a in zip(states,amps):
            if not a:continue
            target=apply_word(word,state)
            if target:
                dest,sign=target;energy+=c*a*amps[index[dest]]*sign
    upper=energy/norm
    if upper<b:raise ValueError('Inconsistent interval')
    return {'lower':str(b),'lower_float':float(b),'upper':str(upper),'upper_float':float(upper),
            'width':str(upper-b),'width_float':float(upper-b),'schur_dimension':len(pivots),
            'qsector':str(resolved_complement_lower() if resolved else (F(19,5) if parity else complement_lower())),'schur_pivots':[str(d) for d in pivots],
            'hamiltonian_bound':True,'scope':'Specified half-filled pair-charge-conserving spin family'}


def run(epsilon,parity=False,resolved=False):
    parity=parity or resolved
    epsilon=F(epsilon);out=ROOT/'results/marginal_schur_transfer'/(str(epsilon).replace('/','_')+('_resolved' if resolved else ('_parity' if parity else '')))
    if out.exists():raise ValueError('Preserve previous run')
    started=time.monotonic();search=lower(epsilon,parity=parity,resolved=resolved)
    from experiments.marginal_symmetry_transfer import integer_upper
    h=model(epsilon)
    certificate={'kind':'matched_spin_schur_resolved_v1' if resolved else ('matched_spin_schur_parity_v1' if parity else 'matched_spin_schur_v1'),'modes':10,'particles':5,'epsilon':str(epsilon),
                 'lower':search['lower'],'hamiltonian':encode(h),'independent_upper':integer_upper(h,10,5)}
    receipt=replay(certificate);receipt['total_seconds']=time.monotonic()-started
    out.mkdir(parents=True)
    for name,obj in [('certificate',certificate),('receipt',receipt),('search',search)]:
        (out/(name+'.json')).write_text(json.dumps(obj,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('lower_float','upper_float','width_float','total_seconds')}),flush=True)
    return receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--epsilon',default='1/1000');parser.add_argument('--verify',type=Path)
    parser.add_argument('--parity',action='store_true')
    parser.add_argument('--resolved',action='store_true')
    args=parser.parse_args()
    if args.verify:print(json.dumps(replay(json.loads(args.verify.read_text())),indent=2))
    else:run(F(args.epsilon),args.parity,args.resolved)
