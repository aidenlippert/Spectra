"""H8 ground energy replay from exact signed-orbit full moments.

Numerical or submitted moment tables are never trusted. The original model,
complete singlet embedding and charge DP are checked; every full moment is
recomputed by the symmetry oracle. Nonorthogonal renewal and Schur arithmetic
then certify the lower bound, and a polynomial physical state gives the upper.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from results.marginal_graded_hubbard8.singlet_valence_embedding.independent_replay import verify as verify_embedding,require,ldl
from results.marginal_graded_hubbard8.discovery.singlet_energy_replay import THEOREM,schur
from experiments.marginal_joint_charge_dp import model,replay as verify_gap
from experiments.marginal_symbolic import decode
from experiments.marginal_spin_reduction import spin_operators,commutator
from experiments.marginal_symmetry_moments import SymmetryMomentOracle,projected_moments,simple
ROOT=Path(__file__).resolve().parents[1]

def response_data(M,coefficients,gamma):
    return block_response_data(M,[coefficients],gamma)

def block_response_data(M,polynomials,gamma):
    K=projected_moments(M);n=len(M[0]);polynomials=[list(map(F,p)) for p in polynomials]
    directions=[(a,i) for a in polynomials for i in range(n)]
    def linear(offset):return [[simple(sum(a[k]*K[k+offset][i][j] for k in range(len(a)))) for j in range(n)] for a,i in directions]
    def quadratic(offset):return [[simple(sum(a[k]*b[l]*K[k+l+offset][i][j] for k in range(len(a)) for l in range(len(b)))) for b,j in directions] for a,i in directions]
    return {'G':M[0],'leakage':K[0],'metric':quadratic(0),'hamiltonian':quadratic(1),
            'squared':quadratic(2),'coupling':linear(0),'h_coupling':linear(1),'gamma':gamma}

def polynomial_upper(M,coefficients):
    def contraction(shift):
        return sum(a*b*M[k+l+shift][i][j] for k,row in enumerate(coefficients)
                   for l,col in enumerate(coefficients) for i,a in enumerate(row) for j,b in enumerate(col))
    norm=contraction(0);require(norm>0,'positive physical polynomial witness norm')
    numerator=contraction(1)
    return F(numerator,norm),norm

def upper_monomials(coefficients):
    """Fixed Chebyshev basis T_k((H-8I)/12), transformed exactly."""
    polynomials=[[F(1)],[F(-2,3),F(1,12)]]
    for k in range(2,len(coefficients)):
        previous=polynomials[-1];before=polynomials[-2];p=[F(0)]*(k+1)
        for j,a in enumerate(previous):p[j]-=F(4,3)*a;p[j+1]+=F(1,6)*a
        for j,a in enumerate(before):p[j]-=a
        polynomials.append(p)
    return [[sum(coefficients[k][j]*polynomials[k][i] for k in range(i,len(coefficients)))
             for j in range(len(coefficients[0]))] for i in range(len(coefficients))]

def replay(c):
    started=time.monotonic()
    require(c.get('kind')=='h8_singlet_moment_energy_research_v1','moment energy certificate kind')
    require(c.get('ground_spin_theorem')==THEOREM,'explicit Lieb ground-spin theorem')
    embedding=verify_embedding(c['embedding']);gap=c['complement']
    n,U,t,gamma,*_=model(gap);require(n==8 and U==4 and t==1,'fixed repulsive connected bipartite H8 model')
    h=json.loads((ROOT/'hamiltonian.json').read_text())
    require(decode(gap['hamiltonian'],16,4)==decode(h['hamiltonian'],16,4),'same original Hamiltonian')
    gap_receipt=verify_gap(gap)
    require(('response_polynomial' in c)!=('response_polynomials' in c),'choose one response polynomial representation')
    raw=c.get('response_polynomials',[c.get('response_polynomial')])
    require(type(raw) is list and 1<=len(raw)<=2,'one or two response polynomial blocks')
    polynomials=[]
    for coefficients in raw:
        require(type(coefficients) is list and 1<=len(coefficients)<=11,'one to eleven response coefficients within moment-order24 budget')
        require(all(type(a) is str and len(a)<=128 for a in coefficients),'bounded rational response coefficients')
        a=list(map(F,coefficients));require(a[-1]!=0 and all(max(abs(x.numerator),x.denominator)<=10**30 for x in a),'response coefficient magnitude or degree')
        polynomials.append(a)
    require(('upper_polynomial_coefficients' in c)!=('upper_chebyshev_coefficients' in c),'choose one upper polynomial representation')
    upper_coeff=c.get('upper_chebyshev_coefficients',c.get('upper_polynomial_coefficients'))
    require(type(upper_coeff) is list and 1<=len(upper_coeff)<=12 and all(type(row) is list and len(row)==14 for row in upper_coeff),'bounded polynomial upper shape within moment-order24 budget')
    require(all(type(x) is int and abs(x)<=10**15 for row in upper_coeff for x in row),'bounded integer upper coefficients')
    require(any(x for row in upper_coeff for x in row),'nonzero upper recipe')
    upper_basis='chebyshev_shift8_scale12' if 'upper_chebyshev_coefficients' in c else 'monomial'
    if 'upper_chebyshev_coefficients' in c:upper_coeff=upper_monomials(upper_coeff)
    require(type(c['lower']) is str and len(c['lower'])<=128,'bounded lower endpoint');lower=F(c['lower'])
    require(lower<gamma,'positive Schur denominator')
    oracle=SymmetryMomentOracle(h);plus,_,z=spin_operators(16)
    require(not commutator(oracle.base.h,plus) and not commutator(oracle.base.h,z),'exact spin symmetry')
    V=[{int(s):x for s,x in v.items()} for v in c['embedding']['basis']]
    maximum_length=max(map(len,polynomials));order=max(2*maximum_length+2,2*len(upper_coeff)-1)
    M,work=oracle.moments(V,order)
    data=block_response_data(M,polynomials,gamma)
    pivots=ldl(schur(data,lower));require(all(x>0 for x in pivots),'positive moment response Schur')
    upper,norm=polynomial_upper(M,upper_coeff);require(lower<=upper,'consistent moment energy interval')
    return {'accepted':True,'lower':str(lower),'upper':str(upper),'width':str(upper-lower),
            'lower_float':float(lower),'upper_float':float(upper),'width_float':float(upper-lower),
            'response_degree':maximum_length-1,'response_block_count':len(polynomials),'upper_degree':len(upper_coeff)-1,'upper_basis':upper_basis,'retained_dimension':14,'response_dimension':14*len(polynomials),
            'polynomial_upper_norm':str(norm),'schur_pivots':list(map(str,pivots)),
            'moment_work':work,'embedding':embedding,'complement':gap_receipt,
            'ground_spin_theorem':THEOREM,
            'complement_transfer':'Complete valence singlet P leaves only ionic D>=1 states in S=0 Q.',
            'replay_seconds':time.monotonic()-started,
            'scope':'Exact original H8 ground interval with explicit Lieb theorem dependency. Full moments are recomputed in a signed finite-group quotient. No submitted moments, full determinant vectors or full sector matrix are used for moment/energy contraction. Finite orbit compression does not prove scalable representability.'}
if __name__=='__main__':print(json.dumps(replay(json.loads(Path(sys.argv[1]).read_text())),indent=2))
