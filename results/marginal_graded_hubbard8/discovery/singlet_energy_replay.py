"""Bounded H8 singlet Schur research certificate. Exact replay CLI.

Ground-spin identification uses Lieb's repulsive half-filled bipartite Hubbard
 theorem (PRL 62,1201); the model hypotheses are checked, the theorem is not
re-proved here. The independent embedding verifier establishes complete valence
S=0 coverage; charge DP bounds its ionic complement. No existing gate is raised.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from results.marginal_graded_hubbard8.singlet_valence_embedding.independent_replay import verify as verify_embedding, require, gram, ldl
from experiments.marginal_joint_charge_dp import replay as verify_gap, model
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_spin_reduction import spin_operators, commutator
from experiments.marginal_symbolic import decode
from experiments.marginal_enlarged_schur import solve_positive

THEOREM='Lieb-1989-PRL-62-1201-repulsive-half-filled-bipartite'
def apply(oracle,v,excluded=()):
    out={}
    for s,a in v.items():
        for t,x in oracle.action(s).items():
            if t not in excluded:out[t]=out.get(t,F(0))+a*x
    return {s:a for s,a in out.items() if a}

def prepare(c):
    require(c.get('kind')=='h8_singlet_polynomial_response_research_v1','research certificate kind')
    require(c.get('ground_spin_theorem')==THEOREM,'explicit ground-spin theorem premise')
    e=c['embedding']; embedding=verify_embedding(e)
    gap=c['complement']; n,U,t,gamma,*_=model(gap)
    require(n==8 and U==4 and t==1,'strict connected repulsive balanced bipartite model')
    h=json.loads((Path(__file__).resolve().parents[1]/'hamiltonian.json').read_text())
    require(decode(gap['hamiltonian'],16,4)==decode(h['hamiltonian'],16,4),'same original Hamiltonian')
    gap_receipt=verify_gap(gap)
    oracle=DeterminantOracle(h); plus,minus,z=spin_operators(16)
    require(not commutator(oracle.h,plus) and not commutator(oracle.h,z),'exact spin symmetry')
    P=set(e['valence_states']); V=[{int(s):F(a) for s,a in v.items()} for v in e['basis']]
    W=[apply(oracle,v) for v in V]
    require(all(not (set(w)&P) for w in W),'valence H block zero')
    coefficients=c.get('response_polynomial')
    require(type(coefficients) is list and len(coefficients)<=8,'response polynomial order budget')
    require(all(type(x) is str and len(x)<=128 for x in coefficients),'bounded rational coefficients')
    coefficients=list(map(F,coefficients))
    require(all(max(abs(x.numerator),x.denominator)<=10**30 for x in coefficients),'coefficient magnitude budget')
    require(not coefficients or coefficients[-1]!=0,'canonical polynomial degree')
    B=[]
    if coefficients:
        for w in W:
            power=w; b={}
            for k,a in enumerate(coefficients):
                for s,x in power.items(): b[s]=b.get(s,0)+a*x
                if k+1<len(coefficients):power=apply(oracle,power,P)
            B.append({s:x for s,x in b.items() if x})
        require(all(x>0 for x in ldl(gram(B,B))),'independent response basis')
    A=[apply(oracle,b,P) for b in B]
    data={'G':gram(V,V),'leakage':gram(W,W),'metric':gram(B,B),'hamiltonian':gram(B,A),
          'squared':gram(A,A),'coupling':gram(B,W),'h_coupling':gram(A,W),'gamma':gamma}
    return data,oracle,V,W,B,{'embedding':embedding,'complement':gap_receipt,
        'complement_transfer':'Within S=0, the verified complete valence kernel is P. Its orthogonal complement contains only ionic D>=1 states, so the original ionic charge-DP bound applies.',
        'ground_spin_hypotheses':'Exact original U=4, t=1 uniform connected open eight-site chain, eight particles, real spin-independent hopping, bipartition sizes 4 and 4.',
        'ground_spin_source':'https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.1201'}

def schur(data,tau):
    gamma=data['gamma']; require(tau<gamma,'positive Schur denominator')
    count=len(data['G']); r=len(data['metric'])
    corr=[[F(0) for _ in range(count)] for _ in range(count)]
    if r:
        D=[[data['squared'][i][j]-(tau+gamma)*data['hamiltonian'][i][j]+tau*gamma*data['metric'][i][j]
            for j in range(r)] for i in range(r)]
        L=[[data['h_coupling'][i][j]-gamma*data['coupling'][i][j] for j in range(count)] for i in range(r)]
        solved=solve_positive(D,L)
        corr=[[sum(L[k][i]*solved[k][j] for k in range(r)) for j in range(count)] for i in range(count)]
    return [[-tau*data['G'][i][j]-(data['leakage'][i][j]-corr[i][j])/(gamma-tau)
             for j in range(count)] for i in range(count)]

def replay(c):
    data,oracle,V,W,B,checks=prepare(c)
    require(type(c['lower']) is str and len(c['lower'])<=128,'bounded lower endpoint')
    lower=F(c['lower']); pivots=ldl(schur(data,lower))
    require(all(x>0 for x in pivots),'positive generalized response Schur')
    upper=oracle.upper(c['independent_upper'])
    require(lower<=upper,'consistent interval')
    return {'accepted':True,'lower':str(lower),'upper':str(upper),'width':str(upper-lower),
            'lower_float':float(lower),'upper_float':float(upper),'width_float':float(upper-lower),
            'retained_dimension':14,'response_dimension':len(B),'polynomial_degree':len(c['response_polynomial'])-1,
            'response_support_sizes':list(map(len,B)),'unique_action_states':len(oracle.cache),
            'referenced_determinants':oracle.referenced_state_count(),
            'schur_pivots':list(map(str,pivots)), 'checks':checks,
            'ground_spin_theorem':THEOREM,
            'scope':'Exact interval for the fixed eight-site half-filled open U4 Hubbard chain, using the stated Lieb theorem. Finite sparse research verifier; no universal chemistry or scaling claim.'}
if __name__=='__main__':
    print(json.dumps(replay(json.loads(Path(sys.argv[1]).read_text())),indent=2))
