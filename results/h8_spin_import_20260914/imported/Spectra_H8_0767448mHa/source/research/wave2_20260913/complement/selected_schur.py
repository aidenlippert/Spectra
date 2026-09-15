"""Exact rational selected-subspace Schur controls for saved H4/H6 fixtures."""
from fractions import Fraction as F
from pathlib import Path
import json, math, hashlib
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.marginal_transfer_verify import apply_word

ROOT = Path(__file__).resolve().parents[3]

def matrix(fixture):
    d=json.loads((ROOT/fixture).read_text()); modes=d['modes']; n=d['particles']
    basis=[s for s in range(1<<modes) if s.bit_count()==n]; ix={s:i for i,s in enumerate(basis)}
    H=[[F(0)]*len(basis) for _ in basis]
    for term in d['hamiltonian']:
        w=tuple(tuple(x) for x in term['word']); c=F(term['coefficient'])
        for s in basis:
            z=apply_word(w,s)
            if z: t,sg=z; H[ix[t]][ix[s]] += c*sg
    return d,basis,H

def glb(M):
    return min(M[i][i]-sum(abs(M[i][j]) for j in range(len(M)) if j!=i) for i in range(len(M)))

def inertia_psd(M, l):
    n=len(M); a=[[F(M[i][j])-(l if i==j else 0) for j in range(n)] for i in range(n)]
    for k in range(n):
        if a[k][k] < 0: return False
        if a[k][k] == 0:
            if any(a[k][j] for j in range(k+1,n)): return False
            continue
        p=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n): a[i][j]-=a[i][k]*a[k][j]/p
    return True

def ldl_lower(M):
    import numpy as np
    ev=float(np.linalg.eigvalsh(np.array(M,dtype=float))[0]); scale=10**10
    l=F(math.floor(ev*scale),scale)
    while not inertia_psd(M,l): l -= F(1,scale)
    return l

def run(fixture, states):
    d,basis,H=matrix(fixture); ix={s:i for i,s in enumerate(basis)}; P=[ix[s] for s in states]
    Q=[i for i in range(len(H)) if i not in P]; A=[[H[i][j] for j in P] for i in P]
    C=[[H[i][j] for j in Q] for i in Q]; B=[[H[i][j] for j in Q] for i in P]
    a=ldl_lower(A); mu=ldl_lower(C); coupling=sum(x*x for row in B for x in row)
    gap=a-mu; radicand=gap*gap+4*coupling; rad_up=F(math.ceil(math.sqrt(float(radicand))*10**10),10**10)
    while rad_up*rad_up < radicand: rad_up += F(1,10**10)
    endpoint=(a+mu-rad_up)/2
    return {'fixture':fixture,'dimension':len(H),'selected_states':states,'p_dim':len(P),'q_dim':len(Q),
      'a_ldl':str(a),'mu_ldl':str(mu),'coupling_frobenius_sq':str(coupling),
      'schur_endpoint':str(endpoint),'status':'accepted_general_block_bound','mu_gt_a':mu>a,
      'full_complement_entries':len(Q)**2,'diagonal_only_mu_diagnostic':str(min(C[i][i] for i in range(len(C))))}

if __name__=='__main__':
    def states(p): return json.loads((ROOT/p).read_text())['independent_upper']['states']
    fixture='results/certificate_scaling/active_space_ladder/h4/fixture.json'
    print(json.dumps({'source_sha256':hashlib.sha256((ROOT/fixture).read_bytes()).hexdigest(),
      'h4_step01_16':run(fixture,states('results/all_angles_20260913/selected_refinement/h4/step_01_upper.json')),
      'h4_final20':run(fixture,states('results/all_angles_20260913/selected_refinement/h4/step_02_upper.json'))},indent=2))
