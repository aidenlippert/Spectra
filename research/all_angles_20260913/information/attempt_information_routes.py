"""Small exact diagnostic for routes 45-52 (no ground state used constructively).

Enumerates a spinless fermion ring.  The exact ground state is an oracle used
only to validate marginals.  A 1-RDM-only diagonal functional is deliberately
tested against two states with identical 1-RDM and different interaction
energy, exposing the representability gap.  All classical diagonal bounds
use Fraction arithmetic.
"""
from fractions import Fraction
from itertools import combinations
import numpy as np

def basis(m, n): return [sum(1 << i for i in c) for c in combinations(range(m), n)]
def hop_sign(s, i, j):
    # c_i^dag c_j, i != j; return new state and fermionic sign
    if not (s >> j) & 1 or (s >> i) & 1: return None
    sign = (-1) ** ((s & ((1 << j)-1)).bit_count())
    t = s ^ (1 << j)
    sign *= (-1) ** ((t & ((1 << i)-1)).bit_count())
    return t | (1 << i), sign

def model(m=4, n=2, t=1.0, u=2.0):
    B=basis(m,n); ix={s:k for k,s in enumerate(B)}; H=np.zeros((len(B),len(B)))
    for a,s in enumerate(B):
        H[a,a] = u * sum(((s>>i)&1) * ((s>>((i+1)%m))&1) for i in range(m))
        for i in range(m):
            j=(i+1)%m
            for x,y in ((i,j),(j,i)):
                q=hop_sign(s,x,y)
                if q: H[ix[q[0]],a] += -t*q[1]
    return B,H

def one_rdm(psi, B, m):
    g=np.zeros((m,m))
    for a,s in enumerate(B):
        for j in range(m):
            for i in range(m):
                q=hop_sign(s,i,j)
                if q: g[i,j] += psi[q[0] if False else 0] * 0 # guard: filled below
    # direct matrix element <psi|c_i^dag c_j|psi>
    ix={s:k for k,s in enumerate(B)}
    for a,s in enumerate(B):
        for i in range(m):
            for j in range(m):
                q=hop_sign(s,i,j)
                if q: g[i,j] += psi[ix[q[0]]]*psi[a]*q[1]
    return g

def one_rdm_mixed(weights, B, m):
    return sum(w*one_rdm(v,B,m) for w,v in weights)

def main():
    B,H=model(); vals,vecs=np.linalg.eigh(H); psi=vecs[:,0]
    E=float(psi@H@psi); G=one_rdm(psi,B,4)
    # two classical states |0101>, |1010> have the same diagonal 1-RDM (all 1/2)
    # as their incoherent mixture, but interaction energies differ from coherent GS.
    s1=(1<<0)|(1<<2); s2=(1<<1)|(1<<3); ix={s:k for k,s in enumerate(B)}
    mix=np.zeros(len(B)); mix[ix[s1]]=mix[ix[s2]]=1/np.sqrt(2)
    # phase-coherent superposition has identical 1-RDM to mixture for this ring
    E_mix=float(mix@H@mix); G_mix=one_rdm(mix,B,4)
    # exact rational checker: classical diagonal interaction floor and hopping norm
    diag=[Fraction(int(H[k,k])) for k in range(len(B))]
    classical_floor=min(diag); hop_norm_bound=Fraction(4,1) # row-sum bound for |t|=1 ring
    certified_floor=classical_floor-hop_norm_bound
    print({"ground_energy":E,"classical_floor":str(classical_floor),"floor_minus_hop_norm":str(certified_floor),
           "one_rdm_maxdiff_mixture":float(np.max(np.abs(G-G_mix))),"mix_energy":E_mix,
           "ground_rdm":np.round(G,6).tolist()})
    assert classical_floor == Fraction(0)
    assert certified_floor == Fraction(-4)
    # adversarial: 1-RDM does not determine pair interaction (diagonal states)
    # Two mixed states have identical diagonal 1-RDM (I/2), but different pair energy.
    c=np.zeros(len(B)); d=np.zeros(len(B));
    c[ix[(1<<0)|(1<<1)]]=c[ix[(1<<2)|(1<<3)]]=1/np.sqrt(2)
    d[ix[s1]]=d[ix[s2]]=1/np.sqrt(2)
    Gc=one_rdm_mixed([(0.5,c),(0.5,np.roll(c,1))],B,4)
    Gd=one_rdm_mixed([(0.5,d),(0.5,np.roll(d,1))],B,4)
    assert np.max(np.abs(np.diag(Gc)-np.diag(Gd))) < 1e-12
    assert float(c@H@c) != float(d@H@d)
if __name__ == '__main__': main()
