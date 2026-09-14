"""H4 rectangular STO-3G fixture and spin-orbital CAR export.

The electronic Hamiltonian is exported in the convention
H = h_pq a†_p a_q + 1/4 (pq||rs) a†_p a†_q a_s a_r.
Numerical integrals are rounded to rational coefficients; the reported bound
is only as rigorous as that coefficient-rounding budget.
"""
from fractions import Fraction as F
from pathlib import Path
import json, numpy as np

def car_ed(hamiltonian, modes=8, particles=4):
    from itertools import combinations
    states=[sum(1<<i for i in c) for c in combinations(range(modes),particles)]; pos={s:i for i,s in enumerate(states)}
    A=np.zeros((len(states),len(states)))
    for w,c in hamiltonian.items():
      for j,s0 in enumerate(states):
       s=s0; sign=1
       for creation,i in reversed(w):
        sign *= -1 if ((s & ((1<<i)-1)).bit_count()&1) else 1
        if creation:
         if s>>i&1: break
         s|=1<<i
        else:
         if not s>>i&1: break
         s&=~(1<<i)
       else: A[pos[s],j]+=float(c)*sign
    return float(np.linalg.eigvalsh(A)[0])

def build():
    from pyscf import gto, scf, fci, ao2mo
    geom=[['H',(0.,0.,0.)],['H',(1.,0.,0.)],['H',(0.,1.5,0.)],['H',(1.,1.5,0.)]]
    mol=gto.M(atom=geom,basis='sto-3g',unit='Angstrom',spin=0,charge=0,verbose=0)
    mf=scf.RHF(mol).run(conv_tol=1e-12,verbose=0)
    hcore=mf.mo_coeff.T@mf.get_hcore()@mf.mo_coeff
    eri=ao2mo.restore(1,ao2mo.kernel(mol,mf.mo_coeff),mol.nao_nr())
    n=mol.nao_nr(); ns=2*n
    hs=np.zeros((ns,ns)); hs[0::2,0::2]=hcore; hs[1::2,1::2]=hcore
    # physicist ERI (pq|rs); antisymmetrize in the last pair.
    ep=eri.reshape(n,n,n,n)
    spin_eri=np.zeros((ns,ns,ns,ns))
    for p in range(ns):
      for q in range(ns):
       for r in range(ns):
        for s in range(ns):
         # PySCF's restored tensor is (p r | q s) in chemist notation;
         # antisymmetrize the ket pair for the a†p a†q a_s a_r convention.
         spin_eri[p,q,r,s]=(ep[p//2,r//2,q//2,s//2]*(p%2==r%2)*(q%2==s%2)
           -ep[p//2,s//2,q//2,r//2]*(p%2==s%2)*(q%2==r%2))
    def rat(x): return F(float(x)).limit_denominator(10**10)
    h={}; raw_float={}
    for p in range(ns):
      for q in range(ns):
       raw_float[((1,p),(0,q))]=float(hs[p,q])
       if abs(hs[p,q])>1e-13: h[((1,p),(0,q))]=rat(hs[p,q])
    for p in range(ns):
     for q in range(ns):
      for r in range(ns):
       for s in range(ns):
        c=spin_eri[p,q,r,s]/4
        raw_float[((1,p),(1,q),(0,s),(0,r))]=float(c)
        if abs(c)>1e-13: h[((1,p),(1,q),(0,s),(0,r))]=rat(c)
    # Independent PySCF FCI benchmark (electronic plus nuclear separately).
    e,ci=fci.direct_spin1.kernel(hcore,eri,n,4,ecore=mol.energy_nuc(),conv_tol=1e-12)
    # Canonical CAR collection is explicitly Hermitian; this also removes
    # duplicate tensor presentations before rational export.
    from experiments.marginal_symbolic import canonical, adj
    raw=canonical({w:F.from_float(v) for w,v in raw_float.items()})
    h=canonical({w:v for w,v in h.items()})
    ah=adj(h)
    h=canonical({w:(h.get(w,0)+ah.get(w,0))/2 for w in set(h)|set(ah)})
    coeff_error=sum(abs(h.get(w,0)-raw.get(w,0)) for w in set(h)|set(raw))
    car_ground=car_ed(h,ns,4)
    return {'modes':ns,'particles':4,'geometry':geom,'basis':'STO-3G','orbital_basis':'RHF canonical MO','nuclear_repulsion':mol.energy_nuc(),'electronic_fci_total':float(e-mol.energy_nuc()),'car_rational_ground':car_ground,'car_vs_fci_error':abs(car_ground-float(e-mol.energy_nuc())),'total_fci':float(e),'integral_rounding_l1_bound':float(coeff_error),'hamiltonian':[{'word':[[c,i] for c,i in w],'coefficient':str(v)} for w,v in sorted(h.items())]}

if __name__=='__main__':
 out=build(); p=Path('results/marginal_molecule'); p.mkdir(parents=True,exist_ok=True); (p/'h4_rectangle_sto3g.json').write_text(json.dumps(out,indent=2)); print(json.dumps({k:out[k] for k in ('modes','particles','nuclear_repulsion','electronic_fci_total','total_fci','integral_rounding_l1_bound')},indent=2))
