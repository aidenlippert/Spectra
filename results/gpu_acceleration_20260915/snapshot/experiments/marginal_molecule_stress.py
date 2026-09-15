"""Stretched H4 square STO-3G cubic marginal certificate."""
from fractions import Fraction as F
from pathlib import Path
import json, math, numpy as np


def parity_generators(h, modes):
    """Independent binary phase symmetries from a GF(2) support nullspace."""
    from experiments.marginal_symbolic import canonical, validate_word
    if type(modes) is not int or modes<1:raise ValueError('Invalid mode count')
    for w in h:validate_word(w,modes,4)
    pivots={}
    for word in canonical(h):
        row=0
        for c,i in word:row^=1<<i
        while row:
            pivot=(row&-row).bit_length()-1
            if pivot in pivots:row^=pivots[pivot]
            else:pivots[pivot]=row;break
    generators=[]
    for free in range(modes):
        if free in pivots:continue
        vector=1<<free
        for pivot in sorted(pivots,reverse=True):
            if (pivots[pivot]&vector).bit_count()%2:vector^=1<<pivot
        generators.append(vector)
    return generators

def build_h():
    from pyscf import gto, scf, fci, ao2mo
    geom=[['H',(0.,0.,0.)],['H',(2.,0.,0.)],['H',(0.,2.,0.)],['H',(2.,2.,0.)]]
    mol=gto.M(atom=geom,basis='sto-3g',unit='Angstrom',spin=0,charge=0,verbose=0)
    mf=scf.RHF(mol).run(conv_tol=1e-12,verbose=0)
    hcore=mf.mo_coeff.T@mf.get_hcore()@mf.mo_coeff
    eri=ao2mo.restore(1,ao2mo.kernel(mol,mf.mo_coeff),mol.nao_nr()); n=mol.nao_nr(); ns=2*n
    hs=np.zeros((ns,ns)); hs[0::2,0::2]=hcore; hs[1::2,1::2]=hcore
    ep=eri.reshape(n,n,n,n); se=np.zeros((ns,ns,ns,ns))
    for p in range(ns):
      for q in range(ns):
       for r in range(ns):
        for s in range(ns):
         se[p,q,r,s]=ep[p//2,r//2,q//2,s//2]*(p%2==r%2)*(q%2==s%2)-ep[p//2,s//2,q//2,r//2]*(p%2==s%2)*(q%2==r%2)
    raw={((1,p),(0,q)):float(hs[p,q]) for p in range(ns) for q in range(ns)}
    raw.update({((1,p),(1,q),(0,s),(0,r)):float(se[p,q,r,s]/4) for p in range(ns) for q in range(ns) for r in range(ns) for s in range(ns)})
    from experiments.marginal_symbolic import canonical
    hraw=canonical({w:F(float(v)) for w,v in raw.items()})
    h=canonical({w:F(float(v)).limit_denominator(10**10) for w,v in raw.items() if abs(v)>1e-13})
    ah=canonical({tuple((1-c,i) for c,i in reversed(w)):v for w,v in h.items()})
    h=canonical({w:(h.get(w,0)+ah.get(w,0))/2 for w in set(h)|set(ah)})
    # independent FCI and exact CAR energy are recorded separately
    ef,ci=fci.direct_spin1.kernel(hcore,eri,n,4,ecore=mol.energy_nuc(),conv_tol=1e-12)
    return mol,geom,h, float(ef-mol.energy_nuc()), float(mol.energy_nuc()), hraw

def car_ed(h,modes=8,particles=4):
 from itertools import combinations
 states=[sum(1<<i for i in c) for c in combinations(range(modes),particles)]; pos={s:i for i,s in enumerate(states)}; A=np.zeros((len(states),len(states)))
 for w,c in h.items():
  for j,s0 in enumerate(states):
   s=s0; sign=1
   for cr,i in reversed(w):
    sign*= -1 if ((s&((1<<i)-1)).bit_count()&1) else 1
    if cr:
     if s>>i&1: break
     s|=1<<i
    else:
     if not s>>i&1: break
     s&=~(1<<i)
   else: A[pos[s],j]+=float(c)*sign
 return float(np.linalg.eigvalsh(A)[0])

def solve_certificate(h):
 from itertools import combinations
 from experiments.marginal_adaptive import assemble_and_solve
 from experiments.marginal_coefficient import export
 masks=parity_generators(h,8)
 blocks={}
 for d in range(4):
  for nc in range(d+1):
   for left in combinations(range(8),nc):
    for right in combinations(range(8),d-nc):
     w=tuple((1,i) for i in left)+tuple((0,i) for i in right)
     # spin conservation: alpha modes are even, beta odd
     q=(sum((2*c-1) for c,i in w if i%2==0),sum((2*c-1) for c,i in w if i%2==1))
     q+=tuple(sum((s>>i)&1 for c,i in w)%2 for s in masks)
     blocks.setdefault(q,[]).append(w)
 blocks=[{'name':str(q),'words':v} for q,v in sorted(blocks.items())]
 print(json.dumps({'parity_masks':masks,'blocks':len(blocks),'largest_block':max(len(b['words']) for b in blocks)}),flush=True)
 proposal,_=assemble_and_solve(h,8,4,blocks,degree=3,groups=[[0,2,4,6],[1,3,5,7]],residual_penalty=False,parity_masks=masks)
 cert,receipt=export(h,8,4,blocks,proposal,denominator=10**9)
 receipt['parity_masks']=masks
 return cert,receipt

def write_fixture():
 mol,geom,h,ef,en,hraw=build_h(); out=Path('results/marginal_molecule_stress'); out.mkdir(parents=True,exist_ok=True)
 from experiments.marginal_symbolic import encode, canonical, adj
 hermitian_raw=canonical(adj(hraw))
 symmetric={w:(hraw.get(w,0)+hermitian_raw.get(w,0))/2 for w in set(hraw)|set(hermitian_raw)}
 den=10**12; rounded={w:F(round(v*den),den) for w,v in symmetric.items() if round(v*den)}
 err=sum(abs(rounded.get(w,0)-hraw.get(w,0)) for w in set(hraw)|set(rounded))
 from experiments.marginal_transfer import upper_ed
 upper=upper_ed(rounded,8,4)
 fixture={'modes':8,'particles':4,'geometry':geom,'basis':'STO-3G','orbital_basis':'RHF canonical MO','nuclear_repulsion':en,'electronic_fci_total':ef,'car_rational_ground':car_ed(rounded),'car_vs_fci_error':abs(car_ed(rounded)-ef),'total_fci':ef+en,'integral_rounding_l1_error':str(err),'numerical_integral_perturbation_l1_upper':str(err),'coefficient_denominator':den,'hamiltonian':[{'word':[[c,i] for c,i in w],'coefficient':str(v)} for w,v in sorted(rounded.items())]}
 (out/'h4_square_sto3g.json').write_text(json.dumps(fixture,indent=2)+'\n')
 cert={'modes':8,'particles':4,'hamiltonian':fixture['hamiltonian'],'independent_upper':upper}; receipt={'status':'upper_only_fixture','fci_electronic':ef,'car_ground':fixture['car_rational_ground']}
 (out/'upper_only_fixture.json').write_text(json.dumps(cert)+'\n'); (out/'upper_only_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({k:fixture[k] for k in ('electronic_fci_total','car_rational_ground','car_vs_fci_error','total_fci')},indent=2)); print(json.dumps({'lower':receipt.get('lower'),'upper':receipt.get('upper'),'status':receipt.get('status')},indent=2))


def run_saved():
 from experiments.marginal_symbolic import decode
 from experiments.marginal_transfer import upper_ed
 from experiments.marginal_transfer_verify import replay
 out=Path(__file__).resolve().parents[1]/'results/marginal_molecule_stress'
 fixture=json.loads((out/'h4_square_sto3g.json').read_text())
 h=decode(fixture['hamiltonian'],8,4);cert,receipt=solve_certificate(h)
 cert['independent_upper']=upper_ed(h,8,4);receipt.update(replay(cert))
 budget=F(fixture['numerical_integral_perturbation_l1_upper'])
 receipt['numerical_integral_perturbation_l1_upper']=str(budget)
 receipt['coefficient_perturbation_enclosed_lower']=str(F(receipt['lower'])-budget)
 receipt['coefficient_perturbation_enclosed_upper']=str(F(receipt['upper'])+budget)
 (out/'accepted_degree3_certificate.json').write_text(json.dumps(cert)+'\n')
 (out/'accepted_degree3_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({k:v for k,v in receipt.items() if not isinstance(v,list)}),flush=True)
 return cert,receipt


if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--solve-certificate',action='store_true')
 args=parser.parse_args()
 if args.solve_certificate:run_saved()
 else:write_fixture()
