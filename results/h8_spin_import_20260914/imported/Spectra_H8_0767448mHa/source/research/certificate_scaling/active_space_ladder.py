"""Generate a straight H4/H6/H8/H10 STO-3G chain fixture ladder.

This is deliberately separate from the earlier square/rectangular H4 fixtures.
It exports the RHF-canonical spin-orbital Hamiltonian with a rational coefficient
enclosure and a cheap Hartree--Fock determinant upper witness.  FCI is optional
and only enabled for H4/H6 validation.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse, json, time, sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from experiments.marginal_symbolic import canonical, adj, encode, hermitian
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_implicit_certificate import rational_text

def build(natoms, output, denominator=10**12, fci_control=False, localized=False):
    import numpy as np
    from pyscf import gto, scf, ao2mo
    started = time.monotonic()
    geometry = [['H', [0., 0., float(F(7*i, 5))]] for i in range(natoms)]
    mol = gto.M(atom=geometry, basis='sto-3g', unit='Angstrom', spin=0, charge=0, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-12, verbose=0)
    if not mf.converged: raise RuntimeError('RHF did not converge')
    orbitals = mol.nao_nr(); modes = 2*orbitals; particles = mol.nelectron
    coeff = mf.mo_coeff
    localization = {'status':'not_run','reason':'canonical ladder baseline'}
    if localized:
      nocc = particles // 2
      # In this collinear chain, diagonalizing the projected z-coordinate is
      # the cheap one-axis Boys first-moment transform.  Apply it separately
      # to occupied and virtual spaces, preserving the HF determinant.
      overlap = mol.intor_symmetric('int1e_ovlp')
      xop = mol.intor_symmetric('int1e_r', comp=3)[2]
      from scipy.linalg import eigh
      uo = eigh(coeff[:, :nocc].T @ xop @ coeff[:, :nocc],
                coeff[:, :nocc].T @ overlap @ coeff[:, :nocc])[1]
      uv = eigh(coeff[:, nocc:].T @ xop @ coeff[:, nocc:],
                coeff[:, nocc:].T @ overlap @ coeff[:, nocc:])[1]
      occ = coeff[:, :nocc] @ uo; vir = coeff[:, nocc:] @ uv
      coeff = np.hstack((occ, vir))
      localization = {'status':'block_boys','occupied_virtual_separate':True,
        'orthogonality_error':float(np.max(np.abs(coeff.T @ overlap @ coeff-np.eye(orbitals))))}
    h1 = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, coeff), orbitals).reshape((orbitals,)*4)
    raw = {}
    for p in range(modes):
      for q in range(modes):
        value = float(h1[p//2,q//2]) if p%2 == q%2 else 0.
        if value: raw[((1,p),(0,q))] = F.from_float(value)
        for r in range(modes):
          for s in range(modes):
            value=(eri[p//2,r//2,q//2,s//2]*(p%2==r%2)*(q%2==s%2)-eri[p//2,s//2,q//2,r//2]*(p%2==s%2)*(q%2==r%2))/4
            if value: raw[((1,p),(1,q),(0,s),(0,r))] = F.from_float(float(value))
    raw=canonical(raw); sym=canonical(adj(raw))
    h={w:F(round((raw.get(w,F(0))+sym.get(w,F(0)))/2*denominator),denominator)
       for w in set(raw)|set(sym) if round((raw.get(w,F(0))+sym.get(w,F(0)))/2*denominator)}
    if not hermitian(h): raise RuntimeError('rational export lost Hermiticity')
    error=sum(abs(h.get(w,F(0))-raw.get(w,F(0))) for w in set(h)|set(raw))
    fixture={'kind':'straight_hydrogen_chain_fixture_v1','natoms':natoms,'modes':modes,'particles':particles,
      'geometry':geometry,'unit':'Angstrom','basis':'STO-3G','orbital_basis':'RHF canonical MO',
      'sector_dimension':comb(modes,particles),'hamiltonian':encode(h),'coefficient_denominator':denominator,
      'raw_float_coefficient_l1_difference':rational_text(error),'nuclear_repulsion':str(mol.energy_nuc()),
      'pyscf_version':__import__('pyscf').__version__,'localized_basis':localization}
    oracle=DeterminantOracle(fixture); hf=(1<<particles)-1
    hfq=oracle.action(hf).get(hf,F(0)); witness={'states':[hf],'amplitudes':[1]}
    if abs(hfq - F.from_float(float(mf.e_tot-mol.energy_nuc()))) > error + F(1,10**9):
      raise RuntimeError('exact HF replay disagrees with RHF electronic energy')
    result={'status':'fixture_and_hf_upper','natoms':natoms,'modes':modes,'particles':particles,
      'sector_dimension':comb(modes,particles),'hamiltonian_terms':len(h),'hf_electronic_upper':rational_text(hfq),
      'hf_float_electronic':float(mf.e_tot-mol.energy_nuc()),'nuclear_repulsion':float(mol.energy_nuc()),
      'raw_coefficient_l1_difference':rational_text(error),'elapsed_seconds':time.monotonic()-started,
      'scope':'Finite-basis rational electronic Hamiltonian; HF is a variational upper witness only.'}
    if fci_control:
      from pyscf import fci
      solver=fci.direct_spin1.FCI(mol); solver.conv_tol=1e-12
      e,_=solver.kernel(h1,eri,orbitals,(particles//2,particles//2),ecore=0.)
      result.update(fci_electronic=float(e), hf_fci_gap=float(hfq)-float(e), fci_scope='validation only')
    output.mkdir(parents=True,exist_ok=True)
    (output/'fixture.json').write_text(json.dumps(fixture,indent=2)+'\n')
    (output/'upper.json').write_text(json.dumps({'independent_upper':witness,'upper':rational_text(hfq)},indent=2)+'\n')
    (output/'generation.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output-root',default=str(ROOT/'results/certificate_scaling/active_space_ladder'))
    ap.add_argument('--atoms',type=int,nargs='+',default=[4,6,8,10]); ap.add_argument('--fci-control',action='store_true')
    ap.add_argument('--localized',action='store_true')
    a=ap.parse_args(); root=Path(a.output_root)
    root.mkdir(parents=True,exist_ok=True)
    allr=[]
    for n in a.atoms:
      out=root/f'h{n}';
      if out.exists(): raise ValueError(f'preserve existing output: {out}')
      allr.append(build(n,out,fci_control=a.fci_control and n<=6, localized=a.localized))
    (root/'ladder.json').write_text(json.dumps({'kind':'straight_hydrogen_chain_ladder_v1','spacing_angstrom':1.4,'entries':allr},indent=2)+'\n')
    print(json.dumps(allr,indent=2))
if __name__=='__main__': main()
