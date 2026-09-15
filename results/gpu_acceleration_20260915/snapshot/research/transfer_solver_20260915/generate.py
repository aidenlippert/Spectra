"""Build the declared molecular model from new integrals, with no FCI input."""
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import time
import numpy as np
from pyscf import gto, scf, ao2mo
from experiments.marginal_symbolic import canonical, adj, encode, hermitian
from experiments.marginal_determinant_tree import DeterminantOracle
from research.transfer_solver_20260915.budget import dump


def build(case):
    started = time.monotonic()
    spec = json.loads((case/'specification.json').read_text())
    mol = gto.M(atom=spec['geometry'], basis=spec['basis'], unit='Angstrom',
                spin=0, charge=spec.get('charge', 0), verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-12, verbose=0)
    if not mf.converged:
        raise RuntimeError('RHF did not converge')
    norb = mol.nao_nr()
    modes, particles = 2*norb, mol.nelectron
    h1 = mf.mo_coeff.T@mf.get_hcore()@mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), norb).reshape((norb,)*4)
    raw = {}
    for p in range(modes):
        for q in range(modes):
            if p % 2 == q % 2 and h1[p//2, q//2]:
                raw[((1, p), (0, q))] = F.from_float(float(h1[p//2, q//2]))
            for r in range(modes):
                for s in range(modes):
                    value = (eri[p//2, r//2, q//2, s//2]*(p % 2 == r % 2)*(q % 2 == s % 2)
                             - eri[p//2, s//2, q//2, r//2]*(p % 2 == s % 2)*(q % 2 == r % 2))/4
                    if value:
                        raw[((1, p), (1, q), (0, s), (0, r))] = F.from_float(float(value))
    raw = canonical(raw)
    conjugate = canonical(adj(raw))
    den = 10**12
    h = {w: F(round((raw.get(w, F(0))+conjugate.get(w, F(0)))*den/2), den)
         for w in set(raw) | set(conjugate)}
    h = {w: c for w, c in h.items() if c}
    if not hermitian(h):
        raise AssertionError('Hermitian rational export failed')
    rounding = sum(abs(h.get(w, F(0))-raw.get(w, F(0))) for w in set(h) | set(raw))
    data = {'kind': 'fresh_molecular_transfer_fixture_v1', 'modes': modes, 'particles': particles,
            'geometry': spec['geometry'], 'basis': spec['basis'], 'unit': 'Angstrom',
            'orbital_basis': 'RHF canonical MO', 'hamiltonian': encode(h),
            'nuclear_repulsion': str(mol.energy_nuc()), 'coefficient_denominator': den,
            'floating_input_coefficient_rounding_l1_Ha': str(rounding),
            'integral_evaluation_error_certified': False, 'sector_dimension': comb(modes, particles)}
    hf = (1 << particles)-1
    hf_energy = DeterminantOracle(data).action(hf).get(hf, F(0))
    if abs(hf_energy-F.from_float(float(mf.e_tot-mol.energy_nuc()))) > rounding+F(1, 10**9):
        raise ValueError('Exact rational HF and integral energies disagree')
    dump(case/'fixture.json', data)
    np.savez_compressed(case/'integrals.npz', h1=h1, eri=eri, mo=mf.mo_coeff)
    dump(case/'generation.json', {'seconds': time.monotonic()-started, 'modes': modes, 'particles': particles,
                                  'hamiltonian_terms': len(h), 'HF_upper_Ha': str(hf_energy),
                                  'RHF_electronic_Ha': float(mf.e_tot-mol.energy_nuc()),
                                  'floating_coefficient_rounding_l1_Ha': str(rounding),
                                  'FCI_used': False, 'HF_determinants_evaluated': 1,
                                  'full_sector_enumerated': False, 'pyscf_version': __import__('pyscf').__version__})


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    build(p.parse_args().case.resolve())
