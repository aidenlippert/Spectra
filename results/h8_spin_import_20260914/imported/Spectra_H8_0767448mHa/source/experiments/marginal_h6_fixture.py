"""H6 chain transfer fixture with separate numerical FCI and exact upper controls."""
from fractions import Fraction as F
import json
from math import comb
from pathlib import Path
import time

from experiments.marginal_symbolic import canonical, adj, encode, hermitian
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_implicit_certificate import rational_text

ROOT = Path(__file__).resolve().parents[1]


def interleaved_state(alpha, beta, orbitals):
    state = sum((1 << (2 * i)) for i in range(orbitals) if alpha & (1 << i))
    state |= sum((1 << (2 * i + 1)) for i in range(orbitals) if beta & (1 << i))
    crossings = sum(1 for i in range(orbitals) if alpha & (1 << i)
                    for j in range(orbitals) if beta & (1 << j) and i > j)
    return state, (-1 if crossings % 2 else 1)


def run():
    import numpy as np
    import pyscf
    from pyscf import gto, scf, fci, ao2mo
    out = ROOT / 'results/marginal_h6'
    if out.exists():
        raise ValueError('Preserve previous H6 transfer fixture')
    started = time.monotonic()
    geometry = [['H', [0., 0., float(F(7 * i, 5))]] for i in range(6)]
    mol = gto.M(atom=geometry, basis='sto-3g', unit='Angstrom', spin=0, charge=0, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-12, verbose=0)
    if not mf.converged:
        raise ValueError('RHF fixture generation did not converge')
    orbitals = mol.nao_nr()
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), orbitals).reshape((orbitals,) * 4)
    modes, particles = 2 * orbitals, mol.nelectron
    raw = {}
    for p in range(modes):
        for q in range(modes):
            value = float(h1[p // 2, q // 2]) if p % 2 == q % 2 else 0.
            if value:
                raw[((1, p), (0, q))] = F.from_float(value)
            for r in range(modes):
                for s in range(modes):
                    value = (eri[p // 2, r // 2, q // 2, s // 2] * (p % 2 == r % 2) * (q % 2 == s % 2)
                             - eri[p // 2, s // 2, q // 2, r // 2] * (p % 2 == s % 2) * (q % 2 == r % 2)) / 4
                    if value:
                        raw[((1, p), (1, q), (0, s), (0, r))] = F.from_float(float(value))
    raw = canonical(raw)
    adjoint = canonical(adj(raw))
    symmetric = {w: (raw.get(w, F(0)) + adjoint.get(w, F(0))) / 2 for w in set(raw) | set(adjoint)}
    denominator = 10**12
    h = {w: F(round(x * denominator), denominator) for w, x in symmetric.items() if round(x * denominator)}
    if not hermitian(h) or any(sum(2*c - 1 for c, i in w if i % 2 == spin) for w in h for spin in (0, 1)):
        raise ValueError('Rational export lost Hermiticity or separate spin-number conservation')
    error = sum(abs(h.get(w, F(0)) - raw.get(w, F(0))) for w in set(h) | set(raw))
    fixture = {'modes': modes, 'particles': particles, 'hamiltonian': encode(h),
               'geometry': geometry, 'unit': 'Angstrom', 'basis': 'STO-3G', 'orbital_basis': 'RHF canonical MO',
               'sector_dimension': comb(modes, particles), 'coefficient_denominator': denominator,
               'raw_float_coefficient_l1_difference': rational_text(error), 'pyscf_version': pyscf.__version__}
    oracle = DeterminantOracle(fixture)
    hf = (1 << particles) - 1
    hf_rational = oracle.action(hf).get(hf, F(0))
    hf_float = float(mf.e_tot - mol.energy_nuc())
    if abs(float(hf_rational) - hf_float) > float(error) + 1e-9:
        raise ValueError('RHF determinant failed the independent integral-convention check')
    solver = fci.direct_spin1.FCI(mol)
    solver.conv_tol = 1e-12
    solver.max_cycle = 100
    fci_electronic, ci = solver.kernel(h1, eri, orbitals, (particles // 2, particles // 2), ecore=0.)
    if not solver.converged:
        raise ValueError('Independent FCI control did not converge')
    strings = list(fci.cistring.make_strings(range(orbitals), particles // 2))
    vector = {}
    for i, alpha in enumerate(strings):
        for j, beta in enumerate(strings):
            state, sign = interleaved_state(int(alpha), int(beta), orbitals)
            amplitude = int(round(float(ci[i, j]) * sign * 10**10))
            if amplitude:
                vector[state] = amplitude
    states = sorted(vector)
    witness = {'states': states, 'amplitudes': [vector[s] for s in states]}
    reference_upper = oracle.upper(witness)
    if abs(float(reference_upper) - float(fci_electronic)) > float(error) + 1e-8:
        raise ValueError('FCI-to-CAR determinant signs or integral convention failed')
    control = {'status': 'upper_only_reference', 'independent_upper': witness,
               'upper': rational_text(reference_upper), 'electronic_fci_numerical': float(fci_electronic),
               'nuclear_repulsion_numerical': float(mol.energy_nuc()), 'rhf_electronic_numerical': hf_float,
               'rhf_rational_diagonal': rational_text(hf_rational),
               'numerical_fci_difference': abs(float(reference_upper) - float(fci_electronic)),
               'scope': 'Separate enumerative FCI validation control; not supplied to sparse reference selection; no lower endpoint or certified integral/basis error'}
    receipt = {'status': 'fixture_and_upper_control_only', 'modes': modes, 'particles': particles,
               'sector_dimension': comb(modes, particles), 'hamiltonian_terms': len(h),
               'reference_upper_support': len(states), 'elapsed_seconds': time.monotonic() - started}
    out.mkdir(parents=True)
    for filename, value in (('fixture.json', fixture), ('reference_upper.json', control), ('generation.json', receipt)):
        (out / filename).write_text(json.dumps(value, indent=2) + '\n')
    print(json.dumps(dict(receipt, reference_upper=float(reference_upper), numerical_fci=float(fci_electronic))), flush=True)


if __name__ == '__main__':
    run()
