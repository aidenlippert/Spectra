"""Neutral CH2 CAS(6,6) from seven STO-3G orbitals, with full core terms."""
from fractions import Fraction as F
from pathlib import Path
import json
import time

from experiments.marginal_symbolic import add, mono, canonical, encode, hermitian

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/ch2_validation_20260913'


def run():
    import numpy as np
    from pyscf import gto, scf, mcscf, ao2mo, fci
    start = time.monotonic()
    geometry = 'C 0 0 0; H 0 0 1.117; H 0 1.047 -0.389'
    mol = gto.M(atom=geometry, basis='sto-3g', unit='Angstrom', charge=0, spin=0, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-12, verbose=0)
    if not mf.converged or mol.nelectron != 8 or mf.mo_coeff.shape[1] != 7:
        raise ValueError('Neutral CH2 RHF input did not match the declared model')
    mc = mcscf.CASCI(mf, 6, 6); mc.ncore = 1
    h1, core = mc.get_h1eff(mf.mo_coeff)
    eri = ao2mo.restore(1, mc.get_h2eff(mf.mo_coeff), 6)
    # Spin-independent spatial rounding, including all ERI symmetries.
    den = 10**12; s = 6
    hq = [[F(round(float((h1[i, j]+h1[j, i])/2)*den), den) for j in range(s)] for i in range(s)]
    pairs = [(i, j) for i in range(s) for j in range(i, s)]; integrals = {}
    for a, (i, j) in enumerate(pairs):
        for k, l in pairs[a:]:
            orbit = {(i,j,k,l),(j,i,k,l),(i,j,l,k),(j,i,l,k),(k,l,i,j),(l,k,i,j),(k,l,j,i),(l,k,j,i)}
            value = F(round(sum(float(eri[x]) for x in orbit)/len(orbit)*den), den)
            for key in orbit:
                integrals[key] = value
    cq = F(round(float(core)*den), den); terms = {(): cq}
    for i in range(s):
        for j in range(s):
            for spin in range(2):
                word = ((1,2*i+spin),(0,2*j+spin)); terms[word] = hq[i][j]
            for k in range(s):
                for l in range(s):
                    for spin in range(2):
                        for other in range(2):
                            word = ((1,2*i+spin),(1,2*k+other),(0,2*l+other),(0,2*j+spin))
                            terms[word] = terms.get(word, F(0))+integrals[i,j,k,l]/2
    h = canonical(terms)
    if not hermitian(h):
        raise ValueError('Spatial rounding lost Hermiticity')
    raw_error = abs(cq-F.from_float(float(core)))+2*sum(abs(hq[i][j]-F.from_float(float(h1[i,j]))) for i in range(s) for j in range(s))
    raw_error += 2*sum(abs(integrals[i,j,k,l]-F.from_float(float(eri[i,j,k,l]))) for i in range(s) for j in range(s) for k in range(s) for l in range(s))
    fixture = {'kind': 'ch2_cas6_6_fixed_geometry_v1', 'modes':12, 'particles':6,
        'physical_electrons':8, 'total_spatial_orbitals':7, 'frozen_core_orbitals':1,
        'active_spatial_orbitals':6, 'geometry':geometry, 'unit':'Angstrom', 'basis':'STO-3G',
        'hamiltonian':encode(h), 'energy_constant_Ha':str(cq), 'constant_includes':'frozen core and nuclear repulsion',
        'nuclear_repulsion_Ha':str(float(mol.energy_nuc())), 'spatial_coefficient_denominator':den,
        'spatial_h1':[[str(x) for x in row] for row in hq],
        'spatial_eri':[str(integrals[i,j,k,l]) for i in range(s) for j in range(s) for k in range(s) for l in range(s)],
        'operator_rounding_envelope_Ha':str(raw_error), 'pyscf_version':__import__('pyscf').__version__,
        'core_method':'PySCF CASCI.get_h1eff and get_h2eff; all frozen-core one-body contributions included.'}
    # Fixed-model numerical comparators; exact certification is a separate run.
    roots = []
    for spin, electrons, solver in [(0,(3,3),fci.direct_spin0.FCI(mol)), (1,(4,2),fci.direct_spin1.FCI(mol))]:
        solver.conv_tol = 1e-12; t = time.monotonic(); e, v = solver.kernel(h1, eri, s, electrons, ecore=float(core))
        ss, multiplicity = fci.spin_op.spin_square(v, s, electrons)
        roots.append({'target_spin':spin,'electrons':list(electrons),'energy_total_Ha':float(e),
            'measured_S2':float(ss),'FCI_seconds':time.monotonic()-t})
    receipt = {'geometry':geometry,'RHF_total_Ha':float(mf.e_tot),'core_constant_Ha':float(core),
        'numerical_FCI_comparators':roots,'operator_rounding_envelope_Ha':str(raw_error),
        'wall_seconds':time.monotonic()-start,'scope':'Physical-model fixture and numerical proposals, not an exact energy certificate.'}
    OUT.mkdir(exist_ok=True); (OUT/'fixture.json').write_text(json.dumps(fixture,indent=2)+'\n')
    (OUT/'generation.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__ == '__main__':
    run()
