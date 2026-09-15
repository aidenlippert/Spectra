"""Conventional all-electron CCSD(T) comparison at the existing CASSCF geometries.

These are numerical, nonvariational single points, not Spectra certificates or
CCSD(T)-optimized adiabatic energies. No measured gap is an optimization input.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, CHEM, PREFIX

BASE = OUT/'ch2_model_study'


def case(basis, spin):
    import numpy as np
    from pyscf import gto, scf, cc, symm
    start = time.monotonic()
    source = BASE/f'{basis}_S{spin}/result.json'
    data = json.loads(source.read_text())
    mol = gto.M(atom=data['geometry_Angstrom'], basis=basis, unit='Angstrom',
        charge=0, spin=2*spin, symmetry='C2v', max_memory=1800, verbose=0)
    if mol.nelectron != 8 or mol.nao_nr() != data['total_spatial_orbitals']:
        raise ValueError('Molecular/basis accounting changed')
    px, py = [int(mol.search_ao_label('C 2p'+a)[0]) for a in ('x', 'y')]
    perpendicular, inplane = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, np.eye(mol.nao_nr())[:, [px, py]])
    if {perpendicular, inplane} != {'B1', 'B2'}:
        raise ValueError('Unexpected planar state symmetry')
    mf = scf.ROHF(mol) if spin else scf.RHF(mol)
    mf.irrep_nelec = {'A1': (3, 2), perpendicular: (1, 0), inplane: (1, 1)} if spin else {'A1': 6, inplane: 2}
    mf.conv_tol = 1e-11
    mf.kernel()
    if not mf.converged:
        raise ValueError('Reference SCF did not converge')
    s2 = float(mf.spin_square()[0])
    if abs(s2-spin*(spin+1)) > 1e-8:
        raise ValueError('Reference has unexpected spin')
    reference = mf.to_uhf() if spin else mf
    model = cc.CCSD(reference, frozen=0)
    model.conv_tol = 1e-10
    model.conv_tol_normt = 1e-8
    model.max_cycle = 100
    model.kernel()
    if not model.converged:
        raise ValueError('CCSD did not converge')
    triples = float(model.ccsd_t())
    dump(BASE/f'{basis}_S{spin}/cc_control.json', {
        'basis': basis, 'target_spin': spin, 'geometry_source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'method': type(model).__name__+'(T)',
        'reference': 'ROHF orbitals converted to unrestricted storage' if spin else 'RHF',
        'reference_spin_squared': s2, 'CC_state_spin_purity_certified': False,
        'correlated_electrons': 8, 'frozen_orbitals': 0,
        'reference_Ha': float(mf.e_tot), 'CCSD_total_Ha': float(model.e_tot),
        'triples_correction_Ha': triples, 'CCSD_T_total_Ha': float(model.e_tot)+triples,
        'converged': bool(model.converged), 'seconds': time.monotonic()-start,
        'geometry_optimized_for': 'CASSCF(6,6), not CCSD(T)',
        'scope': 'Numerical conventional-method comparison; nonvariational, no rigorous total-energy or physical uncertainty interval',
        'source': 'https://pyscf.org/user/cc.html'})


def run():
    steps = [(b, s) for b in ('cc-pvdz', 'cc-pvtz') for s in (0, 1)]
    dump(BASE/'cc_control_protocol.json', {
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'cases': steps, 'per_case_seconds': 240,
        'geometry_cost_from_prior_CASSCF_is_additional': True,
        'all_electrons_correlated': True, 'matching_geometries_to_CASSCF_NEVPT2': True,
        'not_a_prospective_or_certified_physical_prediction': True})
    outcomes = []
    for basis, spin in steps:
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', f'ch2_cc_{basis}_s{spin}',
            '--seconds', '240', '--', CHEM, '-B', '-m', 'research.ch2_model_study_20260915.cc_control',
            'case', '--basis', basis, '--spin', str(spin)], cwd=ROOT)
        outcomes.append({'basis': basis, 'spin': spin, 'exit_code': result.returncode})
    dump(BASE/'cc_control_execution.json', outcomes)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'case'))
    p.add_argument('--basis', choices=('cc-pvdz', 'cc-pvtz'))
    p.add_argument('--spin', type=int, choices=(0, 1))
    a = p.parse_args()
    run() if a.action == 'run' else case(a.basis, a.spin)
