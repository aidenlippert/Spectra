"""Bounded numerical geometry/basis control; never an exact physical certificate."""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/transfer_solver_20260915/ch2_model_study'
BOHR_ANGSTROM = 0.52917721092
CASES = [(basis, spin) for basis in ('cc-pvdz', 'cc-pvtz') for spin in (0, 1)]


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def coordinates(r, theta):
    y, z = r*math.sin(theta/2), r*math.cos(theta/2)
    return [('C', (0., 0., 0.)), ('H', (0., y, z)), ('H', (0., -y, z))]


def chain_gradient(r, theta, gradient):
    dy = gradient[1][1]-gradient[2][1]
    dz = gradient[1][2]+gradient[2][2]
    return [(dy*math.sin(theta/2)+dz*math.cos(theta/2))/BOHR_ANGSTROM,
            r*(dy*math.cos(theta/2)-dz*math.sin(theta/2))/(2*BOHR_ANGSTROM)]


def evaluate(basis, spin, x, log):
    import numpy as np
    from pyscf import gto, scf, mcscf, fci, symm
    from pyscf.mcscf import avas
    started = time.monotonic()
    mol = gto.M(atom=coordinates(*x), basis=basis, unit='Angstrom', spin=2*spin,
                charge=0, symmetry='C2v', verbose=0, max_memory=1800)
    # PySCF can reorient symmetry axes. Identify the out-of-plane p orbital
    # from the actual AO transformation rather than assume its B1/B2 label.
    identity = np.eye(mol.nao_nr())
    px, py = [int(mol.search_ao_label('C 2p'+axis)[0]) for axis in ('x', 'y')]
    perpendicular, inplane = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, identity[:, [px, py]])
    if {perpendicular, inplane} != {'B1', 'B2'}:
        raise ValueError(('Unexpected planar valence symmetry labels', perpendicular, inplane))
    mf = (scf.ROHF(mol) if spin else scf.RHF(mol))
    mf.irrep_nelec = {'A1': (3, 2), perpendicular: (1, 0), inplane: (1, 1)} if spin else {'A1': 6, inplane: 2}
    mf.conv_tol = 1e-11
    mf.kernel()
    if not mf.converged:
        raise ValueError('SCF did not converge')
    ncas, nelecas, mo = avas.avas(mf, ['C 2s', 'C 2p', 'H 1s'], ncore=1,
                                openshell_option=3, verbose=0)
    if (ncas, nelecas, mol.nelectron) != (6, 6, 8):
        raise ValueError(('Unexpected valence active space', ncas, nelecas, mol.nelectron))
    mc = mcscf.CASSCF(mf, 6, (3+spin, 3-spin))
    mc.ncore = 1
    mc.conv_tol = 1e-10
    mc.conv_tol_grad = 1e-5
    mc.max_cycle_macro = 60
    mc.fcisolver.conv_tol = 1e-12
    mc.fcisolver.wfnsym = perpendicular if spin else 'A1'
    mc.kernel(mo)
    ss, multiplicity = fci.spin_op.spin_square(mc.ci, 6, mc.nelecas)
    if not mc.converged or abs(float(ss)-spin*(spin+1)) > 1e-7:
        raise ValueError(('CASSCF convergence or total-spin check failed', mc.converged, float(ss)))
    gradient = mc.nuc_grad_method().kernel()
    gx = np.asarray(chain_gradient(float(x[0]), float(x[1]), gradient))
    row = {'bond_Angstrom': float(x[0]), 'angle_degrees': math.degrees(x[1]),
           'CASSCF_energy_Ha': float(mc.e_tot), 'measured_S2': float(ss), 'out_of_plane_p_irrep': perpendicular,
           'parameter_gradient': gx.tolist(), 'seconds': time.monotonic()-started}
    with log.open('a') as stream:
        stream.write(json.dumps(row)+'\n')
    print(json.dumps(row), flush=True)
    return float(mc.e_tot), gx, mc


def calculate(basis, spin):
    import numpy as np
    from scipy.optimize import minimize
    from pyscf import mrpt, ao2mo
    directory = OUT/f'{basis}_S{spin}'
    directory.mkdir(parents=True, exist_ok=False)
    log = directory/'optimization.jsonl'
    started = time.monotonic()
    cache = {}

    def fun(x):
        key = tuple(map(float, x))
        if key not in cache:
            e, grad, mc = evaluate(basis, spin, x, log)
            cache.clear()
            cache[key] = e, grad, mc
        return cache[key][:2]

    initial = [1.08, math.radians(134)] if spin else [1.11, math.radians(102)]
    opt = minimize(fun, initial, jac=True, method='L-BFGS-B',
                   bounds=[(.9, 1.3), (math.radians(85), math.radians(145))],
                   options={'ftol': 1e-12, 'gtol': 2e-5, 'maxiter': 30, 'maxfun': 50})
    e, grad = fun(opt.x)
    mc = cache[tuple(map(float, opt.x))][2]
    if not opt.success or max(abs(grad)) > 1e-4:
        dump(directory/'failed_geometry.json', {'message': str(opt.message), 'gradient': list(map(float, grad))})
        raise ValueError('Geometry minimization did not meet frozen gradient tolerance')
    correction_start = time.monotonic()
    correction = float(mrpt.NEVPT(mc).kernel())
    h1, core = mc.get_h1eff(mc.mo_coeff)
    eri = ao2mo.restore(1, mc.get_h2eff(mc.mo_coeff), 6)
    np.savez(directory/'active_integrals.npz', h1=h1, eri=eri, core=float(core), mo=mc.mo_coeff)
    dump(directory/'result.json', {'basis': basis, 'total_spin': spin, 'pyscf_state_symmetry': str(mc.fcisolver.wfnsym),
        'physical_state': 'triplet with one out-of-plane p electron' if spin else 'A1 singlet',
        'physical_electrons': 8, 'active_electrons': 6, 'active_spatial_orbitals': 6,
        'CASSCF_inactive_doubly_occupied_orbitals': 1,
        'active_Ms_determinant_dimension': 400 if spin == 0 else 225,
        'active_space_FCI_used': True, 'full_parent_basis_FCI_used': False,
        'total_spatial_orbitals': mc.mo_coeff.shape[1], 'geometry_Angstrom': coordinates(*map(float, opt.x)),
        'bond_Angstrom': float(opt.x[0]), 'angle_degrees': math.degrees(opt.x[1]),
        'geometry_gradient': list(map(float, grad)), 'CASSCF_total_Ha': e,
        'SC_NEVPT2_correction_Ha': correction, 'SC_NEVPT2_total_Ha': e+correction,
        'NEVPT2_seconds': time.monotonic()-correction_start, 'complete_internal_seconds': time.monotonic()-started,
        'geometry_evaluations': sum(1 for _ in log.open()), 'geometry_converged': bool(opt.success),
        'spin_checked_at_every_evaluation': True, 'pyscf_version': __import__('pyscf').__version__,
        'scope': 'Numerical CASSCF-optimized geometries and SC-NEVPT2 single points. No rigorous total-energy bound, zero-point correction, or physical error interval.'})


def run():
    source = Path(__file__)
    protocol = {'frozen_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'cases': CASES,
        'method': 'CASSCF(6,6) C2v geometry optimization; SC-NEVPT2 at those geometries',
        'comparison_observable': 'E_singlet - E_triplet, approximate adiabatic electronic separation',
        'experimental_inferred_Te_kcal_mol': [8.7, .5], 'experimental_raw_splitting_kcal_mol': [9., .09],
        'experimental_source': 'https://experts.umn.edu/en/publications/methylene-a-study-of-the-xsup3supbsub1sub-and-%C3%A3sup1supasub1sub-st/',
        'no_certified_physical_error_bars': True, 'no_prospective_predictions': True,
        'basis_difference_is_sensitivity_not_error_bound': True, 'per_case_timeout_seconds': 360}
    dump(OUT/'protocol.json', protocol)
    outcomes = []
    for basis, spin in CASES:
        if hashlib.sha256(source.read_bytes()).hexdigest() != protocol['source_sha256']:
            raise ValueError('Frozen physical-model source changed')
        cmd = ['/opt/homebrew/Caskroom/miniconda/base/bin/python', '-B', '-S', '-m',
               'research.transfer_solver_20260915.budget', '--name', f'ch2_{basis}_S{spin}', '--seconds', '360', '--',
               str(ROOT/'.venv-molecule/bin/python'), '-B', '-m', 'research.ch2_model_study_20260915.study',
               'case', '--basis', basis, '--spin', str(spin)]
        result = subprocess.run(cmd, cwd=ROOT)
        outcomes.append({'basis': basis, 'spin': spin, 'exit_code': result.returncode})
    dump(OUT/'execution.json', outcomes)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('run', 'case'))
    parser.add_argument('--basis', choices=('cc-pvdz', 'cc-pvtz'))
    parser.add_argument('--spin', type=int, choices=(0, 1))
    args = parser.parse_args()
    run() if args.action == 'run' else calculate(args.basis, args.spin)
