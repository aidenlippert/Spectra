"""Small physical-model validation experiment, separate from exact certificates."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pyscf
import scipy
from pyscf import ao2mo, fci, gto, lib, scf
from scipy.optimize import minimize_scalar


def point(basis, distance):
    start = time.monotonic()
    mol = gto.M(atom=[('H', (0, 0, 0)), ('H', (0, 0, distance))],
                basis=basis, unit='Angstrom', charge=0, spin=0, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    if not mf.converged:
        raise RuntimeError('RHF did not converge')
    orbitals = mf.mo_coeff
    norb = orbitals.shape[1]
    h1 = orbitals.T @ mf.get_hcore() @ orbitals
    h2 = ao2mo.kernel(mol, orbitals)
    solver = fci.direct_spin0.FCI(mol)
    solver.conv_tol = 1e-13
    solver.conv_tol_residual = 1e-10
    solver.max_cycle = 200
    # This two-electron validation has at most 784 determinants. Full diagonalization
    # avoids an iterative-convergence failure without weakening residual acceptance.
    if norb**2 > 1000:
        raise ValueError('This dense validation runner is limited to 1000 determinants')
    solver.pspace_size = norb**2
    energy, vector = solver.kernel(h1, h2, norb, (1, 1))
    if not solver.converged:
        raise RuntimeError('Singlet FCI did not converge')
    effective = fci.direct_spin1.absorb_h1e(h1, h2, norb, (1, 1), .5)
    applied = fci.direct_spin1.contract_2e(effective, vector, norb, (1, 1))
    norm = float(np.linalg.norm(vector))
    residual = float(np.linalg.norm(applied - energy * vector) / norm)
    spin_squared = float(fci.spin_op.spin_square0(vector, norb, (1, 1))[0])
    if not np.isfinite(residual) or residual > 1e-8 or abs(norm - 1) > 1e-10 or abs(spin_squared) > 1e-8:
        raise RuntimeError('FCI numerical residual or normalization gate failed')
    return {'distance_angstrom': float(distance), 'basis': basis, 'spatial_orbitals': norb,
            'fixed_Ms_determinants': norb**2, 'rhf_converged': bool(mf.converged),
            'fci_converged': bool(solver.converged), 'electronic_energy_hartree': float(energy),
            'nuclear_repulsion_hartree': float(mol.energy_nuc()),
            'total_energy_hartree': float(energy + mol.energy_nuc()),
            'fci_residual_norm_hartree': residual, 'vector_norm': norm,
            'spin_squared': spin_squared, 'full_dense_diagonalization': True,
            'wall_seconds': time.monotonic() - start}


def run(reference_path, out):
    start = time.monotonic()
    lib.num_threads(1)
    reference = json.loads(reference_path.read_text())
    if (reference['species'], reference['quantity'], reference['unit']) != (
            'H2', 'equilibrium internuclear distance r_e', 'angstrom'):
        raise ValueError('Reference species, quantity or unit mismatch')
    out.mkdir(parents=True, exist_ok=False)
    results = []
    for basis in ('sto-3g', 'cc-pvdz', 'cc-pvtz'):
        history = []

        def objective(distance):
            record = point(basis, distance)
            history.append(record)
            return record['total_energy_hartree']

        optimum = minimize_scalar(objective, bounds=(.60, .95), method='bounded',
                                  options={'xatol': 1e-7, 'maxiter': 30})
        if not optimum.success:
            raise RuntimeError('Bounded geometry minimization did not converge')
        center = point(basis, optimum.x)
        neighbors = [point(basis, optimum.x + shift) for shift in (-.001, .001)]
        if not all(p['total_energy_hartree'] > center['total_energy_hartree'] for p in neighbors):
            raise RuntimeError('Optimized geometry is not below its tested neighbors')
        delta = float(optimum.x) - float(reference['value'])
        record = {'basis': basis, 'minimum': center, 'neighbor_checks': neighbors,
                  'optimizer_success': bool(optimum.success), 'optimizer_evaluations': int(optimum.nfev),
                  'difference_from_nist_angstrom': delta,
                  'absolute_difference_from_nist_angstrom': abs(delta),
                  'numerical_checks_passed': True, 'experimental_agreement_pass': None,
                  'agreement_note': 'Difference reported without a pass threshold; source uncertainty unspecified.',
                  'evaluations': history}
        results.append(record)
        (out / (basis + '.json')).write_text(json.dumps(record, indent=2) + '\n')
        print(json.dumps({k: record[k] for k in ('basis', 'difference_from_nist_angstrom',
                                                'numerical_checks_passed')}), flush=True)
    receipt = {'reference': reference,
               'reference_sha256': hashlib.sha256(reference_path.read_bytes()).hexdigest(),
               'pyscf': pyscf.__version__, 'scipy': scipy.__version__, 'numpy': np.__version__,
               'threads': 1, 'method': 'Nonrelativistic Born-Oppenheimer singlet numerical FCI geometry',
               'results': results, 'wall_seconds': time.monotonic() - start,
               'scope': 'One species and property in three finite bases. Numerical FCI validation is separate from exact rational certification; no force certificate or experimental error bar is produced.'}
    (out / 'summary.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('DONE', receipt['wall_seconds'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', type=Path, default=Path(__file__).with_name('nist_h2_reference.json'))
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.reference, args.out)
