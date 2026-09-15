"""Numerical FCI reference after discovery, never a constructor input."""
import argparse
from fractions import Fraction as F
import json
from math import comb
from pathlib import Path
import time
import numpy as np
from pyscf import fci
from research.transfer_solver_20260915.budget import dump


def run(case):
    started = time.monotonic()
    data = json.loads((case/'fixture.json').read_text())
    z = np.load(case/'integrals.npz')
    norb, n = data['modes']//2, data['particles']
    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 1e-10
    solver.max_cycle = 200
    solver.max_memory = 2000
    t = time.monotonic()
    e, vector = solver.kernel(z['h1'], z['eri'], norb, (n//2, n//2), ecore=0.)
    solve_seconds = time.monotonic()-t
    s2, multiplicity = fci.spin_op.spin_square(vector, norb, (n//2, n//2))
    record = {'method': 'PySCF direct_spin1 FCI', 'numerical_electronic_energy_Ha': float(e),
              'converged': bool(solver.converged), 'convergence_tolerance_Ha': 1e-10,
              'spin_squared': float(s2), 'spin_multiplicity': float(multiplicity),
              'numerical_solve_seconds': solve_seconds, 'total_seconds': time.monotonic()-started,
              'determinant_coefficients': int(vector.size),
              'expected_Ms0_dimension': comb(norb, n//2)**2,
              'source': 'Newly generated spatial integrals; numerical reference for the pre-rounding model',
              'coefficient_rounding_l1_Ha': data['floating_input_coefficient_rounding_l1_Ha'],
              'rigorous_two_sided_certificate': False,
              'never_used_in_state_or_certificate_discovery': True,
              'accuracy_matched_runtime_claim': False,
              'timing_scope': 'Fixed tight numerical convergence, reported separately from rigorous bounds'}
    if (case/'interval.json').exists():
        result = json.loads((case/'interval.json').read_text())
        lower = json.loads((case/'lower.json').read_text())
        L, U = float(F(result['lower_Ha'])), float(F(result['upper_Ha']))
        rounding = float(F(data['floating_input_coefficient_rounding_l1_Ha']))
        allowance = rounding+1e-7
        record.update(numerical_reference_inside_certificate_with_rounding_and_1e_minus_7_Ha_slack=L-allowance <= e <= U+allowance,
                      upper_minus_reference_mHa=1000*(U-e), reference_minus_lower_mHa=1000*(e-L),
                      certified_width_mHa=result['width_mHa'],
                      singlet_remainder_mHa=1000*float(F(lower['singlet']['residual_l1'])),
                      nonsinglet_lower_Ha=float(F(lower['all_nonsinglets']['lower'])),
                      singlet_lower_Ha=float(F(lower['singlet']['lower'])))
    dump(case/'fci_reference.json', record)
    print(json.dumps(record, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    run(p.parse_args().case.resolve())
