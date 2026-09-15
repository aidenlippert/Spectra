"""Charged finite-sector spectral checks; never used to accept a certificate."""
import json
import time
import numpy as np

from research.positive_response_20260913.block_response import OUT, parameters, basis, exchange, lowering, F


def run():
    start = time.monotonic(); rows = []
    records = json.loads((OUT/'response_certificates.json').read_text())
    for record in records:
        M, q, A0, D0, g, _, l0, l1 = parameters(record['model'])
        if M > 9:
            continue
        before = time.monotonic(); left = basis(M, q); right = basis(M, q-1)
        locations = [{state: i for i, state in enumerate(left)}, {state: len(left)+i for i, state in enumerate(right)}]
        H = np.diag([float(A0)]*len(left)+[float(D0)]*len(right))
        for sector, states, strength in ((0, left, l0), (1, right, l1)):
            for state in states:
                j = locations[sector][state]
                for target, amplitude in exchange(M, {state: F(1)}).items():
                    H[locations[sector][target], j] += float(strength*amplitude)
        for state in left:
            j = locations[0][state]
            for target, amplitude in lowering(M, {state: F(1)}).items():
                i = locations[1][target]; H[i, j] += float(g*amplitude); H[j, i] += float(g*amplitude)
        if not np.array_equal(H, H.T):
            raise AssertionError('Explicit control Hamiltonian is not symmetric')
        energy = float(np.linalg.eigvalsh(H)[0]); b = float(F(record['certificate']['lower']))
        from research.positive_response_20260913.block_response import check
        receipt = check(record['model'], record['certificate']); upper = float(F(receipt['upper']))
        tolerance = 1e-10
        if not b-tolerance <= energy <= upper+tolerance:
            raise AssertionError('Explicit spectrum disagrees with the analytic interval')
        rows.append({'bath_spins': M, 'excitations': q, 'basis_states_enumerated': len(H),
            'full_many_body_matrix_entries': int(H.size), 'numerical_ground_energy': energy,
            'difference_from_exact_upper': energy-upper, 'comparison_tolerance': tolerance,
            'wall_seconds': time.monotonic()-before})
    result = {'cases': rows, 'wall_seconds': time.monotonic()-start,
        'basis_states_enumerated': sum(row['basis_states_enumerated'] for row in rows),
        'scope': 'Numerical finite-sector control only; all matrices and enumeration are charged. No molecular FCI or SDP.'}
    (OUT/'response_numeric_controls.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    run()
