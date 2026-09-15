"""Explicit small spin matrices: charged diagnostic, never the accepting rule."""
from fractions import Fraction as F
from itertools import combinations
import json
import time
import numpy as np

from research.response_consistency_20260913.response import OUT, parameters, check, propose, cases


def matrix(model, couplings=None):
    M, q, _, _, g, _ = parameters(model)
    epsilon, omega, chi = (F(model[k]) for k in ('epsilon', 'omega', 'chi'))
    couplings = [g]*M if couplings is None else couplings
    states = [sum(1 << i for i in occupied) for occupied in combinations(range(M+1), q)]
    lookup = {s: i for i, s in enumerate(states)}; H = np.zeros((len(states), len(states)))
    for j, state in enumerate(states):
        central = state & 1; bath = q-central
        H[j, j] = float(epsilon*central+omega*bath+chi*central*bath)
        for i, coupling in enumerate(couplings, 1):
            if ((state >> i) & 1) != central:
                H[lookup[state ^ 1 ^ (1 << i)], j] = float(coupling)
    return H


def run():
    start = time.monotonic(); records = []
    for model in cases()[:4]:
        before = time.monotonic(); cert = propose(model); accepted = check(model, cert)
        H = matrix(model); energy = float(np.linalg.eigvalsh(H)[0])
        agrees = float(F(accepted['lower']))-1e-12 <= energy <= float(F(accepted['upper']))+1e-12
        if not agrees:
            raise AssertionError('Enumerated diagnostic conflicts with the exact interval')
        records.append({'bath_spins': model['bath_spins'], 'matrix_dimension': len(H), 'numerical_ground_energy': energy,
            'within_exact_interval_up_to_1e_12_float_tolerance': agrees, 'wall_seconds': time.monotonic()-before})
    model = cases()[2]; cert = propose(model); M, q, A0, D0, g, c = parameters(model)
    couplings = [g]*M; couplings[0] = 3*g/2; average = sum(couplings)/M; t = F(cert['trial_t'])
    # The same symmetric trial gives an exact upper for inhomogeneous couplings.
    trial = (A0+t*t*c*D0+2*average*t*c)/(1+t*t*c)
    if trial >= F(cert['lower']):
        raise AssertionError('Intended falsification did not separate')
    try:
        check({**model, 'site_couplings': [str(g) for g in couplings]}, cert)
    except ValueError as error:
        refusal = str(error)
    else:
        raise AssertionError('Homogeneous checker accepted an inhomogeneous input')
    H = matrix(model, couplings)
    changed = {'bath_spins': M, 'couplings': [str(x) for x in couplings], 'matrix_dimension': len(H),
        'numerical_ground_energy': float(np.linalg.eigvalsh(H)[0]), 'exact_symmetric_trial_upper': str(trial),
        'homogeneous_lower': cert['lower'], 'homogeneous_claim_falsified_by_exact_upper': True,
        'homogeneous_checker_refusal': refusal}
    result = {'homogeneous_controls': records, 'fresh_inhomogeneous_counterexample': changed,
        'many_body_basis_states_enumerated': sum(r['matrix_dimension'] for r in records)+len(H),
        'wall_seconds': time.monotonic()-start,
        'scope': 'Numerical spectra are small diagnostic controls. Inhomogeneous invalidity is proved by an exact rational upper, not by floating eigenvalues.'}
    (OUT/'response_controls.json').write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result), flush=True)


if __name__ == '__main__':
    run()
