"""Nonaccepting numerical conditioning probe for failed exact basis exports."""
from pathlib import Path
from fractions import Fraction as F
import gzip
import hashlib
import json
import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/spin_word'


def main():
    cases = {}
    files = [Path(__file__).resolve()]
    for case in ('W_zero', 'W_plus_1'):
        path = BASE / case / 'refined/candidate_ledger.json.gz'
        files.append(path)
        with gzip.open(path, 'rt') as stream:
            data = json.load(stream)
        if not 1 <= len(data['states']) <= 13000:
            raise ValueError('Bounded ledger required')
        matrix = np.array(data['rows']).T
        rhs = np.array(list(map(F, data['rhs'])), float)
        energy = np.array(data['energies'])
        if matrix.shape != (198, len(energy)) or not np.isfinite(matrix).all():
            raise ValueError('Finite198-row ledger required')
        trials = []
        for scale in (1, 1000, 1000000):
            result = linprog(energy, A_eq=scale * matrix[:196], b_eq=scale * rhs[:196],
                             A_ub=scale * matrix[196:], b_ub=scale * rhs[196:], bounds=(0, None),
                             method='highs-ds', options={'presolve': False,
                             'dual_feasibility_tolerance': 1e-10, 'primal_feasibility_tolerance': 1e-10})
            trial = {'scale': scale, 'success': bool(result.success), 'message': result.message}
            if result.success:
                values = np.r_[result.x, rhs[196:] - matrix[196:] @ result.x]
                trial.update(numerical_upper=float(result.fun / 5),
                             support_above_zero=int(sum(values > 0)),
                             support_above_1e_minus14=int(sum(values > 1e-14)),
                             maximum_equality_residual=float(max(abs(matrix[:196] @ result.x - rhs[:196]))),
                             smallest_positive_values=sorted(map(float, values[values > 0]))[:12])
            trials.append(trial)
        cases[case] = trials
    result = {'accepted': False, 'cases': cases,
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
              'scope': 'Numerical conditioning evidence only. Small floating residuals and LP success never replace exact rational consistency or positive-weight physical replay.'}
    (BASE / 'lp_scale_probe.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'accepted': False, 'conditioning_probe_recorded': True}))


if __name__ == '__main__':
    main()
