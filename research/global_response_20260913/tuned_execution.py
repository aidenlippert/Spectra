"""Charged numerical execution of the terminal-aware H6 response variant."""
from fractions import Fraction as F
import json
import time
from research.global_response_20260913 import execution as e, global_program as g
from research.compact_response_20260913 import closure


def run():
    import numpy as np
    start = time.monotonic()
    data, tail, _ = g.load_case('h6'); cert = json.loads((g.OUT/'h6_terminal_tuned.json').read_text())
    proof = g.check(data, tail, cert)
    groups, den, preparation = closure.blocks(data)
    group = next(x for x in groups if 63 in x['states']); states = group['states']; n = len(states)
    M = np.array(group['H'], dtype=float)/float(den)-float(F(cert['target_Ha']))*np.eye(n)
    q = np.array([((s >> 10)&3) == 3 or ((s >> 8)&3) == 3 for s in states], dtype=float); p = 1-q
    v = np.zeros(n); v[states.index(63)] = 1
    counter = e.Counter(M); r = cert['response']; eta = float(F(r['residual_penalty_Ha']))
    t = time.monotonic(); out = e.retained(counter.H, p, q, v, r, counter, True)-eta*counter.project(p, v)
    wall = time.monotonic()-t
    pi = np.flatnonzero(p).tolist(); qi = np.flatnonzero(q).tolist()
    K, _ = e.dense_retained(M, pi, qi, r); K -= eta*np.eye(len(pi))
    oracle = np.zeros(n); oracle[pi] = K@v[pi]
    error = float(np.max(np.abs(out-oracle)))
    if error > 1e-9 or counter.H_rhs != proof['H_actions_per_K_vector']:
        raise ValueError('Tuned physical execution mismatch')
    result = {'H_rhs': counter.H_rhs, 'seconds': wall, 'max_error_vs_spectral_oracle': error,
              'projector_calls': counter.projector_calls, 'recurrence_vector_entries': counter.recurrence_vector_entries,
              'vector_dimension': n, 'dense_H_entries': n*n, 'preparation': preparation,
              'certificate': 'h6_terminal_tuned.json', 'target_Ha': cert['target_Ha'],
              'eta_Ha': str(F(r['residual_penalty_Ha'])), 'eta_mHa': float(1000*F(r['residual_penalty_Ha'])),
              'wall_seconds': time.monotonic()-start,
              'scope': 'Uses a larger response allowance justified by the inherited terminal margin. Same complete 1 mHa target, different uniform response-error budget. Full determinant-block vectors remain diagnostic.'}
    (g.OUT/'tuned_execution.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('H_rhs', 'seconds', 'eta_mHa', 'max_error_vs_spectral_oracle')}))


if __name__ == '__main__': run()
