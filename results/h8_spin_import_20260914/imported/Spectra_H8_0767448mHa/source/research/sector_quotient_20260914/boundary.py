"""Matrix-free boundary-point proposals for a product PSD cone plus free scalars.

min c.x subject to Bx=h. PSD blocks use the Frobenius inner product. The
boundary-point dual augmented-Lagrangian updates never materialize the dense
compact Gram-to-sextic map. This module accepts no energy claims.
"""
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg


def positive(A):
    ev, U = np.linalg.eigh((A+A.T)/2)
    return (U*np.maximum(ev, 0))@U.T


def solve(A, AT, free, rhs, start_Q, start_x, seconds, callback=None, mu=1., max_cg=80,
          diagonal=None, initial_y=None, max_outer=10000, tolerance=1e-8):
    """Inexact dual ALM, with every actual residual and CG failure observable."""
    started = time.monotonic()
    Q = [positive(q) for q in start_Q]
    x = np.array(start_x, copy=True)
    c = np.zeros(free.shape[1]); c[0] = -1.
    y = np.zeros(len(rhs)) if initial_y is None else np.array(initial_y, copy=True)
    Z = [np.zeros_like(q) for q in Q] if initial_y is None else [positive(-a) for a in AT(y)]
    def normal(v):
        return A(AT(v)) + free@(free.T@v)
    operator = LinearOperator((len(rhs), len(rhs)), matvec=normal)
    preconditioner = None if diagonal is None else LinearOperator(operator.shape, matvec=lambda v: v/np.maximum(diagonal, 1e-8))
    stats = []; last_print = started; cg_calls = 0
    for iteration in range(max_outer):
        elapsed = time.monotonic()-started
        if elapsed >= seconds:
            break
        residual = rhs-A(Q)-free@x
        target = free@c-A(Z)+residual/mu
        counter = [0]
        def count(_):
            counter[0] += 1
        # The reachable right side makes the semidefinite normal equation
        # consistent in exact arithmetic. CG status and true residual are
        # retained; failure never certifies infeasibility or convergence.
        y, status = cg(operator, target, x0=y, rtol=min(1e-6, max(1e-10, .01*np.linalg.norm(residual))),
                       atol=1e-11, maxiter=max_cg, M=preconditioner, callback=count)
        cg_calls += counter[0]
        At = AT(y)
        Znew = []; Qnew = []
        for q, a in zip(Q, At):
            W = -a-q/mu
            z = positive(W)
            Znew.append(z)
            Qnew.append((q+mu*(a+z)+q.T+mu*(a.T+z.T))/2)
        x += mu*(free.T@y-c)
        Q, Z = Qnew, Znew
        primal = rhs-A(Q)-free@x
        free_error = free.T@y-c
        psd_error = [a+z for a, z in zip(At, Z)]
        pn = float(np.linalg.norm(primal))
        dn = float(np.sqrt(sum(np.sum(s*s) for s in psd_error)+np.dot(free_error, free_error)))
        dual_objective = float(rhs@y)
        record = {'iteration': iteration+1, 'seconds': time.monotonic()-started,
                  'b': float(x[0]), 'primal_residual_l2': pn, 'dual_residual_l2': dn,
                  'dual_minimization_objective': dual_objective,
                  'objective_gap': float(-x[0]-dual_objective), 'mu': mu,
                  'cg_iterations': counter[0], 'cg_status': int(status), 'total_cg_iterations': cg_calls}
        if iteration == 0 or time.monotonic()-last_print >= 10 or pn < tolerance and dn < tolerance:
            stats.append(record)
            if callback:
                callback(record, Q, x, y, Z)
            last_print = time.monotonic()
        if pn < tolerance and dn < tolerance and abs(record['objective_gap']) < tolerance:
            break
        if (iteration+1) % 25 == 0:
            if pn > 5*dn:
                mu = min(mu*1.5, 1e4)
            elif dn > 5*pn:
                mu = max(mu/1.5, 1e-5)
    return Q, x, y, Z, stats
