"""Memoized schedule wrapper for two nested polynomial Schur actions.

The schedule is deliberately an execution aid, not a new bound.  It preserves
the recursive operator order and memoizes only identical vector/action pairs.
"""
from fractions import Fraction as F
from research.composable_response_20260913 import recursive


def key(v):
    return tuple(v)


class MemoAction:
    """Count and memoize exact oracle applications by (layer, vector)."""
    def __init__(self, action, layer="H"):
        self.action = action
        self.layer = layer
        self.cache = {}
        self.calls = 0
        self.h_rhs = 0

    def __call__(self, v):
        k = (self.layer, key(v))
        if k not in self.cache:
            self.cache[k] = list(self.action(list(v)))
            self.calls += 1
            self.h_rhs += 1
        return list(self.cache[k])


def nested_action(H, P1, Q1, P2, Q2, v, first, second):
    """Reference nested schedule, with no cross-call memoization."""
    K1 = lambda x: recursive.retained_action(H, P1, Q1, x, first)
    return recursive.retained_action(K1, P2, Q2, v, second)


def flattened_action(H, P1, Q1, P2, Q2, v, first, second):
    """Equivalent schedule with exact common-subexpression elimination.

    This is a CSE wrapper, not an independent nonrecursive implementation.
    Every K1 evaluation is still expanded in the same order.  Only repeated
    identical H-right-hand sides are shared; projections and vector algebra
    remain explicit and are counted separately by the caller.
    """
    memo = MemoAction(H)
    def K1(x):
        return recursive.retained_action(memo, P1, Q1, x, first)
    out = recursive.retained_action(K1, P2, Q2, v, second)
    return out, {'unique_H_rhs': memo.calls, 'H_rhs': memo.h_rhs,
                 'cache_entries': len(memo.cache)}


def schedule_counts(k1, k2):
    """Formal counts before CSE, matching the existing oracle recurrence."""
    first = 2*k1 + 1
    return {'k1': k1, 'k2': k2, 'per_K1': first,
            'outer_rhs': (k2-1)*first,
            'full_K2': (2*k2+1)*first}


def factored_response_constant(order, z0):
    """Return c=T_{2k}(z0)/(T_{2k}(z0)+1) for F_k=c*p_{2k}.

    The identity follows from T_k(z)^2=(T_{2k}(z)+1)/2 and the definitions
    p_k=(1-r_k)/d and F_k=2p_k-d p_k^2.  This helper is algebraic only; it
    does not change the H-action schedule.
    """
    from fractions import Fraction
    t0 = Fraction(1); t1 = z0
    for _ in range(2*order):
        t0, t1 = t1, 2*z0*t1-t0
    return t0/(t0+1)
