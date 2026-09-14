"""Exact, finite-dimensional identifiability certificates for linear queries.

All proof arithmetic uses :class:`fractions.Fraction`; floating point is never
used as a tolerance-based proof.  This module intentionally handles small
matrices (at most 32 columns) and does not infer physical parameter constraints.
"""
from dataclasses import dataclass
from fractions import Fraction
from typing import Any


Q = Fraction


def _q(x):
    return x if isinstance(x, Fraction) else Fraction(str(x))


def _matrix(X):
    try: a = [[_q(x) for x in row] for row in X]
    except (TypeError, ValueError): raise ValueError("X must be a rectangular numeric matrix")
    if not a or not a[0]:
        raise ValueError("X must be a nonempty matrix")
    n = len(a[0])
    if n > 32 or len(a) > 4096 or any(len(row) != n for row in a):
        raise ValueError("X must be rectangular with at most 32 columns")
    return a


def _rref(a):
    a = [row[:] for row in a]; m = len(a); n = len(a[0]) if m else 0
    pivots = []; row = 0
    for col in range(n):
        pivot = next((i for i in range(row, m) if a[i][col]), None)
        if pivot is None: continue
        a[row], a[pivot] = a[pivot], a[row]
        z = a[row][col]; a[row] = [v/z for v in a[row]]
        for i in range(m):
            if i != row and a[i][col]:
                z = a[i][col]; a[i] = [u-z*v for u,v in zip(a[i], a[row])]
        pivots.append(col); row += 1
        if row == m: break
    return a, pivots


def _solve(A, b):
    aug = [row[:] + [val] for row, val in zip(A, b)]
    rr, piv = _rref(aug); n = len(A[0])
    if any(all(not x for x in row[:n]) and row[n] for row in rr):
        return None
    out = [Q(0) for _ in range(n)]
    for i, p in enumerate(piv):
        if p < n: out[p] = rr[i][n]
    return out


def nullspace(X):
    """Return exact basis vectors for ``ker(X)``."""
    a = _matrix(X); rr, piv = _rref(a); free = [j for j in range(len(a[0])) if j not in piv]
    basis = []
    for f in free:
        v = [Q(0)] * len(a[0]); v[f] = Q(1)
        for i,p in enumerate(piv): v[p] = -rr[i][f]
        basis.append(v)
    return basis


def verify_row_certificate(X, r, weights):
    """Independent exact verifier for ``X.T @ weights == r``."""
    try: a = _matrix(X); q = [_q(x) for x in r]; w = [_q(x) for x in weights]
    except (TypeError, ValueError, ZeroDivisionError): return False
    if len(q) != len(a[0]) or len(w) != len(a): return False
    return all(sum(a[i][j]*w[i] for i in range(len(a))) == q[j] for j in range(len(q)))


def verify_nullspace_witness(X, v, r=None):
    try:
        a = _matrix(X); z = [_q(x) for x in v]
        if r is not None and len(r) != len(a[0]): return False
        q = None if r is None else [_q(x) for x in r]
    except (TypeError, ValueError, ZeroDivisionError): return False
    if len(z) != len(a[0]) or not any(z): return False
    if any(sum(row[j]*z[j] for j in range(len(z))) for row in a): return False
    return q is None or sum(q[j]*z[j] for j in range(len(q))) != 0


@dataclass
class QueryCertificate:
    identified: bool
    weights: list | None = None
    nullspace_witness: list | None = None
    theta0: list | None = None
    theta_plus: list | None = None
    theta_minus: list | None = None
    error_bound: Fraction | None = None
    reason: str = ""


def identify_query(X, r, y=None, *, margin=Q(1), theta0=None, eta=None):
    """Certify exact or narrowly noisy query identifiability.

    ``identified=True`` denotes row-space membership; it implies uniqueness
    only for exact compatible data, not a positive noisy radius.
    In noisy mode ``theta0`` and compatible data must be supplied;
    the residual check is ``||X theta0-y||_inf <= eta``.  ``y`` is never used
    to manufacture a proof of exact identifiability.
    """
    a = _matrix(X); q = [_q(x) for x in r]
    if len(q) != len(a[0]): raise ValueError("query dimension mismatch")
    if _q(margin) <= 0: raise ValueError("margin must be positive")
    if eta is not None and _q(eta) < 0: raise ValueError("eta must be nonnegative")
    if y is not None and len(y) != len(a): raise ValueError("y length mismatch")
    if theta0 is not None and len(theta0) != len(a[0]): raise ValueError("theta0 length mismatch")
    if eta is not None and (theta0 is None or y is None): raise ValueError("noisy certificate requires theta0 and y")
    if y is not None:
        yy = [_q(x) for x in y]
        if eta is None and _solve(a, yy) is None: raise ValueError("incompatible data: y is not in exact column space")
        if theta0 is not None:
            residual = [sum(row[j]*_q(theta0[j]) for j in range(len(row))) - yy[i] for i,row in enumerate(a)]
            if max(map(abs, residual), default=Q(0)) > (_q(eta) if eta is not None else Q(0)):
                raise ValueError("incompatible data: supplied theta0 exceeds declared eta")
    weights = _solve([[a[i][j] for i in range(len(a))] for j in range(len(a[0]))], q)
    if weights is not None:
        bound = None
        if eta is not None: bound = _q(eta) * sum(abs(x) for x in weights)
        return QueryCertificate(True, weights=weights, error_bound=bound, reason="query lies in row space")
    witness = next((v for v in nullspace(a) if sum(q[j]*v[j] for j in range(len(q)))), None)
    if witness is None: raise ArithmeticError("row-space failure without nullspace witness")
    if eta is not None:
        if theta0 is None or y is None: raise ValueError("noisy certificate requires theta0 and y")
        t0 = [_q(x) for x in theta0]; yy = [_q(x) for x in y]
        residual = [sum(row[j]*t0[j] for j in range(len(t0))) - yy[i] for i,row in enumerate(a)]
        if max(map(abs, residual), default=Q(0)) > _q(eta):
            raise ValueError("incompatible data: supplied theta0 exceeds declared eta")
        delta = _q(margin) / (2 * abs(sum(q[j]*witness[j] for j in range(len(q)))))
        plus = [t0[j] + delta*witness[j] for j in range(len(t0))]
        minus = [t0[j] - delta*witness[j] for j in range(len(t0))]
        return QueryCertificate(False, nullspace_witness=witness, theta0=t0, theta_plus=plus, theta_minus=minus,
                                reason="exact nullspace gives compatible noisy models")
    if y is not None:
        yy = [_q(x) for x in y]; t0 = _solve(a, yy)
        if t0 is None: raise ValueError("incompatible data: y is not in the exact column space of X")
        delta = _q(margin) / (2 * abs(sum(q[j]*witness[j] for j in range(len(q)))))
        return QueryCertificate(False, nullspace_witness=witness, theta0=t0,
                                theta_plus=[t0[j]+delta*witness[j] for j in range(len(t0))],
                                theta_minus=[t0[j]-delta*witness[j] for j in range(len(t0))],
                                reason="compatible exact models disagree on query")
    return QueryCertificate(False, nullspace_witness=witness, reason="query outside row space")


def ambiguity_from_witness(X, y, r, v, theta0, eta, margin=Q(1)):
    """Fast path using a caller-supplied exact nullspace witness."""
    a = _matrix(X); q = [_q(x) for x in r]; z = [_q(x) for x in v]; t = [_q(x) for x in theta0]
    if len(y) != len(a) or len(q) != len(a[0]) or len(z) != len(a[0]) or len(t) != len(a[0]): raise ValueError("dimension mismatch")
    if _q(eta) < 0 or _q(margin) <= 0: raise ValueError("eta must be nonnegative and margin positive")
    if not verify_nullspace_witness(a, z, q): raise ValueError("invalid nullspace witness")
    yy = [_q(x) for x in y]; residual = [sum(row[j]*t[j] for j in range(len(t)))-yy[i] for i,row in enumerate(a)]
    if max(map(abs,residual), default=Q(0)) > _q(eta): raise ValueError("incompatible theta0")
    d = _q(margin)/(2*abs(sum(q[j]*z[j] for j in range(len(q)))))
    return QueryCertificate(False, nullspace_witness=z, theta0=t,
        theta_plus=[t[j]+d*z[j] for j in range(len(t))], theta_minus=[t[j]-d*z[j] for j in range(len(t))],
        reason="exact supplied nullspace gives compatible noisy models")


def verify_compatible_models(X, y, r, theta_plus, theta_minus, eta, margin):
    """Independent direct verifier for two feasible, disagreeing models."""
    try:
        a=_matrix(X); yy=[_q(x) for x in y]; q=[_q(x) for x in r]; p=[_q(x) for x in theta_plus]; m=[_q(x) for x in theta_minus]; e=_q(eta); g=_q(margin)
        if len(yy)!=len(a) or len(q)!=len(a[0]) or len(p)!=len(a[0]) or len(m)!=len(a[0]) or e<0 or g<=0:return False
        rp=[sum(row[j]*p[j] for j in range(len(p)))-yy[i] for i,row in enumerate(a)]
        rm=[sum(row[j]*m[j] for j in range(len(m)))-yy[i] for i,row in enumerate(a)]
        return max(map(abs,rp),default=Q(0))<=e and max(map(abs,rm),default=Q(0))<=e and abs(sum(q[j]*(p[j]-m[j]) for j in range(len(q))))>=g
    except (TypeError, ValueError, IndexError): return False
