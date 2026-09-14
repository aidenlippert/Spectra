"""Request-only prediction guarantees for the unconstrained affine raw ARX class.

The gate consumes training rows and a proposed input sequence.  Future outputs
are deliberately not an argument: an ambiguous query returns a witness.
"""
from fractions import Fraction as F
from .v7_identifiability import ambiguity_from_witness, verify_compatible_models, nullspace


def _q(x):
    return x if isinstance(x, F) else F(str(x))


def _feature(y_history, u_history, next_u, order=8, inputs=2):
    if len(y_history) < order or len(u_history) < order - 1:
        raise ValueError("insufficient ARX history")
    if len(next_u) != inputs or any(len(row) != inputs for row in u_history[-(order-1):]):
        raise ValueError("input width mismatch")
    # Exact v6_headroom.features ordering: y[t],...,y[t-p+1], then
    # u[t],...,u[t-p+1] for each input channel, then intercept.
    ys = [_q(x) for x in y_history[-order:]][::-1]
    past = u_history[-(order-1):][::-1] if order > 1 else []
    us = [list(map(_q, next_u))] + [list(map(_q, row)) for row in past]
    return ys + [us[i][j] for i in range(order) for j in range(inputs)] + [F(1)]


def _training(data):
    if isinstance(data, dict):
        X, y = data["X"], data["y"]
        theta0 = data.get("theta0")
        eta = data.get("eta")
        order = data.get("order", 8)
        inputs = data.get("inputs", 2)
    else:
        X, y = data
        theta0 = eta = None; order, inputs = 8, 2
    return X, y, theta0, eta, order, inputs


def request_prediction(training_data, y_history, u_history, next_u, tolerance, uncertainty,
                       *, theta0=None, witness=None, order=8, input_width=2, horizon=1):
    """Return a certified prediction or abstain with compatible-model witness.

    Only ``training_data``, observed ``history`` and requested inputs cross the
    API boundary.  ``uncertainty`` is an l-infinity training residual enclosure.
    """
    X, y, supplied_theta, supplied_eta, data_order, data_inputs = _training(training_data)
    order = data_order if isinstance(training_data, dict) else order
    input_width = data_inputs if isinstance(training_data, dict) else input_width
    theta0 = theta0 if theta0 is not None else supplied_theta
    eta = _q(uncertainty if uncertainty is not None else supplied_eta)
    if theta0 is None:
        raise ValueError("theta0 is required for a finite noisy guarantee")
    if type(horizon) is not int or horizon < 1:
        raise ValueError('positive integer horizon required')
    if _q(tolerance) <= 0 or eta < 0:
        raise ValueError("positive tolerance and nonnegative uncertainty required")
    if len(y_history) < order or len(u_history) < order - 1:
        return {"status": "no_guarantee", "reason": "insufficient history", "witness": None}
    r = _feature(y_history, u_history, next_u, order, input_width)
    if len(r) != len(X[0]):
        raise ValueError("training design width does not match ARX order/input width")
    queries = [r]
    active_witness = witness if witness is not None and sum(_q(x)*_q(z) for x,z in zip(r,witness)) else None
    basis = [] if active_witness is not None else nullspace(X)
    for step, r in enumerate(queries):
        v = active_witness or next((z for z in basis if sum(r[j]*z[j] for j in range(len(r)))), None)
        if v is not None:
            c = ambiguity_from_witness(X, y, r, v, theta0, eta, margin=4*_q(tolerance))
            cdict = c.__dict__.copy()
            cdict["step"] = step
            cdict["query"] = [str(x) for x in r]
            if not verify_compatible_models(X, y, r, c.theta_plus, c.theta_minus, eta, 4*_q(tolerance)):
                raise AssertionError("invalid ambiguity witness")
            return {"status": "no_guarantee", "reason": "requested query is ambiguous",
                    "witness": cdict, "horizon": horizon}
    if horizon != 1:
        return {"status": "no_guarantee", "reason": "multi-step interval propagation unavailable", "witness": None}
    # A row-space enclosure is centered at w^T y, not at an arbitrary feasible theta.
    from .v7_identifiability import identify_query, verify_row_certificate
    cert = identify_query(X, r, y=y, theta0=theta0, eta=eta, margin=_q(tolerance))
    if not cert.identified or cert.weights is None or not verify_row_certificate(X, r, cert.weights):
        return {"status": "no_guarantee", "reason": "query not certified", "witness": None}
    radius = eta * sum(abs(_q(w)) for w in cert.weights)
    center = sum(_q(w) * _q(yy) for w, yy in zip(cert.weights, y))
    if radius > _q(tolerance):
        return {"status": "no_guarantee", "reason": "support interval exceeds tolerance",
                "center": str(center), "radius": str(radius), "witness": None}
    return {"status": "guaranteed", "predictions": [str(center)],
            "support_error_bound": str(radius), "weights": [str(w) for w in cert.weights], "horizon": 1,
            "scope": "unconstrained affine raw ARX; training residual enclosure only"}


request_gate = request_prediction
