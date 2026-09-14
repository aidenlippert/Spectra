"""Exact spectral inputs for the conditional residual-targeting theorem.

This checks an excited-state threshold, not merely a complement threshold.
The theorem does not bound retained dimension, atom rank, or rational cost.
"""
from fractions import Fraction as F

from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_defect_dicke import reference_complement
from experiments.marginal_general_schur import norm_bound
from experiments.marginal_sector_reference import bracket, jacobi_data
from experiments.marginal_symbolic import add, scale


def reference_excitation(pairs):
    c = reference_complement(pairs)
    diagonal, squared = jacobi_data(2 * pairs, F(1, 5), 0)
    pivots = []
    for i, value in enumerate(diagonal):
        pivot = value - c - (squared[i - 1] / pivots[-1] if i else 0)
        if not pivot:
            raise ValueError('Zero Jacobi pivot: threshold needs refinement')
        pivots.append(pivot)
    if sum(x < 0 for x in pivots) != 1:
        raise ValueError('Reference excited-state threshold not certified')
    # Q is bounded by c; inertia proves exactly one P eigenvalue below c.
    return c, pivots


def parameters(h, pairs, upper):
    upper = F(upper)
    c, pivots = reference_excitation(pairs)
    delta = add(h, scale(hopping_polynomial(2 * pairs, F(1, 5)), -1))
    eta = norm_bound(delta, 'hermitian_pairs')
    beta = c - eta
    if upper >= beta:
        raise ValueError('Upper witness does not lie below excited-state threshold')
    lower = F(bracket(2 * pairs, F(1, 5), 0)['lower']) - eta
    ceiling = F(pairs * (pairs - 1), 2) + F(pairs, 5) + eta
    if upper < lower:
        raise ValueError('Inconsistent supplied upper bound')
    width = ceiling - lower
    factor = 1 - (beta - upper) / (2 * width)
    return {'reference_excitation_lower': str(c), 'reference_inertia_pivots': [str(x) for x in pivots],
            'perturbation_norm': str(eta), 'excited_lower': str(beta),
            'ground_lower': str(lower), 'spectral_upper': str(ceiling),
            'witness_upper': str(upper), 'spectral_width_bound': str(width),
            'contraction_factor': str(factor), 'contraction_factor_float': float(factor),
            'scope': 'Exact Ritz residual targeting, conditional on supplied valid upper; approximate Ritz adds internal-residual and next-vector errors'}


def ceiling(value, denominator=10**24):
    value = F(value)
    return F(-(-value.numerator * denominator // value.denominator), denominator)


def internal_residual(data, coefficients):
    from experiments.marginal_implicit_certificate import rayleigh
    from experiments.marginal_enlarged_schur import solve_positive
    mu = rayleigh(data, coefficients)
    g, a = data['metric'], data['projected_h']
    norm = sum(x * y * g[i][j] for i, x in enumerate(coefficients) for j, y in enumerate(coefficients))
    rhs = [[sum(x * y for x, y in zip(row, coefficients))] for row in a]
    solved = solve_positive(g, rhs)
    residual = sum(x[0] * y[0] for x, y in zip(rhs, solved)) / norm - mu * mu
    if residual < 0:
        raise ValueError('Negative internal residual squared')
    return mu, residual


def context(h, pairs, upper, error_floor):
    error_floor = F(error_floor)
    if error_floor <= 0:
        raise ValueError('Positive convergence error floor required')
    result = parameters(h, pairs, upper)
    result['error_floor'] = str(error_floor)
    return result


def propose_retained_lower(data, coefficients, guarantee):
    from experiments.marginal_implicit_certificate import rayleigh
    budget = (1 - F(guarantee['contraction_factor'])) * F(guarantee['error_floor'])
    # Reserve half the total error budget for the next-vector solve. A simple
    # decimal lower keeps the exact LDL input smaller than the Rayleigh ratio.
    unit = F(1, 10**24)
    if budget < 8 * unit:
        raise ValueError('Requested convergence floor exceeds proposal precision budget')
    candidate = rayleigh(data, coefficients) - budget / 4
    return -ceiling(-candidate)


def check_step(data, coefficients, next_data, next_coefficients, retained_lower, guarantee):
    """Check the error terms; the caller must rebuild the actual residual extension."""
    from experiments.marginal_implicit_certificate import rayleigh, rational_text
    from experiments.marginal_schur_transfer import ldl_pivots
    mu, residual = internal_residual(data, coefficients)
    next_mu = rayleigh(next_data, next_coefficients)
    upper = F(guarantee['witness_upper'])
    if mu > upper or next_mu > upper:
        raise ValueError('Target Rayleigh value exceeds fixed convergence upper')
    lower = F(retained_lower)
    matrix = [[x - lower * y for x, y in zip(a, g)]
              for a, g in zip(next_data['projected_h'], next_data['metric'])]
    if next_mu < lower or ldl_pivots(matrix) is None:
        raise ValueError('Next retained-space lower bound failed exact positivity')
    width = F(guarantee['spectral_width_bound'])
    q = F(guarantee['contraction_factor'])
    floor = F(guarantee['error_floor'])
    # Store short outward rational bounds rather than enormous residual ratios.
    residual_upper = ceiling(residual)
    zeta_upper = ceiling(next_mu - lower)
    error = residual_upper / (2 * width) + zeta_upper
    if error > (1 - q) * floor:
        raise ValueError('Approximate Ritz errors exceed convergence budget')
    return {'current_upper': rational_text(mu), 'next_upper': rational_text(next_mu),
            'internal_residual_squared_upper': str(residual_upper),
            'next_ritz_error_upper': str(zeta_upper), 'retained_lower': rational_text(lower),
            'additive_error_upper': str(error), 'allowed_additive_error': str((1 - q) * floor),
            'current_dimension': data['retained_dimension'], 'next_dimension': next_data['retained_dimension']}
