"""Exact signed-likelihood-moment planner; no posterior-search dependency."""
from fractions import Fraction as F
from experiments.v4_planner import _check


def signed_gram_plan(prior, likelihood, signs):
    m, actions = _check(prior, likelihood, signs)
    weights = [prior[i] * signs[i] for i in range(m)]
    bias = sum(weights)
    weighted = [[weights[i] * likelihood[i][a] for a in range(actions)] for i in range(m)]
    first = [sum(weighted[i][a] for i in range(m)) for a in range(actions)]
    gram = [[F(0) for _ in range(actions)] for _ in range(actions)]
    for a in range(actions):
        for b in range(a, actions):
            gram[a][b] = gram[b][a] = sum(weighted[i][a] * likelihood[i][b] for i in range(m))
    plan = plan_from_moments(bias, tuple(first), tuple(tuple(row) for row in gram))
    bit_count = max(max(x.numerator.bit_length(), x.denominator.bit_length())
                    for x in [bias] + first + [v for row in gram for v in row])
    plan['work'] = {'signed_weight_products': m, 'wL_products': m * actions,
                    'G_products': m * actions * (actions + 1) // 2,
                    'max_rational_bits': bit_count}
    return plan


def plan_from_moments(bias, first, gram):
    """Optimize a signed-score surrogate, even for inconsistent estimates.

    Empirical moments need not be realizable by one joint distribution. The
    output then has a surrogate risk, not a probability without an error bound.
    """
    actions = len(first)
    if not 1 <= actions <= 12 or len(gram) != actions or any(len(row) != actions for row in gram):
        raise ValueError('bounded moment dimensions required')
    values_to_check = [bias] + list(first) + [v for row in gram for v in row]
    if any(type(v) is not F or abs(v) > 1 for v in values_to_check):
        raise ValueError('exact bounded signed moments required')
    if any(gram[a][b] != gram[b][a] for a in range(actions) for b in range(actions)):
        raise ValueError('symmetric moment matrix required')
    values, children = [], []
    for a in range(actions):
        plus = [abs(gram[a][b]) + abs(first[a] - gram[a][b]) for b in range(actions)]
        minus = [abs(first[b] - gram[a][b]) + abs(bias - first[a] - first[b] + gram[a][b]) for b in range(actions)]
        b_plus = min(range(actions), key=lambda b: (-plus[b], b))
        b_minus = min(range(actions), key=lambda b: (-minus[b], b))
        values.append(F(1, 2) - (plus[b_plus] + minus[b_minus]) / 2)
        children.append((b_minus, b_plus))
    action = min(range(actions), key=lambda a: (values[a], a))
    return {'values': tuple(values), 'first_action': action, 'child_actions': tuple(children),
            'moments': {'bias': bias, 'first': tuple(first), 'gram': tuple(tuple(row) for row in gram)}}
