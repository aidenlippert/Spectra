"""Exact two-shot Bayesian design and fixed-policy risk certificates.

Risk certificates establish the risk of specified actions and leaf decisions.
They do not prove first-action optimality. The full teacher enumerates every
first action separately and is charged separately.
"""
from fractions import Fraction as F


def new_work():
    return {'likelihood_weight_products': 0, 'posterior_normalizations': 0,
            'posterior_updates': 0, 'one_step_action_scores': 0,
            'two_step_candidates': 0, 'max_rational_bits': 0}


def _add(work, key, count=1):
    if work is not None: work[key] = work.get(key, 0) + count


def _bits(work, values):
    if work is not None:
        work['max_rational_bits'] = max(work.get('max_rational_bits', 0),
                                       max((max(v.numerator.bit_length(), v.denominator.bit_length()) for v in values), default=0))


def _check(prior, table, signs):
    m = len(prior)
    if not 1 <= m <= 8 or len(table) != m or len(signs) != m:
        raise ValueError('bounded model dimensions required')
    if any(type(s) is not int or s not in (-1, 1) for s in signs):
        raise ValueError('strict integer signs required')
    if any(type(p) is not F or p < 0 for p in prior) or sum(prior) != 1:
        raise ValueError('normalized exact prior required')
    a = len(table[0])
    if not 1 <= a <= 12 or any(len(row) != a or any(type(p) is not F or not 0 <= p <= 1 for p in row) for row in table):
        raise ValueError('bounded exact likelihood table required')
    return m, a


def branch(prior, table, action, outcome, counter=None):
    if type(action) is not int or not 0 <= action < len(table[0]) or type(outcome) is not int or outcome not in (0, 1):
        raise ValueError('invalid action/outcome')
    joint = tuple(p * (table[i][action] if outcome else 1 - table[i][action]) for i, p in enumerate(prior))
    _add(counter, 'likelihood_weight_products', len(prior)); _add(counter, 'posterior_updates')
    probability = sum(joint)
    _bits(counter, joint)
    if probability == 0: return F(0), tuple(F(0) for _ in prior)
    post = tuple(p / probability for p in joint)
    _add(counter, 'posterior_normalizations', len(prior)); _bits(counter, post)
    return probability, post


def terminal_risk(prior, signs):
    return min(sum(p for p, s in zip(prior, signs) if s == sign) for sign in (-1, 1))


def _h1(prior, table, signs, work):
    values = []
    for a in range(len(table[0])):
        risk = F(0)
        for y in (0, 1):
            joint = tuple(p * (table[i][a] if y else 1 - table[i][a]) for i, p in enumerate(prior))
            _add(work, 'likelihood_weight_products', len(prior)); _bits(work, joint)
            risk += terminal_risk(joint, signs)
        values.append(risk)
        _add(work, 'one_step_action_scores')
    return tuple(values)


def h1_values(prior, table, signs, counter=None):
    _check(prior, table, signs)
    return _h1(prior, table, signs, counter)


def action_values_h2(prior, table, signs, counter=None):
    _, actions = _check(prior, table, signs)
    values = []
    for a in range(actions):
        risk = F(0)
        for y in (0, 1):
            probability, post = branch(prior, table, a, y, counter)
            if probability:
                risk += probability * min(_h1(post, table, signs, counter))
        values.append(risk)
        _add(counter, 'two_step_candidates')
    return tuple(values)


def eval_two_step_policy(prior, table, signs, first_action, second_policy='goal_greedy', counter=None):
    _, actions = _check(prior, table, signs)
    if type(first_action) is not int or not 0 <= first_action < actions or second_policy != 'goal_greedy':
        raise ValueError('invalid fixed policy')
    work = new_work() if counter is None else counter
    branches, risk = [], F(0)
    for y in (0, 1):
        probability, post = branch(prior, table, first_action, y, work)
        if not probability:
            branches.append({'outcome': y, 'probability': '0', 'second_action': None,
                             'second_values': [], 'leaf_predictions': []})
            continue
        values = _h1(post, table, signs, work)
        chosen = min(range(actions), key=lambda a: (values[a], a))
        predictions = []
        for y2 in (0, 1):
            mass = [p * (table[i][chosen] if y2 else 1 - table[i][chosen]) for i, p in enumerate(post)]
            _add(work, 'likelihood_weight_products', len(prior)); _bits(work, mass)
            plus = sum(p for p, s in zip(mass, signs) if s == 1)
            minus = sum(mass) - plus
            predictions.append(1 if plus >= minus else -1)
        risk += probability * values[chosen]
        branches.append({'outcome': y, 'probability': str(probability), 'second_action': chosen,
                         'second_values': [str(v) for v in values], 'leaf_predictions': predictions})
    _add(work, 'two_step_candidates')
    return {'first_action': first_action, 'branches': branches, 'risk': str(risk),
            'scope': 'specified two-shot actions and leaf decisions under admitted likelihoods'}


def verify_fixed_policy_certificate(cert, prior, table, signs, counter=None):
    """Direct joint-outcome integration; never calls the planner/generator."""
    try:
        _, actions = _check(prior, table, signs)
        first, branches = cert['first_action'], cert['branches']
        if type(first) is not int or not 0 <= first < actions or type(branches) is not list or len(branches) != 2:
            return False
        risk = F(0)
        for y, item in enumerate(branches):
            mass1 = tuple(p * (table[i][first] if y else 1 - table[i][first]) for i, p in enumerate(prior))
            _add(counter, 'likelihood_weight_products', len(prior)); _bits(counter, mass1)
            probability = sum(mass1)
            if type(item['outcome']) is not int or item['outcome'] != y or item['probability'] != str(probability): return False
            if not probability:
                if item['second_action'] is not None or item['leaf_predictions'] != [] or item['second_values'] != []: return False
                continue
            second, predictions = item['second_action'], item['leaf_predictions']
            if type(second) is not int or not 0 <= second < actions or len(predictions) != 2 or any(type(s) is not int or s not in (-1, 1) for s in predictions): return False
            for y2 in (0, 1):
                mass2 = tuple(p * (table[i][second] if y2 else 1 - table[i][second]) for i, p in enumerate(mass1))
                _add(counter, 'likelihood_weight_products', len(prior)); _bits(counter, mass2)
                risk += sum(p for p, s in zip(mass2, signs) if s != predictions[y2])
        return cert['risk'] == str(risk)
    except (AttributeError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return False


def compile_fixed_actions(prior, table, signs, first_action, child_actions, counter=None):
    """Compile already selected actions without searching any alternatives.

    This lets a moment planner use its child choices directly. Leaf decisions
    and risk are obtained from joint masses, without posterior divisions.
    """
    _, actions = _check(prior, table, signs)
    if (type(first_action) is not int or not 0 <= first_action < actions
            or len(child_actions) != 2
            or any(type(a) is not int or not 0 <= a < actions for a in child_actions)):
        raise ValueError('invalid action tree')
    branches, risk = [], F(0)
    for y, chosen in enumerate(child_actions):
        joint = tuple(p * (table[i][first_action] if y else 1 - table[i][first_action])
                      for i, p in enumerate(prior))
        _add(counter, 'likelihood_weight_products', len(prior)); _bits(counter, joint)
        probability = sum(joint)
        if not probability:
            branches.append({'outcome': y, 'probability': '0', 'second_action': None,
                             'second_values': [], 'leaf_predictions': []})
            continue
        predictions = []
        for y2 in (0, 1):
            mass = tuple(p * (table[i][chosen] if y2 else 1 - table[i][chosen])
                         for i, p in enumerate(joint))
            _add(counter, 'likelihood_weight_products', len(prior)); _bits(counter, mass)
            plus = sum(p for p, s in zip(mass, signs) if s == 1)
            minus = sum(mass) - plus
            predictions.append(1 if plus >= minus else -1)
            risk += min(plus, minus)
        branches.append({'outcome': y, 'probability': str(probability),
                         'second_action': chosen, 'second_values': [],
                         'leaf_predictions': predictions})
    return {'first_action': first_action, 'branches': branches, 'risk': str(risk),
            'scope': 'specified two-shot actions and leaf decisions under admitted likelihoods'}


def verify_policy_certificate(cert, prior, table, signs, counter=None):
    """Stronger audit: fixed-policy risk plus optimal second-shot choices."""
    if not verify_fixed_policy_certificate(cert, prior, table, signs, counter): return False
    try:
        for y, item in enumerate(cert['branches']):
            probability, post = branch(prior, table, cert['first_action'], y, counter)
            if not probability: continue
            values = _h1(post, table, signs, counter)
            if item['second_values'] != [str(v) for v in values] or values[item['second_action']] != min(values): return False
        return True
    except (KeyError, TypeError, ValueError, IndexError):
        return False
