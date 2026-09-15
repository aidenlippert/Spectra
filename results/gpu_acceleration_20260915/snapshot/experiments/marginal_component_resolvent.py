"""Exact second-order resolvent bounds using small reference-block polynomials."""
from fractions import Fraction as F

from experiments.marginal_general_schur import moment_recurrence, resolvent_self_energy


def multiply(a, b):
    result = [[F(0)] * len(b[0]) for _ in a]
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            if x:
                for j, y in enumerate(b[k]):
                    if y:
                        result[i][j] += x * y
    return result


def transpose(a):
    return [list(row) for row in zip(*a)]


def full_recurrence(group):
    """Annihilate the entire block, including directions reached through R."""
    if 'full_recurrence' not in group:
        a = group['matrix']
        def apply(columns):
            return [{i: value for i, row in enumerate(a)
                     if (value := sum(x * column.get(j, F(0)) for j, x in enumerate(row)))}
                    for column in columns]
        group['full_recurrence'] = moment_recurrence(
            [{i: F(1)} for i in range(len(a))], apply, max_degree=12)
    return group['full_recurrence']


def polynomial_inverse(group, shift):
    a = group['matrix']
    recurrence = full_recurrence(group)
    inverse = resolvent_self_energy(recurrence['annihilator'], recurrence['moments'], shift, size=len(a))
    shifted = [[x - shift * (i == j) for j, x in enumerate(row)] for i, row in enumerate(a)]
    product = multiply(shifted, inverse)
    if product != [[F(i == j) for j in range(len(a))] for i in range(len(a))]:
        raise ValueError('Block polynomial failed exact inverse identity')
    return inverse


def block_action(data, shift, rhs):
    result, offset = [], 0
    for group in data['groups']:
        size = len(group['matrix'])
        result.extend(multiply(polynomial_inverse(group, shift), rhs[offset:offset + size]))
        offset += size
    return result


def second_order_correction(data, lower):
    """W^T A^-1 W <= W^T G W - X^T R X + (RX)^T T(RX).

    A=C+R-bI, G=(C-bI)^-1, T=(C-(b+eta)I)^-1, X=GW.
    The caller must certify every C block > (b+eta)I first.
    """
    eta = data['off_block_norm']
    if any(group['lower'] <= lower + eta for group in data['groups']):
        raise ValueError('Second-order response requires strict shifted block gaps')
    w = data['q_coupling']
    x = block_action(data, lower, w)
    y = multiply(data['off_block_matrix'], x)
    leading = multiply(transpose(w), x)
    linear = multiply(transpose(x), y)
    remainder = multiply(transpose(y), block_action(data, lower + eta, y))
    return [[a - b + c for a, b, c in zip(row, linear[i], remainder[i])]
            for i, row in enumerate(leading)]
