"""Integer residual verification of approximate unit-LDL positivity factors."""
from fractions import Fraction as F


def rounded_division(numerator, denominator):
    if denominator <= 0:
        raise ValueError('Positive rounding denominator required')
    q, r = divmod(numerator, denominator)
    return q + int(2*r > denominator or (2*r == denominator and q % 2))


def scaled_matrix(a, scale):
    n = len(a)
    if (type(scale) is not int or not 1 <= scale <= 10**64 or not 1 <= n <= 256
            or any(len(row) != n for row in a)):
        raise ValueError('Invalid bounded factor matrix or integer scale')
    if any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Symmetric matrix required')
    values = [[F(x) * scale for x in row] for row in a]
    if any(x.denominator != 1 for row in values for x in row):
        raise ValueError('Factor scale must represent the input matrix exactly')
    return [[x.numerator for x in row] for row in values]


def propose(a, scale=10**24):
    values = scaled_matrix(a, scale)
    lower, diagonal = [], []
    for i, row in enumerate(values):
        current = []
        for j in range(i):
            numerator = row[j] * scale**2 - sum(current[k] * lower[j][k] * diagonal[k] for k in range(j))
            current.append(rounded_division(numerator, scale * diagonal[j]))
        pivot = rounded_division(row[i] * scale**2 - sum(x*x*d for x, d in zip(current, diagonal)), scale**2)
        if pivot <= 0:
            raise ValueError('Fixed-point proposal found a nonpositive pivot')
        lower.append(current)
        diagonal.append(pivot)
    return {'scale': scale, 'lower': lower, 'diagonal': diagonal}


def verify(a, factor):
    if type(factor) is not dict:
        raise ValueError('Missing fixed-point factor')
    scale = factor.get('scale')
    values = scaled_matrix(a, scale)
    n = len(values)
    lower, diagonal = factor.get('lower'), factor.get('diagonal')
    if (type(lower) is not list or len(lower) != n
            or any(type(row) is not list or len(row) != i or any(type(x) is not int for x in row) for i, row in enumerate(lower))
            or type(diagonal) is not list or len(diagonal) != n
            or any(type(x) is not int or x <= 0 for x in diagonal)):
        raise ValueError('Malformed unit-lower factor or nonpositive diagonal')
    rows = [row + [scale] for row in lower]
    errors = [0] * n
    for i in range(n):
        for j in range(i + 1):
            error = abs(values[i][j] * scale**2 - sum(rows[i][k] * diagonal[k] * rows[j][k] for k in range(j + 1)))
            errors[i] += error
            if i != j:
                errors[j] += error
    row_weights, column_weights = [], [0] * n
    for i in range(n):
        value = sum(abs(x) * row_weights[j] for j, x in enumerate(lower[i]))
        row_weights.append(scale + (value + scale - 1) // scale)
    for i in reversed(range(n)):
        value = sum(abs(lower[j][i]) * column_weights[j] for j in range(i + 1, n))
        column_weights[i] = scale + (value + scale - 1) // scale
    residual = F(max(errors), scale**3)
    inverse_product = F(max(row_weights) * max(column_weights), scale**2)
    margin = F(min(diagonal), scale) / inverse_product - residual
    if margin <= 0:
        raise ValueError('Factor residual exceeds the verified positive margin')
    return {'margin': margin, 'residual_norm_bound': residual,
            'inverse_norm_squared_bound': inverse_product, 'dimension': n}
