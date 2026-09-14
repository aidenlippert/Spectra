"""Bounded fraction-free reduction for exact candidate-basis completion."""


def reduce_basis(columns, rhs, additions):
    """Return D, transformed rhs/additions with selected columns D*[I;0].

    Rows may be permuted. All inputs and every division are exact integers.
    The omitted earlier diagonal columns have the common final value D.
    """
    m, n = len(rhs), len(columns)
    if not 1 <= n <= m <= 85 or len(additions) > 85:
        raise ValueError('Bounded basis dimensions required')
    if any(len(c) != m for c in columns + additions):
        raise ValueError('Column dimension mismatch')
    data = columns + [rhs] + additions
    if any(type(v) is not int or abs(v).bit_length() > 4096 for c in data for v in c):
        raise ValueError('Bounded integer columns required')
    a = [list(row) for row in zip(*data)]
    previous = 1
    width = len(a[0])
    for k in range(n):
        pivot_row = next((i for i in range(k, m) if a[i][k]), None)
        if pivot_row is None:
            raise ValueError('Dependent selected columns')
        a[k], a[pivot_row] = a[pivot_row], a[k]
        pivot = a[k][k]
        if abs(pivot).bit_length() > 50000:
            raise ValueError('Intermediate integer budget exceeded')
        for i in range(m):
            if i == k:
                continue
            factor = a[i][k]
            for j in range(k + 1, width):
                numerator = pivot * a[i][j] - factor * a[k][j]
                value, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError('Nonexact fraction-free division')
                a[i][j] = value
            a[i][k] = 0
        previous = pivot
    return previous, [row[n] for row in a], [row[n + 1:] for row in a]
