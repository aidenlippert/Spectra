"""Exact joint row compilation in Q[x]/(x*x-x) on a two-mode chart."""
from fractions import Fraction as F


def add(*polys):
    return tuple(sum((p[i] for p in polys), F(0)) for i in (0, 1))


def scale(poly, value):
    return poly[0]*value, poly[1]*value


def multiply(a, b):
    return a[0]*b[0], a[0]*b[1]+a[1]*b[0]+a[1]*b[1]


ONE = (F(1), F(0))
ZERO = (F(0), F(0))


def complement(poly):
    return add(ONE, scale(poly, -1))


def absolute(poly):
    first = abs(poly[0])
    return first, abs(poly[0]+poly[1])-first


def exponential(base, exponent):
    a, b = map(F, exponent)
    if base <= 0 or a.denominator != 1 or b.denominator != 1:
        raise ValueError('Positive rational base and integer Boolean exponent required')
    first = base**int(a)
    return first, base**int(a+b)-first


def compile_line(oracle, mask, bits):
    """Return exact Q-row polynomial and honest local operation counts.

    None means this branch is not the supported two-mode, one-electron chart.
    The caller supplies a closed physical branch. No completed determinant
    state, determinant action, or whole-row endpoint lookup is constructed.
    Scalar abs/exp operations do evaluate two local polynomial endpoints.
    """
    free = [i for i in range(oracle.modes) if not mask & (1 << i)]
    if len(free) != 2 or free[0] % 2 != free[1] % 2:
        return None
    spin = oracle.spin_masks[free[0] % 2]
    if oracle.target-(bits & spin).bit_count() != 1:
        return None
    if oracle.close(mask, bits) != (mask, bits) or oracle.counts(mask, bits)[1] == 0:
        raise ValueError('Closed nonempty physical ionic chart required')
    n = [(F(bool(bits & (1 << i))), F(0)) for i in range(oracle.modes)]
    n[free[0]], n[free[1]] = (F(0), F(1)), (F(1), F(-1))
    q = [add(n[2*i], n[2*i+1], (F(-1), F(0))) for i in range(oracle.sites)]
    source_p = ONE
    for charge in q:
        source_p = multiply(source_p, complement(multiply(charge, charge)))
    source_q = complement(source_p)
    diagonal = ZERO
    for support, coefficient in oracle.oracle.diagonal.items():
        term = ONE
        for i in range(oracle.modes):
            if support & (1 << i):
                term = multiply(term, n[i])
        diagonal = add(diagonal, scale(term, coefficient))
    row = diagonal; kernels = {}
    counts = {'transition_groups_used': 0, 'metric_scalar_endpoint_evaluations': 0,
              'amplitude_scalar_endpoint_evaluations': 0, 'metric_factors_touched': 0}
    for (c, a), amplitude in oracle.groups.items():
        event = ONE
        for i in range(oracle.modes):
            bit = 1 << i
            if c & bit:
                event = multiply(event, complement(n[i]))
            elif a & bit:
                event = multiply(event, n[i])
        if event == ZERO:
            continue
        amp = ZERO
        for support, coefficient in amplitude.items():
            amp = add(amp, scale(n[support.bit_length()-1] if support else ONE, coefficient))
        if amp == ZERO:
            continue
        delta = tuple(((c >> (2*i)) & 3).bit_count()-((a >> (2*i)) & 3).bit_count()
                      for i in range(oracle.sites))
        if delta not in kernels:
            changed = [add(q[i], (F(delta[i]), F(0))) for i in range(oracle.sites)]
            target_p = ONE
            for charge in changed:
                target_p = multiply(target_p, complement(multiply(charge, charge)))
            ratio = ONE
            for label, factor in oracle.local_factors:
                if not any(delta[i] for i, _ in label):
                    continue
                before = after = ONE
                for i, power in label:
                    for _ in range(power):
                        before = multiply(before, q[i]); after = multiply(after, changed[i])
                ratio = multiply(ratio, exponential(factor, add(after, scale(before, -1))))
                counts['metric_factors_touched'] += 1
                counts['metric_scalar_endpoint_evaluations'] += 2
            kernels[delta] = multiply(complement(target_p), ratio)
        penalty = multiply(multiply(event, absolute(amp)), kernels[delta])
        row = add(row, scale(penalty, -1))
        counts['transition_groups_used'] += 1
        counts['amplitude_scalar_endpoint_evaluations'] += 2
    # H_ss <= sum|h_word|, and a row lower is no larger than H_ss. This
    # artificial valence value cannot lower the minimum over nonempty Q.
    exclusion_value = sum(map(abs, oracle.oracle.h.values()), F(0))
    row = add(multiply(source_q, row), scale(source_p, exclusion_value))
    counts.update(variable_mode=free[0], eliminated_mode=free[1], feasible_assignments=2,
                  ionic_assignments=oracle.counts(mask, bits)[1], metric_delta_kernels=len(kernels))
    return row, counts
