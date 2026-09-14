"""Exact joint row bounds on one-electron/one-hole occupation simplices.

Coordinates use the vertex basis of Q[x_i]/(x_i*x_j-delta_ij*x_i,
 sum(x_i)-1). The k entries explicitly evaluate k scalar endpoints;
 this algebra does not by itself reduce the number of physical assignments.
"""
from fractions import Fraction as F


def chart_modes(oracle, mask, bits):
    """Return (unresolved same-spin modes, one_hole), or None."""
    free = tuple(i for i in range(oracle.modes) if not mask & (1 << i))
    if len(free) < 2 or len({i % 2 for i in free}) != 1:
        return None
    need = oracle.target-(bits & oracle.spin_masks[free[0] % 2]).bit_count()
    if need not in (1, len(free)-1):
        return None
    return free, need != 1


def compile_simplex(oracle, mask, bits):
    """Return exact Q-row endpoint coordinates and local operation counts.

    No determinant integers or determinant actions are constructed. Each group
    and delta kernel is processed jointly over the k coordinate entries.
    """
    chart = chart_modes(oracle, mask, bits)
    if chart is None:
        return None
    free, one_hole = chart
    if oracle.close(mask, bits) != (mask, bits) or oracle.counts(mask, bits)[1] == 0:
        raise ValueError('Closed nonempty physical simplex required')
    k = len(free)
    zero, one = (F(0),)*k, (F(1),)*k
    def add(a, b): return tuple(x+y for x, y in zip(a, b))
    def multiply(a, b): return tuple(x*y for x, y in zip(a, b))
    def scale(a, c): return tuple(c*x for x in a)
    def complement(a): return tuple(1-x for x in a)
    n = [(F(bool(bits & (1 << i))),)*k for i in range(oracle.modes)]
    for vertex, i in enumerate(free):
        unit = tuple(F(j == vertex) for j in range(k))
        n[i] = complement(unit) if one_hole else unit
    q = [add(add(n[2*i], n[2*i+1]), (-F(1),)*k) for i in range(oracle.sites)]
    def valence(charges):
        out = one
        for charge in charges:
            out = multiply(out, complement(multiply(charge, charge)))
        return out
    source_p = valence(q)
    row = zero
    for support, coefficient in oracle.oracle.diagonal.items():
        term = one
        for i in range(oracle.modes):
            if support & (1 << i): term = multiply(term, n[i])
        row = add(row, scale(term, coefficient))
    kernels = {}
    counts = {'transition_groups_used': 0, 'metric_scalar_endpoint_evaluations': 0,
              'amplitude_scalar_endpoint_evaluations': 0, 'metric_factors_touched': 0}
    for (c, a), amplitude in oracle.groups.items():
        event = one
        for i in range(oracle.modes):
            if c & (1 << i): event = multiply(event, complement(n[i]))
            elif a & (1 << i): event = multiply(event, n[i])
        if event == zero: continue
        amp = zero
        for support, coefficient in amplitude.items():
            amp = add(amp, scale(n[support.bit_length()-1] if support else one, coefficient))
        if amp == zero: continue
        delta = tuple(((c >> (2*i)) & 3).bit_count()-((a >> (2*i)) & 3).bit_count()
                      for i in range(oracle.sites))
        if delta not in kernels:
            changed = [add(q[i], (F(delta[i]),)*k) for i in range(oracle.sites)]
            ratio = one
            for label, factor in oracle.local_factors:
                if not any(delta[i] for i, _ in label): continue
                before = after = one
                for i, power in label:
                    for _ in range(power):
                        before = multiply(before, q[i]); after = multiply(after, changed[i])
                exponents = tuple(y-x for x, y in zip(before, after))
                if any(power.denominator != 1 for power in exponents):
                    raise ValueError('Integer occupation exponents required')
                ratio = multiply(ratio, tuple(factor**int(power) for power in exponents))
                counts['metric_factors_touched'] += 1
                counts['metric_scalar_endpoint_evaluations'] += k
            kernels[delta] = multiply(complement(valence(changed)), ratio)
        penalty = multiply(multiply(event, tuple(abs(x) for x in amp)), kernels[delta])
        row = add(row, scale(penalty, -1))
        counts['transition_groups_used'] += 1
        counts['amplitude_scalar_endpoint_evaluations'] += k
    # R_Q <= H_ss <= sum|h_word|, so this artificial P value preserves min_Q.
    bound = sum(map(abs, oracle.oracle.h.values()), F(0))
    row = add(multiply(complement(source_p), row), scale(source_p, bound))
    counts.update(chart_modes=k, one_hole=one_hole, feasible_assignments=k,
                  ionic_assignments=oracle.counts(mask, bits)[1], metric_delta_kernels=len(kernels))
    return row, counts
