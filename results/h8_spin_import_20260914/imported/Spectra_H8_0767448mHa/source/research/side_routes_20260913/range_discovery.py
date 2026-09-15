"""Numerical proposals for finite-range witnesses; exact verifier is separate."""
from fractions import Fraction as F
from math import lcm
import numpy as np


def site_features(model):
    m = model['modes']
    fields = list(map(F, model['fields']))
    scale = lcm(*(a.denominator for a in fields))
    values = [int(a * scale) for a in fields]
    center = sorted(values)[m // 2]
    features = [(i % 2, int(i in (0, m - 1)), values[i] - center) for i in range(m)]
    if max(abs(x) for row in features for x in row) > 32:
        raise ValueError('Field feature exceeds discovery runner budget')
    return features


def rational_amplitude(model, radius, parameters, denominator=65536):
    if radius not in range(4) or len(parameters) != radius + 3:
        raise ValueError('Invalid discovery parameters')
    weights = [F(max(1, int(round(float(np.exp(a)) * denominator))), denominator) for a in parameters]
    return amplitude_from_weights(model, radius, list(map(str, weights))), list(map(str, weights))


def amplitude_from_weights(model, radius, raw_weights):
    weights = list(map(F, raw_weights))
    if len(weights) != radius + 3 or any(w <= 0 for w in weights):
        raise ValueError('Invalid exact amplitude parameters')
    m = model['modes']
    sites = []
    for row in site_features(model):
        y = F(1)
        for weight, exponent in zip(weights[radius:], row):
            y *= weight ** exponent
        sites.append(str(y))
    amplitude = {'sites': sites, 'pairs': [[str(weights[d - 1])] * (m - d) for d in range(1, radius + 1)]}
    return amplitude


class NumericalChain:
    """Precompile transitions once, then evaluate proposal objectives in NumPy."""

    def __init__(self, model, radius):
        m, n = model['modes'], model['particles']
        if radius not in range(4) or radius >= m:
            raise ValueError('Invalid discovery radius')
        self.radius = radius
        self.dimension = radius + 3
        self.offset = float(F(model.get('offset', 0)))
        features = np.array(site_features(model), dtype=float)
        t, v, h = [np.array(list(map(F, model[k])), dtype=float) for k in ('hopping', 'interaction', 'fields')]
        memory = min(m - 1, 2 * radius + 1)
        self.layers = []
        keys = [(0, 0)]
        for i in range(m):
            sources, destinations, fulls, probabilities = [], [], [], []
            newkeys = {}
            for src, (count, tail) in enumerate(keys):
                for bit in (0, 1):
                    number = count + bit
                    if number > n or number + m - i - 1 < n:
                        continue
                    full = (tail << 1) | bit
                    key = (number, full & ((1 << memory) - 1))
                    dst = newkeys.setdefault(key, len(newkeys))
                    p = np.zeros(self.dimension)
                    p[radius:] = 2 * bit * features[i]
                    for d in range(1, radius + 1):
                        if i >= d:
                            p[d - 1] = 2 * bit * ((tail >> (d - 1)) & 1)
                    sources.append(src); destinations.append(dst); fulls.append(full); probabilities.append(p)
            fulls = np.array(fulls, dtype=np.int64)
            diagonal = h[i] * (fulls & 1)
            rate_features, rates = [], []
            for bond in range(m - 1):
                end = min(m - 1, bond + 1 + radius)
                if end != i:
                    continue
                bit = lambda j: (fulls >> (end - j)) & 1
                a, b = bit(bond), bit(bond + 1)
                diagonal = diagonal + v[bond] * a * b
                g = np.zeros((len(fulls), self.dimension))
                g[:, radius:] = (b - a)[:, None] * (features[bond] - features[bond + 1])
                after = lambda j: b if j == bond else a if j == bond + 1 else bit(j)
                for d in range(1, radius + 1):
                    for p in sorted({bond - d, bond, bond + 1 - d, bond + 1}):
                        q = p + d
                        if 0 <= p < q < m and {p, q} != {bond, bond + 1}:
                            g[:, d - 1] += after(p) * after(q) - bit(p) * bit(q)
                rate_features.append(g)
                rates.append(t[bond] * (a != b))
            self.layers.append((np.array(sources), np.array(destinations), np.array(probabilities),
                                diagonal, rate_features, rates, len(newkeys)))
            keys = list(newkeys)
        self.transitions = sum(len(a[0]) for a in self.layers)

    def interval(self, parameters):
        parameters = np.asarray(parameters)
        if parameters.shape != (self.dimension,) or not np.all(np.isfinite(parameters)):
            raise ValueError('Invalid numerical proposal')
        minimum = np.array([self.offset])
        Z, S = np.array([1.]), np.array([self.offset])
        for source, dest, p, diagonal, rate_features, rates, size in self.layers:
            local = diagonal.copy()
            for g, rate in zip(rate_features, rates):
                local -= rate * np.exp(g @ parameters)
            weight = np.exp(p @ parameters)
            nxt = np.full(size, np.inf)
            np.minimum.at(nxt, dest, minimum[source] + local)
            z = np.bincount(dest, weights=Z[source] * weight, minlength=size)
            s = np.bincount(dest, weights=(S[source] + Z[source] * local) * weight, minlength=size)
            scale = np.max(z)
            minimum, Z, S = nxt, z / scale, s / scale
        lower, upper = float(np.min(minimum)), float(np.sum(S) / np.sum(Z))
        if not np.isfinite(lower + upper) or upper < lower - 1e-8:
            raise ArithmeticError('Invalid numerical screening interval')
        return lower, upper
