"""Non-enumerating exact SOS map tasks and a separate small ladder oracle."""
from fractions import Fraction as F
from itertools import combinations
import random

from experiments.marginal_hunt_car import mul, mono
from experiments.marginal_symbolic import word_product


def oracle(left, right):
    return dict(word_product(tuple(left), tuple(right)))


def dagger(word):
    return tuple((1 - c, i) for c, i in reversed(word))


def independent_action(word, bits):
    """Literal ladder application; used only as a small algebra oracle, never for solving H."""
    sign = 1
    for creation, mode in reversed(word):
        occupied = bool(bits & (1 << mode))
        if occupied == bool(creation):
            return None, 0
        if (bits & ((1 << mode) - 1)).bit_count() % 2:
            sign = -sign
        bits ^= 1 << mode
    return bits, sign


def ladder_identity(left, right, polynomial):
    modes = sorted({i for _, i in left + right} | {i for w in polynomial for _, i in w})
    if len(modes) > 6:
        raise ValueError("Small independent oracle envelope exceeded")
    relabel = {m: i for i, m in enumerate(modes)}
    word = tuple((c, relabel[i]) for c, i in left + right)
    p = {tuple((c, relabel[i]) for c, i in w): v for w, v in polynomial.items()}
    for bits in range(1 << len(modes)):
        end, value = independent_action(word, bits)
        expected = {end: value} if value else {}
        actual = {}
        for w, a in p.items():
            state, b = independent_action(w, bits)
            if b:
                actual[state] = actual.get(state, 0) + a * b
        if {s: a for s, a in actual.items() if a} != expected:
            return False
    return True


def validate_polynomial(p):
    if not isinstance(p, dict):
        raise ValueError("Polynomial must be a dictionary")
    for w, c in p.items():
        if not isinstance(w, tuple) or not isinstance(c, (int, F)) or not c or F(c).denominator != 1:
            raise ValueError("Expected exact nonzero integral CAR coefficients")
        if any(type(flag) is not int or flag not in (0, 1) or type(i) is not int or i < 0 for flag, i in w):
            raise ValueError("Invalid ladder operator")
        if tuple(sorted(w, key=lambda x: (1 - x[0], x[1]))) != w or len(set(w)) != len(w):
            raise ValueError("Noncanonical polynomial")


def product_cases(seed=0, count=400):
    rng = random.Random(seed)
    cases = [((), ()), (((0, 0),), ((1, 0),)), (((1, 5), (1, 5)), ())]
    for _ in range(count):
        modes = rng.sample(range(100), rng.randint(1, 6))
        cases.append(tuple(tuple((rng.randrange(2), rng.choice(modes)) for _ in range(rng.randrange(4)))
                           for side in range(2)))
    return cases


def task_words(spec):
    """The same pair-supported mixed dictionaries as Spectra, before floating spin reduction.

    Different support graphs model different locality/collective regimes. No occupation
    configurations, Hamiltonian diagonalization or global cubic coefficient map is constructed.
    """
    spatial = spec["spatial"]
    shape = spec["family"]
    if shape == "collective":
        supports = list(combinations(range(spatial), 2))
    elif shape == "chain":
        supports = [tuple(range(i, min(spatial, i + spec.get("width", 3))))
                    for i in range(spatial - spec.get("width", 3) + 1)]
    elif shape == "star":
        supports = [(0, i, j) for i, j in combinations(range(1, spatial), 2)]
    elif shape == "disjoint":
        supports = [tuple(range(i, min(i + 3, spatial))) for i in range(0, spatial, 3)]
    else:
        raise ValueError("Unknown declared support graph")
    words = set()
    for positions in supports:
        modes = [2 * p + spin for p in positions for spin in (0, 1)]
        words.update(((0, i),) for i in modes)
        for i, j in combinations(modes, 2):
            for k in modes:
                words.add(((1, k), (0, j), (0, i)))
    weight = spec.get("weight", -1)
    words = [w for w in words if sum((2 * c - 1) * (1 if i % 2 == 0 else -1) for c, i in w) == weight]
    words.sort(key=lambda w: (len(w), w))
    if spec.get("adjoint"):
        words = [dagger(w) for w in words]
    # Order is part of the input and may not be assumed canonical by a retained algorithm.
    random.Random(spec.get("permutation_seed", 0)).shuffle(words)
    offset = spec.get("mode_offset", 0)
    stride = spec.get("mode_stride", 1)
    return [tuple((c, offset + stride * i) for c, i in w) for w in words]


def gram_reference(words, normal=oracle):
    result = {}
    adjoints = [dagger(w) for w in words]
    for i in range(len(words)):
        for j in range(i, len(words)):
            p = normal(adjoints[i], words[j])
            if i != j:
                for w, c in normal(adjoints[j], words[i]).items():
                    p[w] = p.get(w, 0) + c
            for w, c in p.items():
                if c:
                    result.setdefault(w, {})[(i, j)] = int(c)
    return result


def validate_map(value, words):
    if not isinstance(value, dict):
        raise ValueError("Map must be dict canonical_word -> dict (i,j): integer")
    for word, row in value.items():
        validate_polynomial({word: 1})
        if not isinstance(row, dict) or not row:
            raise ValueError("Empty/malformed coefficient row")
        for pair, coefficient in row.items():
            if not isinstance(pair, tuple) or len(pair) != 2 or any(type(i) is not int for i in pair):
                raise ValueError("Invalid Gram coordinate")
            if not 0 <= pair[0] <= pair[1] < len(words) or type(coefficient) is not int or not coefficient:
                raise ValueError("Invalid exact map coefficient")


def map_cost(value):
    return {"coefficient_rows": len(value), "nonzero_entries": sum(map(len, value.values()))}


BASE_SOURCE = '''def propose(payload):
    words = payload["words"]
    adjoints = [tuple((1-c, i) for c, i in reversed(w)) for w in words]
    result = {}
    for i in range(len(words)):
        for j in range(i, len(words)):
            p = oracle(adjoints[i], words[j])
            if i != j:
                for w, c in oracle(adjoints[j], words[i]).items():
                    p[w] = p.get(w, 0) + c
            for w, c in p.items():
                if c:
                    result.setdefault(w, {})[(i,j)] = int(c)
    return result
'''
