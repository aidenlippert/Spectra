"""Exact finite-shot coupled Clifford experiment family.

The planner receives all candidate models and their likelihoods, never the
selected hidden model. target_sign is a prediction label, not a physical gate.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
import numpy as np
from experiments.pauli import multiply


def _word(word, n=None):
    if type(word) is not str or not 2 <= len(word) <= 4 or (n is not None and len(word) != n) or any(p not in 'IXYZ' for p in word) or all(p == 'I' for p in word):
        raise ValueError('nonidentity Pauli word on 2..4 qubits required')


@dataclass(frozen=True)
class Model:
    gates: tuple
    target_sign: int

    def __post_init__(self):
        if type(self.gates) is not tuple or not 1 <= len(self.gates) <= 2 or type(self.target_sign) is not int or self.target_sign not in (-1, 1):
            raise ValueError('one/two gates and a binary target sign required')
        n = len(self.gates[0][0])
        for gate in self.gates:
            if type(gate) is not tuple or len(gate) != 2: raise ValueError('invalid gate')
            p, sign = gate
            _word(p, n)
            if type(sign) is not int or sign not in (-1, 1) or sum(x != 'I' for x in p) > 2:
                raise ValueError('signed one/two-site gates required')


@dataclass(frozen=True)
class Experiment:
    prep: str
    measurement: str
    cost: int

    def __post_init__(self):
        _word(self.prep)
        _word(self.measurement, len(self.prep))
        if type(self.cost) is not int or self.cost < 1: raise ValueError('positive preparation/readout unit cost required')


def conj(q, p, sign):
    _word(q); _word(p, len(q))
    if type(sign) is not int or sign not in (-1, 1): raise ValueError('invalid gate sign')
    left, out = multiply(p, q)
    right, _ = multiply(q, p)
    if left == right: return q, 1
    phase = sign * 1j * left
    if phase.imag != 0 or phase.real not in (-1, 1): raise AssertionError('Clifford image must be signed Hermitian Pauli')
    return out, int(phase.real)


def transform(q, model):
    _word(q, len(model.gates[0][0]))
    word, coefficient = q, 1
    # gates are listed in chronological Schrödinger execution order.
    for p, sign in reversed(model.gates):
        word, factor = conj(word, p, sign)
        coefficient *= factor
    return word, coefficient


def expectation(model, experiment):
    word, coefficient = transform(experiment.measurement, model)
    return coefficient if word == experiment.prep else 0


def likelihood_table(models, actions, visibilities):
    if not 1 <= len(models) <= 8 or not 1 <= len(actions) <= 12 or len(visibilities) != len(actions):
        raise ValueError('bounded models/actions and one visibility per action required')
    if any(not isinstance(v, F) or not 0 <= v <= 1 for v in visibilities):
        raise ValueError('exact visibility in [0,1] required')
    return tuple(tuple((1 + visibilities[j] * expectation(model, action)) / 2
                       for j, action in enumerate(actions)) for model in models)


def sample(model, experiment, shots, visibility, rng):
    if type(shots) is not int or not 0 <= shots <= 100000: raise ValueError('bounded shot count required')
    p = likelihood_table((model,), (experiment,), (visibility,))[0][0]
    if p.denominator > np.iinfo(np.int64).max: raise ValueError('exact sampler denominator cap')
    return tuple(int(x < p.numerator) for x in rng.integers(0, p.denominator, size=shots))


def make_family(n, seed=0, prefix=False):
    if type(n) is not int or not 2 <= n <= 4 or type(seed) is not int or seed < 0 or type(prefix) is not bool:
        raise ValueError('bounded dimension, nonnegative seed and boolean prefix required')
    rng = np.random.default_rng(seed)
    pool = [''.join(x) for x in product('IXYZ', repeat=n) if sum(p != 'I' for p in x) == 2]
    rng.shuffle(pool)
    selected = pool[:4]
    known_prefix = pool[4]
    models = tuple(Model((((known_prefix, 1),) if prefix else ()) + ((p, sign),), sign)
                   for p in selected for sign in (-1, 1))
    locals_ = [''.join(axis if i == site else 'I' for i in range(n))
               for site in range(n) for axis in 'XZ']
    actions = [Experiment(q, q, 2) for q in locals_]
    # One direct correlation probe per unsigned hypothesis, using the public
    # candidate set. All four are represented before deterministic deduplication.
    for minus, plus in zip(models[::2], models[1::2]):
        found = False
        for q in locals_:
            prep, _ = transform(q, plus)
            experiment = Experiment(prep, q, sum(p != 'I' for p in prep) + 1)
            if expectation(plus, experiment) != expectation(minus, experiment):
                actions.append(experiment)
                found = True
                break
        if not found: raise AssertionError('candidate sign has no local readout probe')
    unique = {(a.prep, a.measurement): a for a in actions}
    return models, tuple(unique.values())
