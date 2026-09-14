"""Conventional finite-hypothesis inference from variance-threshold records."""
from fractions import Fraction as F
from math import erf, erfc, sqrt, isfinite
from experiments.v5_frame import BASE_LAMBDAS, dot, frame, null_probe


class Belief(dict):
    def __init__(self, values, workcount=0):
        super().__init__(values)
        self.workcount = workcount
        self.posterior_state_count = len(values)


def _validate(prior):
    if not isinstance(prior, dict) or not prior:
        raise ValueError('nonempty finite belief required')
    if any(type(h) is not tuple or len(h) > 3 or any(type(s) is not int or s not in (-1, 1) for s in h) for h in prior):
        raise ValueError('bounded sign hypotheses required')
    if len({len(h) for h in prior}) != 1 or len(prior) > 8:
        raise ValueError('common bounded hypothesis length required')
    if any(not isfinite(float(w)) or w < 0 for w in prior.values()) or abs(sum(prior.values()) - 1) > 1e-12:
        raise ValueError('finite normalized probability weights required')


def prior_empty():
    return Belief({(): 1.})


def extend_belief(prior):
    _validate(prior)
    if len(next(iter(prior))) == 3:
        raise ValueError('three-mode family exhausted')
    return Belief({h + (s,): float(w) / 2 for h, w in prior.items() for s in (-1, 1)},
                  getattr(prior, 'workcount', 0))


def nominal_likelihood(signs, alpha, probe, outcome):
    if type(alpha) is not F or not 1 <= alpha <= 2 or type(outcome) is not int or outcome not in (0, 1):
        raise ValueError('valid exact scale and binary outcome required')
    if type(probe) is not tuple or len(probe) != 4 or any(type(x) is not F for x in probe) or dot(probe, probe) != 1:
        raise ValueError('exact unit probe required')
    modes, _ = frame(signs)
    variance = F(1) + sum(alpha * lam * dot(probe, mode) ** 2 for lam, mode in zip(BASE_LAMBDAS, modes))
    argument = 20 / sqrt(2 * float(variance))
    # Zero mean in every hypothesis. The outcome is a magnitude threshold,
    # not a sign test and not a fictitious hypothesis-dependent mean.
    return erfc(argument) if outcome else erf(argument)


def update_belief(prior, prefix_length, alpha, probe, outcome):
    _validate(prior)
    if type(prefix_length) is not int or not 1 <= prefix_length <= 3 or any(len(h) != prefix_length for h in prior):
        raise ValueError('hypothesis length must match active source prefix')
    weights = {h: w * nominal_likelihood(h, alpha, probe, outcome) for h, w in prior.items()}
    mass = sum(weights.values())
    if not isfinite(mass) or mass <= 0:
        raise ValueError('impossible or invalid observation')
    return Belief({h: w / mass for h, w in weights.items()}, getattr(prior, 'workcount', 0) + len(prior))


def map_prefix(prior):
    _validate(prior)
    return max(prior, key=lambda h: (prior[h], h))


def choose_probe(prior):
    return null_probe(map_prefix(prior))


def current_target_decision(prior):
    _validate(prior)
    if not next(iter(prior)):
        raise ValueError('no acquired target')
    plus = sum(w for h, w in prior.items() if h[-1] == 1)
    return 1 if plus >= F(1, 2) else -1
