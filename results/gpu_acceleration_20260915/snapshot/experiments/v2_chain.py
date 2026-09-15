"""Exact bounded enumeration of the parity-source complementary-risk theorem."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
from experiments.v2_policy import certificate, verify


def _validate(r, k, v):
    if type(r) is not int or not 1 <= r <= 6 or type(k) is not int or not 0 <= k <= r:
        raise ValueError('enumeration limited to 1 <= r <= 6 and 0 <= k <= r')
    if not isinstance(v, F) or not 0 <= v <= 1:
        raise ValueError('exact visibility in [0,1] required')


def source_eigenvalues(r, v=F(1, 2)):
    _validate(r, 0, v)
    return (((1 + v) / 2 ** r, 2 ** (r - 1)),
            ((1 - v) / 2 ** r, 2 ** (r - 1)))


def predicted_risk(r, k, v=F(1, 2)):
    _validate(r, k, v)
    return F(1, 2) - v / 2 ** (r - k + 1)


def enumerate_risk(r, k, v=F(1, 2), with_count=False):
    _validate(r, k, v)
    hypotheses = [a for a in product((0, 1), repeat=r) if all(a[i] == 0 for i in range(k))]
    action_risks = []
    additions = 0
    # All actions are checked, including those contradicting already known axes.
    for action in product((0, 1), repeat=r):
        joint = {y: {s: F(0) for s in (-1, 1)} for y in (-1, 1)}
        for axes in hypotheses:
            for sign in (-1, 1):
                plus = (1 + sign * v * (axes == action)) / 2
                joint[1][sign] += plus / (2 * len(hypotheses))
                joint[-1][sign] += (1 - plus) / (2 * len(hypotheses))
                additions += 2
        action_risks.append(sum(min(joint[y].values()) for y in (-1, 1)))
    risk = min(action_risks)
    return (risk, additions) if with_count else risk


def verify_chain(r, k, v=F(1, 2)):
    return enumerate_risk(r, k, v) == predicted_risk(r, k, v)


def write_summary(path=None):
    rows = []
    for r in range(1, 7):
        risks = []
        eigenvalues = source_eigenvalues(r)
        assert all(value >= 0 for value, _ in eigenvalues)
        assert sum(value * multiplicity for value, multiplicity in eigenvalues) == 1
        for k in range(r + 1):
            risk, operations = enumerate_risk(r, k, with_count=True)
            assert risk == predicted_risk(r, k)
            risks.append(risk)
            rows.append({'r': r, 'k': k, 'risk': str(risk), 'verified': True,
                         'likelihood_mass_additions': operations, 'actions_checked': 2 ** r})
        gains = [risks[k] - risks[k + 1] for k in range(r)]
        assert all(gains[k + 1] == 2 * gains[k] for k in range(r - 1))
    lower_bound_checks = []
    for t in range(6):
        c = certificate('none', t)
        assert verify(c)
        lower = max(F(0), F(1, 2) - F(t, 8))
        assert F(c['value']) >= lower
        lower_bound_checks.append({'r': 2, 'shots': t, 'coupling_error_lower_bound': str(lower),
                                   'exact_optimal_error': c['value']})
    out = {'rows': rows, 'adaptive_lower_bound_checks': lower_bound_checks,
           'marginal_doubling_verified': True,
           'claim': 'bounded exact enumeration corroborates the general written theorem; no physical validation'}
    path = Path(path) if path else Path(__file__).resolve().parents[1] / 'results/v2/chain_certificates.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + '\n')
    return out


if __name__ == '__main__':
    print(json.dumps(write_summary(), indent=2))
