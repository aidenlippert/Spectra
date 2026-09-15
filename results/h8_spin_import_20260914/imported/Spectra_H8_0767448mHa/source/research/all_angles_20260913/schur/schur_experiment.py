"""Bounded exact Schur-complement experiment with rational matrices."""
from fractions import Fraction as F
from itertools import combinations
import json
import math
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.marginal_transfer_verify import apply_word


def schur_lower(H, trial):
    """Return theta, residual^2, exact complement min bound, Schur endpoint.

    The complement bound is elementary Gershgorin after an exact orthogonal
    basis change; examples here use trial basis vectors, so no irrational
    change of basis is needed.
    """
    n = len(H)
    i = trial
    comp = [j for j in range(n) if j != i]
    theta = F(H[i][i])
    x = sum(F(H[j][i]) ** 2 for j in comp)
    # symmetric complement; Gershgorin lower endpoint is exact rational.
    mu = min(F(H[j][j]) - sum(abs(F(H[j][k])) for k in comp if k != j)
             for j in comp)
    if mu <= theta:
        return {"theta": theta, "residual_sq": x, "mu": mu,
                "endpoint": None, "status": "refuse_mu_le_theta"}
    endpoint = theta - x / (mu - theta)
    return {"theta": theta, "residual_sq": x, "mu": mu,
            "endpoint": endpoint, "status": "accepted"}


def eigen_min(H):
    # Validation oracle only; scipy is intentionally optional.
    try:
        import numpy as np
        return float(np.linalg.eigvalsh(np.array(H, dtype=float))[0])
    except Exception:
        return None


def run():
    cases = {
        "parent_like": ([[F(0), F(1, 10), F(0)],
                         [F(1, 10), F(2), F(1, 5)],
                         [F(0), F(1, 5), F(3)]], 0),
        # Deliberate failure: trial is an excited eigenvector.
        "excited_eigenvector": ([[F(0), F(0)], [F(0), F(1)]], 1),
        # Complement Gershgorin is too weak despite a useful residual.
        "naive_norm_failure": ([[F(1), F(1, 2), F(1, 2)],
                                [F(1, 2), F(3, 2), F(1)],
                                [F(1, 2), F(1), F(3, 2)]], 0),
    }
    out = {}
    for name, (H, trial) in cases.items():
        r = schur_lower(H, trial)
        r["exact_lambda_min_oracle"] = eigen_min(H)
        r["naive_norm_endpoint"] = float(r["theta"]) - math.sqrt(float(r["residual_sq"]))
        out[name] = {k: (str(v) if isinstance(v, F) else v) for k, v in r.items()}
    fixture = Path(__file__).resolve().parents[3] / "results/certificate_scaling/active_space_ladder_boys/h4/fixture.json"
    data = json.loads(fixture.read_text())
    basis = [s for s in range(1 << data["modes"]) if s.bit_count() == data["particles"]]
    idx = {s: i for i, s in enumerate(basis)}
    H = [[F(0) for _ in basis] for _ in basis]
    for term in data["hamiltonian"]:
        word = tuple(tuple(x) for x in term["word"]); coeff = F(term["coefficient"])
        for s in basis:
            result = apply_word(word, s)
            if result:
                target, sign = result
                if target in idx: H[idx[target]][idx[s]] += coeff * sign
    # RHF determinant from the saved independent upper-state support; no FCI
    # eigenvector is used. Full complement is 69x69; diagonal-only is cheaper.
    h4 = schur_lower(H, idx[15])
    h4["dimension"] = len(basis); h4["trial_state"] = 15
    comp = [j for j in range(len(H)) if j != idx[15]]
    h4["full_complement_entries"] = len(comp) ** 2
    h4["diagonal_mu"] = min(H[j][j] for j in comp)
    h4["oracle_min_eigenvalue"] = eigen_min(H)
    h4["h4_interval_reference"] = "results/certificate_scaling/active_space_ladder_references/h4/upper.json"
    out["molecular_H4_RHF_determinant"] = {k: (str(v) if isinstance(v, F) else v) for k, v in h4.items()}
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
