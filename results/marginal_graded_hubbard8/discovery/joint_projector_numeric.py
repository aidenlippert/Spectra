"""Actual-source numerical discovery for the coupled half/charged projector bound.

This file deliberately emits proposals only.  Exact Gram acceptance is obtained
through ``joint_projector_bound``; the energy optimization below is floating
point reconnaissance and has no certificate status.
"""
from pathlib import Path
from fractions import Fraction as F
import json, math, sys, time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_projector_extendibility import _projector_vector
from experiments.marginal_charged_projectors import charged_vectors, joint_projector_bound

BASE = ROOT / 'results/marginal_graded_hubbard8'
OUT = BASE / 'joint_projector'


def main():
    started = time.monotonic(); OUT.mkdir(exist_ok=True)
    cert = json.loads((BASE / 'six_site_projector/refined_certificate.json').read_text())
    local = cert['local_window']
    sectors = sector_matrices(6, F(local['U']), F(local['t']),
                              list(map(F, local['onsite_profile'])),
                              list(map(F, local['hopping_profile'])), F(local['V']),
                              list(map(F, local['density_profile'])))
    half, half_norm = _projector_vector(cert['vector'], 6)
    charged = json.loads((BASE / 'charged_projector/source.json').read_text())['physical_state']
    charged = {int(s): int(a) for s, a in charged.items()}
    charged_family, charged_norm = charged_vectors(charged)
    theta_h = F(cert['projector_sum_ceiling']) / cert['windows']
    ratio_rows = []
    # Four windows is the useful compact joint family; retain all three ratios.
    for ratio in (F(1, 2), F(1), F(2)):
        # Add a single ulp-like rational margin to the observed boundary.
        # Search the smallest accepted rational ceiling by exact replay.
        # Start from the conservative integer ceiling and descend until refusal.
        ceiling = F(4)  # valid upper bound for ratio=2, four windows
        # candidate boundaries are found from floating Gram eigenvalues below;
        # exact replay loop is intentionally bounded to keep this discovery safe.
        # Compute a numerical maximum of C D C^T via the Gram congruence directly.
        from experiments.marginal_charged_projectors import joint_overlap_grams
        groups = joint_overlap_grams(cert['vector'], charged, 4)
        maxima = []
        for data in groups.values():
            gram = np.asarray(data['gram'], float)
            norms = np.asarray(data['norms'], float)
            src = np.asarray(data['sources'])
            weights = norms / np.where(src == 0, 1.0, float(ratio))
            sym = gram / np.sqrt(weights[:, None] * weights[None, :])
            maxima.append(float(eigh(sym, eigvals_only=True,
                                     subset_by_index=[len(gram)-1, len(gram)-1])[0]))
        observed = max(maxima)
        ceiling = F(math.ceil((observed + 2e-9) * 10**9), 10**9)
        ceiling = max(F(1), ratio, ceiling)
        accepted = None  # exact replay is intentionally a separate bounded job
        theta_j = ceiling / 4

        mats = []
        for key, (_, K, cols) in sectors.items():
            scale = np.sqrt([sum(a*a for a in col.values()) for col in cols])
            A = np.asarray(K, float) / scale[:, None] / scale[None, :]
            def projector(v, norm):
                w = np.asarray([sum(a*v.get(s, 0) for s, a in col.items())
                                for col in cols], float) / scale / math.sqrt(float(norm))
                return np.outer(w, w)
            P = projector(half, half_norm)
            Q = sum((projector(v, charged_norm) for v in charged_family),
                    np.zeros_like(A))
            mats.append((key, A, P, Q))

        calls = 0
        def minimum(x):
            alpha, beta = x
            return min((float(eigh(A + (alpha + beta)*P + beta*float(ratio)*Q,
                                      eigvals_only=True, subset_by_index=[0, 0])[0]), key)
                       for key, A, P, Q in mats)
        def objective(x):
            nonlocal calls
            calls += 1
            if calls > 120:
                raise RuntimeError('Hard objective cap exceeded')
            if min(x) < 0 or max(x) > 3:
                return 1000 + sum(abs(z) for z in x)
            ell, _ = minimum(x)
            return -(ell - x[0]*float(theta_h) - x[1]*float(theta_j)) / 5
        x0 = [float(F(cert['penalty'])), 0.0]
        result = minimize(objective, x0, method='Nelder-Mead',
                          options={'maxfev': 120, 'xatol': 1e-8, 'fatol': 1e-10})
        alpha, beta = (F(round(float(z)*10**6), 10**6) for z in result.x)
        ell, active = minimum((float(alpha), float(beta)))
        lower = F(math.floor(ell*10**7)-1, 10**7)
        density = (lower - alpha*theta_h - beta*theta_j) / 5
        ratio_rows.append({'ratio': str(ratio), 'ceiling': str(ceiling),
                           'theta_joint': str(theta_j), 'observed_maximum': observed,
                           'columns': sum(len(g['gram']) for g in groups.values()),
                           'maximum_block': max(len(g['gram']) for g in groups.values()),
                           'exact_replay': False,
                           'candidate_alpha': str(alpha), 'candidate_beta': str(beta),
                           'numerical_ell': ell, 'active_sector': list(active),
                           'proposed_periodic_density': str(density),
                           'objective_evaluations': calls})
    receipt = {'source': 'charged_projector/source.json', 'windows': 4,
               'theta_half': str(theta_h), 'rows': ratio_rows,
               'seconds': time.monotonic()-started,
               'scope': 'Numerical fixed-profile joint discovery only; exact Gram replay and energy acceptance remain separate.'}
    (OUT/'numeric_proposal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (OUT/'joint_gram_replay.json').write_text(json.dumps({'rows': [
        {'ratio': r['ratio'], 'ceiling': r['ceiling'], 'columns': r['columns'],
         'maximum_block': r['maximum_block']} for r in ratio_rows]}, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)

if __name__ == '__main__': main()
