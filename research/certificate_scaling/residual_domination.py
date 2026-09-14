"""Exact rational local residual bounds for number-conserving CAR remainders.

For a one-body residual R=dGamma(A), the coefficient-l1 bound charges every
Hermitian hopping word separately.  A row-sum bound instead uses
||A||_infty and ||dGamma(A)||_{N} <= N||A||_infty.  Both bounds are rational,
local, and require no global fixed-N matrix.  Small Fock enumeration is used
only as an independent audit of the inequality, never as the verifier.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def l1_bound(edges):
    return 2 * sum(abs(c) for _, _, c in edges)


def row_bound(modes, particles, edges):
    rows = [F(0) for _ in range(modes)]
    for i, j, c in edges:
        rows[i] += abs(c); rows[j] += abs(c)
    return particles * max(rows, default=F(0))


def action(modes, particles, edges):
    """Independent audit matrix for dGamma(A) on the fixed sector."""
    states = list(combinations(range(modes), particles)); where = {s:k for k,s in enumerate(states)}
    mat = [[F(0) for _ in states] for _ in states]
    for col, state in enumerate(states):
        occ = set(state)
        for i, j, c in edges:
            for src, dst in ((j, i), (i, j)):
                if src not in occ or dst in occ: continue
                ns = tuple(sorted((occ-{src})|{dst}))
                # creation after annihilation gives the usual CAR sign.
                left = sum(1 for x in state if x < src)
                removed = tuple(x for x in state if x != src)
                left += sum(1 for x in removed if x < dst)
                sign = F(-1 if left % 2 else 1)
                mat[where[ns]][col] += c*sign
    return mat


def gershgorin_norm(mat):
    # For a real symmetric matrix, max absolute row sum is a safe spectral norm.
    return max((sum(abs(x) for x in row) for row in mat), default=F(0))


def case(name, modes, particles, edges):
    naive = l1_bound(edges); structured = row_bound(modes, particles, edges)
    audited = gershgorin_norm(action(modes, particles, edges))
    if audited > structured:
        raise AssertionError((name, audited, structured))
    return {'name': name, 'modes': modes, 'particles': particles,
            'edges': [[i, j, str(c)] for i, j, c in edges],
            'naive_l1': str(naive), 'local_row_bound': str(structured),
            'audit_fixed_sector_row_bound': str(audited),
            'improvement_ratio_l1_over_local': str(naive/structured) if structured else 'inf',
            'verifier': 'rational row sums; fixed-sector matrix only audit'}


def molecular_h4_case():
    """Apply the rule to the actual quantized H4 factor fixture.

    The full residual is retained in the coefficient-l1 budget; only its
    degree-two one-body part is regrouped.  This avoids pretending that the
    higher-body tail is controlled by a one-body theorem.
    """
    from experiments.marginal_coefficient import coefficient_rows, gram_map, dagger
    from research.certificate_scaling.low_rank_compression import quantized_spectral_factor
    from experiments.marginal_symbolic import word_product
    import numpy as np
    fixture = Path(__file__).resolve().parents[2] / 'results/marginal_molecule_stress/h4_square_degree3_certificate.json'
    data = json.loads(fixture.read_text()); modes = int(data['modes']); particles = int(data['particles'])
    lookup = {w:i for i,w in enumerate(coefficient_rows(modes))}; denom = int(data['denominator'])
    total_l1 = F(0); total_structured = F(0); blocks = 0; improved = 0
    def exact_coefficients_correct(words, matrix):
        out = {}
        for i, left in enumerate(words):
            for j, right in enumerate(words):
                value = matrix[i][j]
                if not value: continue
                for word, sign in word_product(dagger(left), right):
                    out[word] = out.get(word, F(0)) + value * sign
        return out
    for block in data['blocks']:
        words = [tuple(tuple(x) for x in w) for w in block['words']]
        fac = np.asarray(block['factor'], dtype=float) / denom
        gram = fac.T @ fac
        vals, vecs = np.linalg.eigh(gram); rank = min(1, len(vals))
        qint, qden = quantized_spectral_factor(gram, rank)
        qfactor = [[F(int(qint[i,j]), qden) for j in range(qint.shape[1])] for i in range(qint.shape[0])]
        qgram = [[sum((qfactor[i][k]*qfactor[j][k] for k in range(qint.shape[1])), F(0)) for j in range(len(words))] for i in range(len(words))]
        original = exact_coefficients_correct(words, [[sum((F(int(block['factor'][k][i]))*F(int(block['factor'][k][j])) for k in range(len(block['factor']))), F(0))/denom**2 for j in range(len(words))] for i in range(len(words))])
        compressed = exact_coefficients_correct(words, qgram)
        delta = {w: compressed.get(w,F(0))-original.get(w,F(0)) for w in set(compressed)|set(original)}
        l1 = sum((abs(v) for v in delta.values()), F(0)); total_l1 += l1
        one = {}; rest = F(0)
        for w,v in delta.items():
            if len(w) == 2 and w[0][0] == 1 and w[1][0] == 0:
                one[(w[0][1],w[1][1])] = one.get((w[0][1],w[1][1]),F(0))+v
            else: rest += abs(v)
        tr = sum((one.get((i,i),F(0)) for i in range(modes)),F(0))/modes
        rows = [sum((abs(one.get((i,j),F(0))-(tr if i==j else F(0))) for j in range(modes)),F(0)) for i in range(modes)]
        grouped = rest + min(particles,modes-particles)*max(rows,default=F(0)) + abs(tr*particles)
        total_structured += grouped
        blocks += 1
        improved += grouped < l1
    return {'name':'actual_h4_rank1_all_blocks', 'fixture':str(fixture), 'blocks':blocks,
            'naive_total_l1':str(total_l1), 'structured_onebody_total':str(total_structured),
            'improved_blocks':improved, 'improvement_ratio':str(total_l1/total_structured) if total_structured else 'inf',
            'scope':'one-body component grouped; higher-body residual remains coefficient-l1'}


def run():
    # Asymmetric held-out: degree is concentrated but coefficients are spread.
    held_out = [(0,1,F(1)), (0,2,F(3,4)), (1,3,F(1,2)), (2,3,F(1,4))]
    # Adverse: a star saturates the degree bound, so no claimed universal gain.
    adverse = [(0,1,F(1)), (0,2,F(1)), (0,3,F(1))]
    out = {'scope': 'one-body CAR residual domination; local rational verifier',
           'cases': [case('asymmetric_held_out',4,2,held_out),
                     case('adverse_star',4,2,adverse), molecular_h4_case()],
           'complexity': {'verifier': 'O(|E| + M)', 'memory': 'O(M)',
                          'global_fock_dimension_constructed': False,
                          'audit_only': 'O(binomial(M,N)^2) for the two four-mode checks'}}
    dest = Path(__file__).resolve().parents[2] / 'results/certificate_scaling/residual_domination'
    dest.mkdir(parents=True, exist_ok=True)
    (dest/'receipt.json').write_text(json.dumps(out, indent=2)+'\n')
    return out


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
