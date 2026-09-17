"""Executable seeds; every mutation becomes a separately retained source."""


def constructor(digits=12, subset=False, guard="1/10000000", margin="8/10000", backend="reference"):
    return f'''def propose(payload):
    rows = []
    for block in payload["blocks"]:
        energy, vector = eigenpair(block["H"], payload["den"], {subset})
        lower = rational(round(energy * 10000000000), 10000000000) - rational("{margin}")
        factor = cholesky(block["H"], payload["den"], lower + rational("{guard}"), {digits})
        rows.append({{"alpha_count": block["alpha"], "lower_Ha": str(lower),
                     "factor_denominator": 10 ** {digits}, "factor": factor,
                     "upper_vector": [int(round(x * 10000000000)) for x in vector]}})
    return {{"blocks": rows}}
'''


def adaptive_constructor(acquired=False):
    checker = "acquired_margin" if acquired else "reference_margin"
    return f'''def propose(payload):
    rows = []
    for block in payload["blocks"]:
        energy, vector = eigenpair(block["H"], payload["den"], True)
        lower = rational(round(energy * 10000000000), 10000000000) - rational("8/10000")
        for digits, fd in [(7, 10000000), (9, 1000000000), (12, 1000000000000)]:
            factor = cholesky(block["H"], payload["den"], lower + rational("1/1000000"), digits)
            if {checker}(block["H"], payload["den"], lower, factor, fd) >= 0:
                break
        rows.append({{"alpha_count": block["alpha"], "lower_Ha": str(lower),
                     "factor_denominator": fd, "factor": factor,
                     "upper_vector": [int(round(x * 10000000000)) for x in vector]}})
    return {{"blocks": rows}}
'''


GRAM_REFERENCE = '''def propose(payload):
    factor = payload["factor"]
    return [[sum(x * y for x, y in zip(a, b)) for b in factor] for a in factor]
'''

GRAM_TRIANGULAR = '''def propose(payload):
    factor = payload["factor"]
    n = len(factor)
    result = [[0] * n for i in range(n)]
    for i in range(n):
        for j in range(i + 1):
            result[i][j] = sum(x * y for x, y in zip(factor[i], factor[j]))
            result[j][i] = result[i][j]
    return result
'''

GRAM_FLINT = '''def propose(payload):
    return integer_gram(payload["factor"])
'''

GRAM_BAD = '''def propose(payload):
    factor = payload["factor"]
    return [[0] * len(factor) for row in factor]
'''

# These are known mathematical tools, not novelty claims. Source is proposed,
# checked by Lean, then its actual axiom report is stored with the method.
LEAN_PAIR = '''import Std
theorem cross_term_bounds (x y : Int) :
    -(x*x + y*y) ≤ 2*(x*y) ∧ 2*(x*y) ≤ x*x + y*y := by
  have hs (z : Int) : 0 ≤ z*z := by
    by_cases hz : 0 ≤ z
    · exact Int.mul_nonneg hz hz
    · exact Int.mul_nonneg_of_nonpos_of_nonpos (by omega) (by omega)
  have hp := hs (x+y)
  have hm := hs (x-y)
  simp only [Int.add_mul, Int.mul_add] at hp
  simp only [Int.sub_mul, Int.mul_sub] at hm
  have hc : y*x = x*y := Int.mul_comm y x
  omega
#print axioms cross_term_bounds
'''

LEAN_COMPOSITION = '''import Std
theorem compose_lower_errors (value lower first second : Int)
    (h : lower ≤ value + first + second) :
    lower - first - second ≤ value := by omega
#print axioms compose_lower_errors
'''
