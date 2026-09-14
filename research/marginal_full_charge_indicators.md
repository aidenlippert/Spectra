# Complete indicator consistency: stronger bounds and a classical/quantum separation

Enforcing the remaining six empty/double-site indicator constraints strengthens
all four exact energy certificates. It also resolves a representation question:
the accepted local mixtures now have stationary classical extensions after
binary projection, but their mixed signed-charge moments still violate
translation consistency. Classical indicator extendibility therefore does not
establish extendibility of the underlying quantum marginal.

All results below concern open half-filled one-dimensional chains with
U=4, t=1, V=1/2 and N=1,000,000.

| W | Accepted lower/site, approximately | Recomputed physical upper/site | Gain over preceding lower |
|---|---:|---:|---:|
| 0 | -0.6432237895841160 | -0.6106763470511881 | 4.059518217284e-5 |
| +1/10 | -0.6438715054058292 | -0.6114511605830021 | 4.069229552820e-5 |
| -1/10 | -0.6428444456974698 | -0.6099015335193740 | 4.005864816572e-5 |
| +1 | -0.6617035788877922 | -0.6184244823693281 | 3.626556315036e-5 |

The physical uppers were recomputed afresh, including the exact 24-site
contraction enclosed before million-site transfer, and remain unchanged.
Each new lower was accepted against freshly reconstructed integer matrices.
No GPU ran.

## The full indicator correction space

Let q_i=n_up,i+n_down,i-1 and p_i=q_i^2. Each p_i is a binary indicator of
an empty or doubly occupied site. Every function of five p_i can be written
in square-free monomials. The reflection-odd subspace has dimension twelve:
two single-site differences, four pair differences, four triple differences
and two quadruple differences.

The single-site and pair directions were already enforced. The new six
operators are

    Y_S = product(p_i, i in S) - product(p_(4-i), i in S),
    S = {0,1,2}, {0,1,3}, {0,1,4}, {0,2,3},
        {0,1,2,3}, {0,1,2,4}.

A v9 energy certificate stores their rational coefficients in
`higher_charge_indicator_telescope`. Each contributes Y_left-Y_right on six
sites, so its translated sum cancels exactly before opening the physical
Hamiltonian. At W=0 the target remains the original nearest-neighbor model.

An exact rank calculation checks the entire 32-pattern indicator space and
all 4096 physical six-site determinants. The twelve indicator directions are
independent. The actual combined diagonal span—nine sparse shapes, six
quadratic directions, four charge-square pairs and six higher indicators—has
rank 25. Thus the new six directions add independent constraints even after
the actual sparse shapes are included.

The sparse correction still has 64 entries. The local PSD decomposition still
has 94 blocks, total dimension 4096 and maximum dimension 200. No PSD size cap
was increased.

## Numerical limitations and exact family ceilings

The first four searches reached the 120-iteration limit, despite improving
the proposed lower bounds. The original dual sampling radius 0.002 then
reported three sampled LPs infeasible and one unknown status. Those results
were preserved; they are not mathematical infeasibility certificates.

A second bounded optimization pass improved each lower further. W=1 declared
convergence after 94 evaluations; the other three again reached 120 iterations.
Each individual search retained the 250-matrix-evaluation hard cap. The first
and second passes used 240 evaluations in total per original target and 214
for W=1. Optimizer convergence is not claimed for the three original targets.

Widening the numerical dual sampling radius to 0.02 yielded exact rational
30-vector proposals for every target. Additional radii 0.01 (original targets)
and 0.005 (W=1) were tested. The selected 0.02 witnesses were tighter for the
original targets; 0.005 was slightly tighter for W=1. All selected witnesses
were subsequently accepted by independent exact physical replay.

The v5 family verifier checks all twelve indicator moments exactly, including
the old six, along with all six quadratic moments, nine sparse moments,
nearest and range-two profile conditions, trace, positivity and both fixed
projector fidelity inequalities. The source cap rises from 24 to 30 only in
this new mode. Older formats retain their existing caps.

| W | Periodic gap from accepted lower to full-indicator family ceiling/site |
|---|---:|
| 0 | 1.7942465185608435e-6 |
| +1/10 | 1.8776663979155050e-6 |
| -1/10 | 1.6584185985964500e-6 |
| +1 | 7.714695899635311e-7 |

These are conservative exact limits on further improvement within the fixed
projector sources, ratio, ceilings and sparse span, allowing arbitrary
coefficients of the full indicator and quadratic correction spaces and
reflected mean-correct profiles. They certify proximity to the family ceiling
at the displayed tolerances. They do not identify the exact optimum, certify
numerical convergence, or cap other source choices or larger correction spaces.
The family ceilings are not physical ground-energy upper bounds.

## An exact classical extension, with a remaining quantum obstruction

The accepted dual gives a nonnegative six-bit probability distribution for
(p0,...,p5). A separate exact construction verifies all 32 equalities between
its five-bit prefix and suffix marginals. If pi(x) is that common marginal,
set the transition from x to its shifted successor y to

    P(x,y) = probability(the corresponding six-bit word) / pi(x).

Every pi(x) is positive in the four actual witnesses. The emitted sparse
32-state transition matrices have nonnegative entries and rows summing to one.
Exact checks verify shift compatibility, pi P=pi and recovery of all 64
six-bit probabilities. Iterating this stationary kernel yields a classical
order-five Markov extension. All 32 stationary states and all 64 six-bit
patterns have positive probability in these witnesses. Independent checks
cover an i.i.d. binary model, an alternating model and rejection of an
inconsistent distribution.

This extension certifies only the empty/double indicators. The exact probe
of the complete 52-element particle-hole-even reflection-odd signed-charge
polynomial space finds sixteen zero moments and 36 nonzero moments for every
witness. The twelve indicator directions and four other quadratic charge
directions are now canceled; mixed signed-charge/indicator constraints remain.

For W=0,+1/10,-1/10, the strongest normalized separator is

    Y = q0 q4 (p2 p3 - p1 p2).

Its local telescope has norm one and Y has only 64 nonzero diagonal entries.
Its expectations are approximately 0.0001442513924, 0.0001292162586 and
0.0001599940010. For W=1, the strongest is

    Y = q0 q4 (p3 - p1),

with local norm two, 128 nonzero entries and expectation
-0.0003867305860. These nonzero expectations rule out a stationary extension
matching the full six-site quantum marginal: translated copies of each
Y_left-Y_right must have zero expectation in such an extension.

Thus a positive local quantum mixture can have a fully extendible classical
binary projection while still failing necessary consistency of its signed
charges. This is a concrete limitation of the indicator representation.
It does not rule out extending a single marginal into an inhomogeneous state,
or assert a new general representability theorem.

Testing these mixed signed-charge operators is the next bounded direction.
The 64-entry separator fits the raw sparse cap by itself, but combining it
with the existing nine shapes may exceed their union cap; its compact
formula avoids assuming that union is still allowed. No lower improvement
using these new mixed corrections has been tested yet.

## Artifacts and validation

Outputs are under `results/marginal_graded_hubbard8/full_charge_indicators/`.
The root `combined_summary.json` uses the final `polished/` directory inside
each target. Those directories contain accepting energy and v5 family
receipts, `charge_polynomial_overlap.json`, and
`indicator_markov_extension.json` with the exact stationary measure and
transition edges. All 90 source/input hash entries per final target were
checked. `full_indicator_rank.json` records the two rank checks.

The first accepted energy results remain in the parent target directories.
Their failed initial dual searches are retained in `initial_dual_diagnostic.json`.
Both later dual sampling candidates are retained in each `polished/`
directory. Changed earlier source versions are preserved under
`pre_indicator_sources/`.

Focused validation passed: 22 energy/correction tests in 9.84 seconds and
nine family tests in 10.89 seconds. These check the independent full-Fock
expansion, exact PSD acceptance and rejection, all indicator symmetries and
periodic cancellation, malformed inputs, refusal of an older dual under the
new constraints and preservation of older source caps. The full regression passed: 584 tests and 96 subtests in 501.18 seconds.
The central validation receipt is `results/marginal_final_validation.json`.

The broader goal remains active. General quantum representability, arbitrary
molecular or long-range transfer, higher dimensions and requested-accuracy
scalability remain unproved. The W=1 result remains a finite one-dimensional
transfer example with a wider physical interval than the original model.
