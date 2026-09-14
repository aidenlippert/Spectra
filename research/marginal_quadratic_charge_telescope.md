# Stronger certificates from compact quadratic charge telescopes

A two-product charge correction strengthens the exact lower certificate on
all three one-dimensional targets without increasing the largest local PSD
block. For the original million-site model, the open-chain energy per site
is enclosed by approximately [-0.6433546571919653, -0.6106763470511881].
The preceding interval was [-0.6434299980096883, -0.6106763470511881].

| Range-two coupling W | Accepted lower/site | Recomputed physical upper/site | Gain over preceding lower |
|---|---:|---:|---:|
| 0 | -0.6433546571919653 | -0.6106763470511881 | 7.534081772296e-5 |
| +1/10 | -0.6440112230530247 | -0.6114511605830021 | 5.589655587208e-5 |
| -1/10 | -0.6429668879422543 | -0.6099015335193740 | 9.607866993368e-5 |

All cases use U=4, t=1, V=1/2 and N=1,000,000. Physical uppers were
recomputed afresh, with the exact 24-site contraction enclosed before the
million-site transfer. The uppers agree exactly with the previous results.
The gain closes approximately 0.2300%, 0.1714%, and 0.2897% of the respective
preceding interval widths. These are modest, rigorously verified gains.

## What the new correction does

Let q_i=n_up,i+n_down,i-1. The five-site operator

    Y = q0 q3 - q1 q4

adds the six-site term

    gamma (Y_left-Y_right) = gamma (q0 q3 - 2 q1 q4 + q2 q5).

Its translated sum is exactly zero. It changes the local PSD certificate
while preserving the physical Hamiltonian, including when W=0. No extra
opening-boundary penalty is needed: the correction cancels on the periodic
chain before the existing physical opening deduction is applied.

The v7 certificate format stores this operator as a rational coefficient
and a canonical pair label. The new bounded helper supports all six
reflection-odd quadratic charge polynomials on five sites. In expanded
form this particular Y has 416 nonzero diagonal entries; the formula has
only two products. The separately supplied sparse diagonal still has nine
shapes and 64 nonzero entries. All 94 local Fock blocks, their combined
dimension 4096, and their maximum dimension 200 remain unchanged.

The winning coefficients gamma are -492/3125, -66823/500000 and
-180299/1000000 for W=0,+1/10,-1/10. The free range-two profile parameter
and the previous nearest profiles, projector penalties and sparse shape
coefficients were optimized simultaneously. The numerical searches took
81, 72 and 82 matrix evaluations under the existing 250-evaluation cap.
Fresh CAR reconstruction agreed with the affine numerical model within
3.02e-14. Exact lower acceptance subsequently used all 151 PSD checks,
with optional integer congruence witnesses. No GPU ran.

## Exact limit of the expanded family

The preceding free-profile dual mixtures have a nonzero expectation of the
new correction and cannot cap the expanded family. A v3 family verifier
now checks all six quadratic moments exactly, in addition to the six
nearest-profile gradients, nine sparse-shape moments, free range-two
profile moment, trace, positivity and fixed-projector fidelity inequalities.
Five quadratic directions are implied by the profile constraints; only
one independent scalar row is added. A separate exact all-Fock CAR replay
checks these identities and proves rank six for the six quadratic telescope
columns. The identities are 2(A-B), D-E, R2, R3, 2B and E, where A,B
are onsite-profile derivatives and D,E are nearest-density derivatives.
The bounded mixture size rises from
19 to 20 sources for nine sparse shapes. Old family versions retain their
existing source caps.

Fresh 20-vector rational mixtures were accepted for all three targets:

| W | Periodic gap from accepted lower to full-quadratic family ceiling/site |
|---|---:|
| 0 | 5.406925868495681e-8 |
| +1/10 | 6.144059153559798e-8 |
| -1/10 | 1.233973613464014e-7 |

These limits cover arbitrary coefficients of all six quadratic telescopes,
all reflected mean-correct nearest and range-two profiles, and the fixed
nine-shape sparse span, for the existing projector sources, ratio and
ceilings. They are upper limits on attainable LOWER certificates. They
are not physical ground-energy upper bounds. Changing the sources, ratio,
correction span or support can escape these caps.

## The next exact consistency obstruction

A separate exact probe evaluated the complete 52-element particle-hole-even,
reflection-odd basis of functions of five site charges. Each single-site
exponent is in {0,1,2}; this is complete because q takes values {-1,0,1}.
The six degree-two moments vanish exactly. Every one of the other 46 basis
moments is nonzero in each accepted dual mixture.

Ranking by absolute expectation divided by the local operator norm selects

    Y = q3^2 q4^2 - q0^2 q1^2.

Its Y_left-Y_right has norm 2. Its expectations are approximately
-0.0035524850, -0.0037840449 and -0.0033375248 for W=0,+1/10,-1/10.
Because q_i^2 is the indicator of an empty or doubly occupied site, this
is a compact consistency test for neighboring charge-excitation pairs.
Y has 384 nonzero diagonal entries but only two products of charge squares.
The probe independently checks the six-site expansion, reflection and
particle-hole identities over the full determinant spaces.

These are violated necessary translation-consistency conditions. They do
not establish a new physical positivity inequality or an additional energy
improvement. Optimizing the quartic correction remains the next bounded
experiment. General representability remains unresolved even if every
charge-diagonal moment is eventually enforced; spin-resolved and
noncommuting operators lie outside this probe.

## Artifacts and validation

Artifacts are in `results/marginal_graded_hubbard8/quadratic_charge_telescope/`.
Each target directory contains the v7 energy certificate, congruence
witnesses, accepting `range_two_replay.json`, accepting v3
`range_two_family_limit_replay.json`, and `charge_polynomial_overlap.json`.
The root `combined_summary.json` compares all three targets. All 83 recorded
source/input hash entries per target were checked against current files.
Pre-change sources are preserved under `pre_quadratic_sources/`.

Focused validation passed: 22 family/correction tests in 18.37 seconds,
plus the earlier 22 correction/range-two tests in 9.95 seconds. These include
full-Fock acceptance and rejection against an independent diagonal minimum,
exact periodic cancellation, malformed coefficient refusals, rejection of
a previously accepted dual under the new moment condition, and preservation
of the older mixture caps. The full regression passed: 556 tests and 96 subtests in 386.06 seconds.
The central validation receipt is `results/marginal_final_validation.json`.

The result extends the compact certificate method on the original model
and two finite range-two transfer examples. Arbitrary molecular and
long-range interactions, higher dimensions, general representability and
requested-accuracy scalability remain unproved.
