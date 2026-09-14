# Complete signed-charge certificates expose an off-diagonal obstruction

The complete particle-hole-even, reflection-odd signed-charge correction
space strengthens all four exact energy lower bounds. More significantly,
the accepted local mixtures now admit stationary classical extensions of
their symmetry-averaged signed-charge laws, yet still fail an explicitly
verified hopping consistency condition. This isolates an off-diagonal
obstruction beyond classical charge consistency.

All targets are open half-filled chains with U=4, t=1, V=1/2 and
N=1,000,000. The interaction W couples charges two sites apart.

| W | Accepted lower/site, approximately | Recomputed physical upper/site | Gain over preceding lower |
|---|---:|---:|---:|
| 0 | -0.6431820851905259 | -0.6106763470511881 | 4.170439359016e-5 |
| +1/10 | -0.6438293270332963 | -0.6114511605830021 | 4.217837253284e-5 |
| -1/10 | -0.6428033207035662 | -0.6099015335193740 | 4.112499390368e-5 |
| +1 | -0.6616379487206246 | -0.6184244823693281 | 6.563016716760e-5 |

Every lower passed fresh exact all-Fock replay. Physical uppers were
recomputed, including the exact 24-site contraction within its 160-bit
enclosure before million-site transfer, and remain unchanged. No GPU ran.

## A complete charge function space

Each site charge q_i takes values {-1,0,1}. On five sites there are 243 charge
patterns. Projection onto functions even under q -> -q and odd under spatial
reflection has dimension

    (243 + 1 - 27 - 9) / 4 = 52.

The four traces count the patterns fixed by identity, particle-hole inversion,
reflection, and their composition. A new exact helper represents these
functions using 52 disjoint pattern classes: each has value +1 on a canonical
pattern and its charge inverse, -1 on their reflections, and zero elsewhere.
There are 208 nonzero pattern occurrences in total.

An exact identity check reconstructs all 52 charge-polynomial basis elements
in this pattern basis on every physical five-site determinant. The previous
sixteen charge/indicator directions are contained in this space; the remaining
36 are independent of them and of the actual nine sparse shapes. The combined
five-site diagonal rank is 61, checked exactly on all 1024 determinants.

The v10 energy format adds a bounded `signed_charge_telescope` field. Its
six-site contribution is Y_left-Y_right, whose translated sum cancels exactly.
The numerical chart fixes the old diagonal profiles and compact corrections,
then varies all 52 signed-charge pattern coefficients, nine sparse-shape
coefficients, two hopping-profile parameters and two projector penalties.
Fixing the redundant old charge coordinates does not shrink the covered
operator span. Forty-four pattern coefficients became nonzero in each result.

The existing sparse field remains at 64 entries. The 94 local PSD blocks,
their total dimension 4096 and maximum dimension 200 are unchanged. The
Hamiltonian at W=0 remains the original nearest-neighbor model.

## What was and was not resolved numerically

All four searches reached their 240-iteration limit. Each retained the hard
250-matrix-evaluation cap and passed fresh physical reconstruction within
3.56e-15 of the affine chart. Optimizer convergence was not established.
The disjoint coordinates were chosen to reduce coordinate overlap; no
controlled conditioning or speedup claim is made.

The first dual searches, using only nearby low-energy vectors, reported the
sampled LP infeasible for all four targets. These are preserved numerical
failures, not exact infeasibility certificates. Adding all 4096 physical
occupation determinants produced bounded pools of 4548 candidates. This
includes trivially feasible witnesses and avoids restricting the dual to
a small neighborhood of the current primal search.

The selected numerical supports had fewer than 66 independent columns. A
bounded exact rectangular solve handles those supports, checks every equation
and refuses inconsistent or dependent candidates. The production verifier
then independently reconstructs all physical expectations. It checks all
52 pattern moments, all earlier moments, trace, positivity, profile conditions
and both projector fidelity inequalities. Its source cap is 66; the accepted
witnesses use 59, 58, 59 and 58 vectors, respectively.

| W | Periodic gap from accepted lower to full-charge family ceiling/site |
|---|---:|
| 0 | 3.905668618334712e-5 |
| +1/10 | 4.256669255112121e-5 |
| -1/10 | 5.164377015198710e-5 |
| +1 | 8.273009059499841e-5 |

These exact ceilings cover the full signed-charge function space, the fixed
nine-shape sparse span and reflected mean-correct profiles, with the existing
projector sources, ratio and overlap ceilings. They are limits on attainable
LOWER certificates, not physical energy upper bounds. Their gaps are much
larger than in preceding smaller families. The current numerical lower bounds
must not be described as tightly optimized or as exhausting this larger family.

A separate CLARABEL/CVXPY check was attempted on W=0, with bounded active
cones and solver settings. It exited with code 137 without producing a
proposal. Its failure cause is not established. The last live observation was
at 285 seconds; no final elapsed time or solver convergence is inferred.
The process was gone before the wall guard could act. This failed attempt
supports no improvement or performance conclusion and was not used in any
accepted certificate.

## Classical charge extension after explicit symmetry averaging

The raw dual witnesses need not themselves be particle-hole or reflection
symmetric. Their raw five-charge prefix/suffix distributions differ, with
L1 discrepancies approximately 0.1121, 0.09565, 0.1142 and 0.1095. The extension
construction therefore explicitly averages each six-charge distribution over
identity, charge inversion, spatial reflection and their composition.

For these averaged distributions, all 243 prefix/suffix equalities hold
exactly. Conditional six-charge probabilities define sparse transitions on
243 five-charge states. Exact replay checks stochastic rows, shift
compatibility, stationarity and reconstruction of all 729 six-charge
probabilities. Zero-mass states receive harmless shift-compatible rows.
The four witnesses have 205, 203, 207 and 205 positive stationary states,
with 443, 439, 435 and 435 positive six-charge patterns.

Iterating the stationary kernel gives a classical order-five Markov extension
of the averaged signed-charge law. Independent checks cover an i.i.d. ternary
model and rejection of inconsistent prefix/suffix data. This result concerns
classical charges after the stated averaging. It does not certify fixed total
particle number, spin distributions, fermionic coherences or a quantum extension.

## A verified off-diagonal separator

Define the spin-summed hopping observable

    B(i,j) = sum_spin(c_i^dagger c_j + c_j^dagger c_i).

The five-site operator Y=B(0,3)-B(1,4) yields

    T = Y_left-Y_right = B(0,3)-2B(1,4)+B(2,5).

The twelve CAR words of T were evaluated exactly on all 4096 determinants.
Hermiticity and invariance under physical particle-hole, reflection and spin
exchange transformations were checked. Its periodic translated sum is zero.
Its norm is exactly eight: the triangle bound for six disjoint mode-pair
hoppings gives eight, and an emitted integer 64-component vector attains
the +8 eigenvalue exactly.

The accepted full-charge duals have these nonzero moments:

| W | Exact-replayed expectation of T, approximately | Absolute expectation / norm |
|---|---:|---:|
| 0 | 0.002147081989623470 | 0.000268385248702934 |
| +1/10 | 0.002572210325880284 | 0.000321526290735036 |
| -1/10 | 0.002690734705416714 | 0.000336341838177089 |
| +1 | 0.019802676676793113 | 0.002475334584599139 |

Because T is invariant under the averaging transformations, these violations
persist for the corresponding symmetry-averaged quantum density matrices.
A stationary quantum extension matching such a local density matrix would
have zero expectation of Y_left-Y_right. Thus the classical signed-charge
extension exists while the averaged quantum marginal still fails a necessary
off-diagonal consistency condition.

This does not refute extending a single marginal into an inhomogeneous state,
or solve general representability. It identifies a compact next correction:
a mean-zero range-three hopping profile [gamma,-2 gamma,gamma] on the local
six-site window. No stronger energy certificate using this operator has yet
been optimized or accepted. Its compatibility with the existing symmetry
blocks is checked at the operator level, but a production certificate path
and independent tests are still required before claiming an energy gain.

## Artifacts and validation

The root is `results/marginal_graded_hubbard8/full_signed_charge/`.
`combined_summary.json` compares all targets. Each target directory contains
the v10 energy certificate, congruence hints, accepted energy and v6 family
receipts, `charge_markov_extension.json`, and `coherent_overlap.json` with the
norm-attaining vector. All 120 recorded source/input hashes per target were
checked. Root `span_rank.json` and `basis_identity.json` record the exact
basis results. Earlier changed sources and the pre-anchor dual script are
preserved under `pre_signed_sources/`.

The first failed dual attempts remain in `initial_dual_diagnostic.json`.
The failed conic attempt is recorded under W_zero, including its execution
exit code and lack of a proposal.

The full regression passed: 596 tests and 96 subtests in 454.40 seconds.
Focused validation passed 20 energy/correction tests in 13.22 seconds and
12 family/correction tests in 10.94 seconds. Tests cover full-Fock physical
expansion, PSD acceptance and rejection, canonical pattern/symmetry checks,
malformed input refusal, rejection of the preceding indicator-only witness
under the new conditions, and preservation of older source caps.

General quantum representability, arbitrary molecular or long-range transfer,
higher dimensions and requested-accuracy scalability remain unproved. The
goal remains active, with both a quantitative optimization gap and a concrete
coherence constraint still to address.
