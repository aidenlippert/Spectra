# Interaction geometry: channels and orbital cuts

## Coverage

This run covers routes 13–16 on the repository's H4/H6/H8 active-space
molecular fixtures. It forms the quartic coefficient tensor
`T[(p,q),(r,s)]` directly from fixture Hamiltonian rows, computes its SVD,
and rationalizes rank-k factors at denominator `10^8`. The complete receipt is
`results/all_angles_20260913/interaction_geometry/channel_rank_metrics.json`.

## Findings

Canonical fixtures have 8, 12, and 16 spin orbitals (H4/H6/H8). On this
three-size tested ladder, their
quartic channel singular tails (Frobenius norm) at ranks 1/4/8/16 are:

| fixture | r=1 | r=4 | r=8 | r=16 |
|---|---:|---:|---:|---:|
| H4 | 1.9555 | 1.3688 | 1.1592 | 0.5815 |
| H6 | 2.5960 | 2.1679 | 1.7739 | 1.4485 |
| H8 | 3.1186 | 2.7675 | 2.2458 | 1.9905 |

The tail increases with system size at fixed rank on this tested ladder; this
is evidence against (rather than an asymptotic theorem ruling out) a fixed
number of global density channels. At rank 16 the exact
rationalized tensor L1 residuals are 6.10 (H4), 26.83 (H6), and 60.17 (H8),
which are safe coefficient residual bounds for this tensor only. They are not
SOS Gram bounds and do not certify an energy interval. Full-rank rounding
residuals are below `6e-6` in all cases.

The Boys/localized H4/H6 runs give nearly identical singular values, but this
should be read as an empirical observation for these reshaped arrays. A
general orbital rotation is unitary on one-particle space, yet the displayed
normal-ordered upper-index tensor and its coordinatewise SVD need not transform
unitarily without antisymmetrizing/completing the pair basis. The rational L1
residual changes: at rank 16 it improves H4
from 6.10 to 5.85 and H6 from 26.83 to 26.31; at rank 1 it worsens H4 and
improves H6. This is a real localization tradeoff in a coordinate-sensitive
safe bound, not a reduction in intrinsic channel rank.

## Hypothesis and negative control

Collective Coulomb factors plus local residual certificates could help only if
the long-range collective part is removed in a basis where the *coefficient*
L1 tail is subextensive. These data show the unqualified global-rank version
fails: both numerical tail and rational safe residual increase with H size.
The useful next test is geometry-aware block/chordal factorization, with each
block residual certified separately and summed; integral compression must not
be treated as SOS compression.

The exact-safe ingredient used here is simple: every retained rank-k tensor
factor is rounded to integers over `10^8`, reconstructed exactly as a rational
pair tensor, mapped back to every normal-ordered quartic word (including all
index permutations present in the fixture), combined by exact `Fraction`
arithmetic, and then summed coefficientwise in absolute value. This is a
valid induced L1 bound for the omitted quartic operator component because each
CAR monomial has norm at most one. It is not a full-Hamiltonian bound: any
quadratic/scalar residual and any SOS-factor transport cost must be added.

As positive control, full-rank rational reconstruction has residual below
`6e-6` on every fixture. As adverse control, a random Frobenius-matched rank-1
channel is substantially worse; its exact residual is recorded under
`adverse_random_rank1_tensor_l1` in the JSON receipt. These controls test the
reconstruction and show that SVD ordering, rather than arbitrary factor choice,
is doing the useful compression.
