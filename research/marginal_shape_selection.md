# Stronger certificates by replacing diagonal shapes

The new exact million-site U4,t1,V1/2 lower bound is
`-3217627421190763/5000000000000000 = -0.6435254842381526` per site.
With the historical physical upper `-0.6106763470511881`, the interval width
is approximately0.0328491371869645. This improves the preceding accepted
lower by0.00013407268524356 and closes0.40649% of that preceding interval.

The construction still uses nine diagonal shapes and now uses exactly64
nonzero entries. No production verifier or matrix-size cap changed.

## Entry-count correction

The preceding `symmetric_diagonals_9/extra_31` certificate actually had60
entries. Its certificate and both exact energy receipts were correct;
earlier summaries incorrectly said64. The source catalog confirms that
the additional shape31 has four entries. That left room for a four-entry
addition, which was tested here. The affected summaries have been corrected.

## Five bounded trials

The previous exact dual mixture identified violated consistency directions.
All proposals used the same projector sources and ceilings, the same
250-evaluation cap, and a fresh CAR matrix reconstruction. The largest
fresh/affine discrepancy among these trials was2.14e-14.

| Trial | Proposed periodic lower/site | Matrix evaluations | Exact acceptance |
|---|---:|---:|---|
| Add four-entry shape8 | -0.643642560995 | 64 | Numerical only |
| Replace shape0 with11 | -0.643584209949 | 70 | Energy and family limit |
| Replace shape0 with20 | -0.643599582065 | 73 | Numerical only |
| Then replace shape3 with20 | -0.643522984238 | 73 | Energy and family limit |
| Instead replace shape3 with0 and8 | -0.643552485919 | 79 | Numerical only |

The last two trials start from the accepted shape0-to11 replacement.
The best final indices are `[1,2,4,5,6,7,31,11,20]`. Their nine coefficients,
the six local-profile parameters and the two projector penalties are stored
as exact rationals in the certificate. These five trials do not establish
the best support among all120 available directions.

The final energy replay reconstructed every Gram and all94 local Fock
blocks, accepting151 PSD checks in15.478 seconds with optional integer
congruence witnesses. Singular-safe fallback remains available. Exact
source/input hashes were checked after acceptance.

## The changed family's limit

A new18-vector local PSD mixture cancels all six profile derivatives and
all nine selected diagonal moments exactly. Its physical replay caps the
fixed family's attainable periodic lower certificates at approximately
`-0.6435228599405247`, leaving a gap of1.242976278974152e-7 above the accepted
periodic lower. This limit is conditional on the fixed sources, ratio,
ceilings and shape span; it is not a physical ground-energy upper bound.

The family-limit replay took1.961 seconds. The new mixture still violates93
of120 tested consistency directions. The largest remaining moments are
shape0 at approximately0.00192879, shape12 at-0.00185326 and shape8
at-0.00156874. They identify further separators, not accepted energy gains.
The support now genuinely fills64 entries, so further additions require
replacement or a different bounded representation.

## Negative-coupling transfer

The same final nine shapes and projector sources were applied to U4,t1,V=-1/2;
only scalar profiles, penalties and shape coefficients were retuned.
The exact million-site interval is approximately
`[-0.5510143532722551, -0.5240746138143995]` per site.

The lower improves on the previously accepted eight-shape transfer bound
by6.481601897812e-5, closing0.24002% of its interval. The physical upper was
recomputed with the matched transfer-state recipe, including a24-site exact
contraction inside its outward enclosure; it matches the preceding upper.
Combined exact replay took18.050 seconds. This is parameter transfer within
the one-dimensional nearest-neighbor family, not generic molecular transfer.

## Evidence and remaining scope

All outputs are under
`results/marginal_graded_hubbard8/joint_projector/signed_density/shape_selection/`.
`selection_summary.json` distinguishes all numerical trials from exact ones.
The best directory is `swap_3_for_20/`, with:

* `profile_joint_r1_2_certificate_accelerated_replay.json`: exact lower.
* `diagonal_family_limit_replay.json`: exact family cap and gap.
* `remaining_diagonal_overlap_moments.json`: exact missing constraints.
* `combined_summary.json`: comparison to the turn's starting interval.
* `transfer_V-1_2/matched_transfer_replay.json`: exact transfer interval.

Thirty-six focused tests and eight subtests passed in6.70 seconds. Production
modules are unchanged from the earlier504-test regression and14-test family
supplement. Only discovery/replay drivers and reports changed; original
driver versions used by earlier receipts are preserved in this directory.

No GPU was launched or left running. Generic molecular and long-range
transfer, higher dimensions, general representability and end-to-end
requested-accuracy scalability remain unproved. The finite improvements
and fast fixed-size verification do not solve general quantum chemistry.
