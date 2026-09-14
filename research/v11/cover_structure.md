# V11 cover-structure headroom audit

**Superseded negative inference:** the factual V8 artifact description below
is retained, but the conclusion against a small family/one-sweep construction
was refuted by a root experiment. A union of the two ordinary partitions and
one sweep produces a valid14-group bound0.000998348, which meets the actual
0.001 requirement. Matching the old0.000965 bound is unnecessary. This does not
yet establish complete-cost headroom. See the authoritative
[next preflight and correction](NEXT_COVER_PREFLIGHT.md).

## Scope and evidence

This is a read-only audit of the V8 exact overlap-cover artifact against the
V11 conventional lower-HS rejection gate.  It uses `research/ACTIVE_STATUS.md`,
`results/v8/fractional_probe.json`, `results/v8/REPORT.md`,
`experiments/v8_fractional_cover.py`, and the V11 baseline profile.  The V8
artifact is a development result, not a held-out result; no novelty claim is
made here.

## The actual near-threshold winner

The only residual that crossed the tolerance at a lower degree is row 10 of
`fractional_probe.json`: `n=3`, XXZ, `gamma=2`, `T=1/2`, order 11, 30 residual
terms.  Its integrated cover bound is
`19761823751632003/20480000000000000000` (about `0.000965`).  The ordinary
partition alternatives on this residual are approximately `0.001783` (l1),
`0.001105` (first-fit), and `0.001017` (weighted), so the crossing is a narrow
bound improvement rather than a large margin.

The accepted cover is not a small reusable motif.  It has 61 groups, 295
nonzero coefficient incidences in the exact witness (the proposer accounting
reports 69 candidate groups, 319 incidences before the retained witness), 575
exact pair checks in the final checker, 356 square operations, and 295 exact
coefficient additions.  Groups are mostly 3–5 terms, with repeated labels
across groups.  The coefficients are large, unrelated rationals produced by
coordinate minimization and reconciliation; they are not a fixed small set of
fractions or a translational pattern.  The first groups already mix labels such
as `IIZ,IXX,IYX,YIY,YZX`, then `XZY,YZY,ZXX,ZYX,ZZI`, with later groups reusing
the same labels in different combinations.  This is a dense fractional
incidence pattern, not a closed-form splitting identity.

The construction itself used 16 cyclic sweeps and 5,104 coordinate updates,
5,704 candidate-family pair checks, and 0.01648 seconds; its internal exact
cover check took about 0.00331 seconds.  The complete cover path took about
0.01932 seconds for this residual.  By comparison, ordinary partition
construction took 0.000566 seconds.  Thus the accepted lower degree removes one
Taylor layer but pays roughly 34x the partition construction time on the
winning residual.  The complete V8 replay measured the overlap arm at 2.09–
2.25x slower than the integer Taylor baseline, despite order 12 to 11.

## Structural conclusion

The only mathematically visible structure is the supplied anticommuting graph:
each group is a clique, and overlap lets one coefficient be split among several
cliques.  The objective is a convex norm sum, so the split depends on the
current magnitudes of *all* other coefficients in the incident groups.  The
recorded coordinate trace decreases gradually (`5.0515` to `4.9272` over 16
sweeps), which is evidence of numerical optimization rather than convergence
to a fixed symbolic rule.  Replacing it with equal splitting, sign splitting,
or one pass of a local formula would not be justified by the certificate: the
checker requires exact reconstruction and each group’s exact squared-norm
upper bound.

The existing C5 half-cover identity is a known baseline and does not explain
this certificate.  Nor does Y-parity/bipartiteness supply the missing
coefficient split; it only restricts the operator graph.  Across all 48 V8
residuals, overlap tightened every norm, but only this one crossed tolerance,
while the full arm remained slower.  This gives no evidence that a fixed motif
would recur with enough margin to amortize its search or checker cost.

## Bounded candidate decision

The initial audit proposed no static candidate. Its proposed condition that a
candidate match the row-10 witness was unnecessarily strong; only the fixed
tolerance is required. The
actual accepted witness supplies neither: its 61-group incidence pattern and
large rational coefficients depend on the residual.  A reusable 3-term or
5-term overlap template would therefore be an unverified heuristic and would
need the same expensive exact checking, while losing the only observed
degree-11 crossing is likely.

The V11 lower-HS gate can safely reject cover proposals whose conventional
lower bound is already outside the useful near-threshold window.  It cannot
make the overlap construction cheap.  On the measured evidence, complete cost
remains the critical gate: the cover's norm tightness does not realize a net
benefit, and there is no supported small-motif/closed-form splitting procedure
to carry forward.
