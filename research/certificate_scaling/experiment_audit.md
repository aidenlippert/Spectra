# Certificate-scaling experiment audit

This note defines the evidence needed before calling a new structural rule a
scaling result. It is deliberately stricter than showing a smaller exported
file on one favorable instance: discovery cost, replay cost, omitted
operator error, and transfer must all be measured.

## Baselines to preserve

The current baselines are useful but bounded.

* The direct H6 spin construction starts from the Hamiltonian alone, certifies
  an original-H interval of `1.931992366627265e-10 Ha`, and uses 400 explicit
  spin-zero source configurations, complement blocks of 200 and 168, 28
  response directions, and a 1,643,705-byte certificate. Its roughly 223 s
  construction time is not an asymptotic benchmark; explicit sector
  enumeration remains exponential.
* The H6 transfer reaches widths `2.05690e-10` and `1.54576e-10 Ha` at the
  recorded perturbations without factoring the joined 368-state block. This
  demonstrates a transfer mechanism, not generic scalability.
* The matched M=12 Reynolds construction has 30 charge orbits, 216 reduced
  PSD blocks, 2,599 PSD scalar variables, largest block 14, and an exact
  replay width `2.0645970172874137e-6` in model units. It is evidence for one
  symmetry family, not a generic chemistry result.

Every new result should report the same quantities, plus peak memory, wall
time split into dictionary generation/optimization/exact replay, generated
operator count, maximum block dimension, and whether any full sector or dense
reference was enumerated.

## Pass/fail gates

**G1 — Independent truth on small instances.** On sizes where the full
restricted problem can be solved, compare the restricted certificate with a
full reference run. The reference may be used for evaluation only, never as an
input to dictionary selection, coefficient proposal, or stopping. A proposed
certificate fails if its exact replay does not reproduce the claimed lower
endpoint or if the omitted-block improvement exceeds its claimed remainder.

**G2 — Same accuracy, stated units.** Compare methods at the same additive
energy width in Hartree (or explicitly stated model units), with both lower
and independently generated variational upper endpoints. Do not use relative
“fraction of optimum” or a coefficient sparsity percentage as a substitute.
For a family, report total error and per-site error; a fixed total target and
a fixed per-site target are different experiments.

**G3 — No hidden exponential work.** Count every determinant/configuration,
operator word, matrix element, eigensolve dimension, and candidate rejected
by the selector. A compact final certificate fails the scaling claim if
construction first enumerates the full fixed-number sector, constructs a dense
FCI/vector oracle, or scans an exponentially large dictionary. Small-instance
full references are allowed only in the separate G1 lane.

**G4 — Discovery is part of the result.** Measure the complete path from bare
Hamiltonian and physical metadata to the accepted certificate. Post-hoc
compression of a full optimum is a compression result, not an efficient
discovery result. Any structural rule must be fixed before held-out runs, or
all per-instance choices must be counted as search cost.

**G5 — Held-out chemistry transfer.** Fit or select the rule on one geometry,
coupling range, and basis. Apply it unchanged to held-out geometries and at
least one perturbed or asymmetric Hamiltonian, with a small Fe–S test if the
claim is chemically motivated. Report failures. Symmetry-matched M=12 or
H6-chain transfer alone cannot pass this gate.

**G6 — Omitted-operator certificate.** Hide predetermined operator blocks,
solve the restricted problem, and output an a priori upper bound on the
improvement available from the hidden blocks. Check the bound against the
full solve on small cases. The gate fails even when the restricted interval is
tight if the hidden improvement exceeds the predicted remainder.

**G7 — Replay and provenance.** A fresh exact checker must validate the
Hamiltonian identity, positivity/factor residuals, endpoint arithmetic, and
upper witness without solver state or FCI metadata. Record certificate bytes,
hashes, and all source counts. Unsupported or missing checks are indeterminate,
not passing.

## What would count as a real win

A credible scaling claim needs G1–G7 on a predeclared benchmark table and
must show that generated local blocks, peak memory, and discovery work grow
substantially more slowly than the full rank-2 construction at fixed additive
accuracy. A compact replay alone is insufficient. The strongest theorem target
is a computable omitted-operator bound tied to locality, symmetry, a gap or
correlation assumption, and the observed residual, together with an adaptive
algorithm that terminates when the bound meets the requested width.

The current artifacts support ambitious work toward that target. They do not
yet establish a universal structural condition, polynomial discovery, or
generic FeMoco/superconductivity capability. Competitor comparisons should be
made against measured implementations and stated certificate budgets rather
than claims that another group has universally stalled.

## Audit of the first three probes

The probes are useful as falsification and instrumentation, but none currently
passes the scaling gates. `operator_pricing.py` prices every unused one-body
candidate at every round and diagonalizes a projected trial for each one. Its
reported dictionary sizes are only 30 and 50, so the adaptive path is itself a
quadratic candidate scan at this scale. It also constructs the full fixed-
particle Hamiltonian and diagonalizes it to obtain the reference. That is
acceptable for G1, but the receipt must keep the reference work separate from
the discovery work. The experiment is an upper-subspace heuristic: it has no
SOS positivity or omitted-operator bound, and the one-body dictionary cannot
support a claim about general Gram-direction discovery.

`locality_compression.py` does provide a valid conservative lower-bound
construction for the stated Hubbard family, assuming the boundary hopping
operator norm bound. Its full-chain eigensolve is clearly labeled as an audit
reference. The central result is a negative one: the `2|t|b` penalty grows with
boundaries, so fixed cluster size does not give fixed absolute error. A future
claim must not report the exact diagonalization as part of certificate
generation and must test a boundary correction against this baseline.

`low_rank_compression.py` is only a numerical proxy. Its random normalized map
is not a CAR operator dictionary, and no exact coefficient identity, PSD repair,
or energy remainder is checked. The experiment can compare spectral decay with
entry sparsification, but it cannot support a chemistry certificate or a
discovery theorem until the map is replaced by an actual replayable Gram-to-CAR
coefficient map and the repair cost is bounded.

## Audit of the follow-up CAR and quantized probes

The follow-up `operator_pricing_sos.py` does use the repository's CAR Gram map,
but it remains a floating-point SDP proposal. Its `residual()` function
reconstructs coefficients from floating Gram matrices and reports a NumPy
maximum norm; it does not replay the rounded integer factors written by the
certificate exporter, and it does not convert the residual into a rational
lower-bound allowance. A residual below `1e-10` therefore cannot be called an
exact certificate. In addition, the greedy loop silently catches every
exception while scoring candidates, so an implementation error and a genuinely
infeasible block are indistinguishable. The receipt should report these as
explicit statuses and separate solver failures from infeasibility.

The updated low-rank fixture path is a meaningful improvement because it uses
the actual CAR `gram_map` on an H4 export. It still forms Gram matrices by
converting integer factors to floating arrays, diagonalizes in floating point,
and evaluates the truncated residual numerically. The `exact_full_factor_replay`
entry only replays the same factor against itself; it does not prove that any
truncated Gram has a PSD repair or an energy error bound. The result supports a
negative rank-frontier observation, not exact quantized compression.

The locality output field `certified_bracket_width` is misleading. The penalty
is an exact operator-norm bound for the omitted hopping, but the decoupled
cluster energies are obtained by numerical eigensolvers, so the lower endpoint
is not certified unless those local minima receive outward-rounded lower bounds.
The field should be renamed to `omitted_penalty_width` (or the implementation
should add rational local lower certificates), and `certificate_status` should
not imply more than “numerical cluster minima plus exact operator penalty.”

An independent small CAR check now confirms that the corrected exact mapper in
`low_rank_compression.py` agrees with the repository `gram_map` on every
degree-4 row for a mixed rational Gram (with the expected extra higher-degree
words retained by the exact mapper). The prior H4 artifact with the incorrect
left-word convention is explicitly invalidated and must not be used as evidence.
The corrected H4 path still needs its rank-cap campaign receipts and exact
checker results before any compression gain is reported.

The earlier per-block residual-domination probe is superseded for whole-
certificate claims. Its one-body row-sum argument is sound in its stated
scope, but the reported H4 gain (about 0.74%, one block out of 25) should be
treated as historical diagnostic data. The full-certificate replay now
reconstructs the complete residual first, charges all higher-body words by
exact coefficient-ℓ1, and applies the centered one-body bound only to the
portion it actually recognizes.

The full replay's self-tests pass under `python -S`, and an independent replay
of `results/marginal_reynolds/m4_d4_exact/certificate.json` produced a
rational structured lower endpoint and no fixed-sector enumeration. The
remaining limitation is theorem scope: this is an exact bound for a supplied
finite certificate and its recognized one-body residual, not a proof that
future compressed certificates exist or that discovery scales efficiently.

## Final campaign audit

The authoritative campaign report and interval summary are internally
consistent on the headline byte reductions: `(175912-119096)/175912` is
32.30% for square H4, and `(129032-112873)/129032` is 12.52% for rectangle
H4. The denominator-100,000 widths are `0.0012728552101868398 Ha` and
`0.0007310853907606697 Ha`; both are below the stated `0.0015 Ha` target. The
denominator-1,000,000 widths are also below target. The summary records 70
fixed-N sector states for the independently recomputed variational upper in
each case, while the exact lower replay is rational.

The geometry separation is correctly stated. The remote `heldout_gpu` folder
contains another representation of the same square geometry and is not
held-out molecular transfer. The rectangle is the only limited geometry
transfer reported for the denominator rule. Thus the 32.30% and 12.52% figures
pass as finite-basis certificate compression results, not as a held-out
discovery or active-space scaling result.

The report correctly labels post-processing as post-processing: existing
certificate discovery and upper validation costs remain, and no universal or
generic chemistry claim is made. The 37 fresh lower replays and lifecycle
record support the stated validation scope. Remaining blockers are discovery
cost, omitted-dictionary scaling, and larger or genuinely held-out active
spaces, all of which remain unmeasured.
