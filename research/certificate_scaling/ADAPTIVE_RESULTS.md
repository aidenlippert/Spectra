# Adaptive sparse-factor discovery — 2026-09-12

The new search escapes the proved equal-weight quadratic pair restriction on
square H4. It does not yet beat full quadratic SDP, achieve the 0.0015 hartree
interval target, or establish an efficiently checkable structural condition
for compact certificates.

## What changed

`adaptive_factor_pricing.py` starts from Hamiltonian coefficients, mode count,
and particle number. It solves a restricted coefficient LP, builds moment
matrices from its dual, and proposes unequal-weight square factors on two,
four, eight, sixteen, or unrestricted operator-word supports. New directions
are not filtered by Hamiltonian overlap. The cubic lane adds higher-degree
operator words and body-two number-ideal multipliers.

Every accepted lower endpoint comes from a fully rationalized factorization,
number-ideal identity, and exact CAR residual replay. Eigenvalues and LP
objectives guide discovery; they are not accepted physical certificates.
Known upper witnesses are attached only after discovery finishes and are
re-evaluated against each certificate's Hamiltonian.

Pricing still constructs the complete coefficient-to-Gram maps and complete
numerical dual moment matrices. Greedy support selection has no omitted-family
optimality guarantee. This is a direct-discovery experiment, not a proof of
polynomial scaling. No existing solved factors, upper states, or many-body
basis enumeration enter discovery. The independent H4 upper controls use
70 states; the H6 control is a previously generated FCI-derived 200-state
witness. Their generation is not claimed to scale.

## Main comparison

All 32 remote experiments finished and their exported certificates passed
exact lower-and-upper interval replay. None passes the 0.0015 hartree target.
Five local exploratory certificates also replay successfully. Best remote
results across both batches:

| Fixture | Best adaptive interval width (Ha) | Interval bytes | Existing full quadratic SDP width (Ha) | SDP interval bytes |
|---|---:|---:|---:|---:|
| Square H4 | 0.141082726552 | 64,988 | 0.005036643277 | 25,251 |
| Rectangle H4 | 0.029386140487 | 58,421 | 0.000857081242 | 19,694 |
| H6 | 3.486464894363 | 127,415 | 0.013029347927 | 127,594 |

The H4 winners use quadratic width-four factors; H6's winner uses the local
cubic dictionary with width-eight factors. These are finite, stored
Hamiltonians. The adaptive method does not improve the accuracy/size tradeoff
over the existing SDP controls. The cubic expansion and longer optimized runs
did not change that conclusion.

The square H4 width-four lower endpoint is -3.471471331702718, above the
previous exact equal-weight-family ceiling -4.459007306545172. That is a
certified escape from this specific family obstruction. Unequal-weight
pairs alone improve the lower endpoint only to -4.45803126892303 in this
budget. This does not imply a new obstruction for all two-word factors.

Both full pair-cone numerical controls returned `optimal_inaccurate`; their
rationalized primal certificates are valid but do not establish the optimal
value or a ceiling for that cone.

## Optimization ablation

The LP dominates elapsed time in many successful searches. We tested retaining
only one coefficient from each Hermitian pair, with residual weight two off
the diagonal, and removing inactive atoms every five rounds. Exact replay
always reconstructs the full polynomial.

Square H4, quadratic width four, 80 pricing rounds, batch 64:

| Rows | Pruning | Exact lower (Ha) | Wall seconds |
|---|---|---:|---:|
| Full | Off | -3.4714713317 | 71.97 |
| Half | Off | -3.4722209197 | 69.02 |
| Full | On | -4.4590087797 | 7.88 |

The half-row-plus-pruning run completed 200 rounds in 23.40 seconds but only
reached -4.3807513836. It added 13,312 cumulative atoms and removed 12,596.
These are single runs, not repeat timing statistics. Reduced LP size is not
an end-to-end discovery win when bound quality collapses. Pruning remains off
by default. Dual nonuniqueness and reintroducing removed constraints are
plausible mechanisms; the ablation implicates aggressive pruning much more
strongly than row reduction.

## Verification and provenance

Six focused controls pass, covering pricing signs, exact exported bounds,
higher-degree residual rows, initial full/half LP equivalence, Hermitian
multipliers, and pruning bookkeeping. Research scripts changed; production
code did not. No broad production regression claim is made.

The remote campaign comprises 22 adaptive searches and two full pair-cone
controls in the original batch, followed by six optimized searches and two
ablation runs. All subprocesses exited successfully; several searches stopped
at a time/LP limit and returned their best earlier certificate. Successful
export is not solver convergence or cone optimality. All 130 downloaded
result-file hashes and all four remote result archive hashes matched.
`results/certificate_scaling/adaptive_intervals/campaign_summary.json` records
all 32 intervals and the best-by-fixture selection. Source archive SHA256:

- v1: `3a898a344f77a0d965484f5e0b4a8a269e1fa36112c63127383b3943c81d829b`
- v2: `308913f89908bf04085cabcd6f5c60c71b7501ee9231f596dd62e9dd362a6f07`

The first source archive and manifests are in
`results/lambda_runs/adaptive_pricing/`; the second frozen source and ablations
are in `results/lambda_runs/adaptive_pricing_v2/`. Never substitute the current
working source for either executed archive. Exact interval certificates and
replay receipts are under `results/certificate_scaling/adaptive_intervals/`.

Both existing A10 hosts were operated through SSH/API, with three independent
jobs per lane. These LP/SDP jobs use host CPUs; no GPU acceleration is claimed.
Hosts remain active intentionally at the user's request, at $2.58/hour
combined. No instance launch or teardown was needed for this campaign.

## Next scientific step

Keep the historical constraints when improving the restricted solver. Test
whether small PSD blocks with jointly adjustable coefficients make faster
accuracy progress than adding fixed rank-one directions. Measure accuracy,
certificate size, and total discovery cost against the existing full-SDP
controls on both H4 geometries and H6. An eventual compactness claim must
also account for the complete pricing maps or replace them with a certified
bounded-cost pricing mechanism. The current evidence provides neither a
competitive advantage over prior SOS/v2RDM work nor a solution of FeMoco or
finite-temperature superconductivity.
