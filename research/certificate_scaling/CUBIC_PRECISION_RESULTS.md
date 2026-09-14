# Cubic precision campaign: H6 and H8 reach the interval target

The frozen straight-chain H4, H6 and H8 fixtures now have independently
replayable ground-energy intervals narrower than 0.0016 Ha. H10 does not.
This is a numerical/certification advance, not the structural scaling theorem
in the active goal. Full cubic dictionaries and their discovery costs remain.

## Best uncompressed ladder

These are RHF-canonical STO-3G hydrogen chains at 1.4 Angstrom spacing,
CAS(n,n), with rational electronic Hamiltonians. Energies omit nuclear
repulsion. Bounds concern the supplied finite-basis Hamiltonian, not basis
error, experimental accuracy, reaction barriers or finite-temperature physics.

| Fixture | Modes / particles | Certified width (Ha) | Target | Original SOS rows | Additional residual factors | Gram entries built | Original discovery + export (s) |
|---|---:|---:|---|---:|---:|---:|---:|
| H4 | 8 / 4 | 0.0000001530024 | pass | 499 | 0 | 17,776 | 7.259 |
| H6 | 12 / 6 | 0.00005499439 | pass | 1,116 | 298 | 217,268 | 259.309 |
| H8 | 16 / 8 | 0.0007505851 | pass | 2,574 | 696 | 1,206,464 | 721.563 |
| H10 | 20 / 10 | 0.02334631 | fail | 6,233 | 1,350 | 4,674,900 | 1,215.936 |

The exact Fraction endpoints and hashes are in
`results/certificate_scaling/cubic_precision/intervals/final_h{4,6,8,10}.json`.
H4 uses the original coefficient-L1 replayer; H6/H8/H10 use the additional
spectral residual witnesses. H4 ran Clarabel locally; the other selected
runs used SCS on the two existing A10 hosts, on CPU. Budgets differ, so the
four timings are not an asymptotic exponent estimate or a fair solver contest.

Original SOS nonzeros are 15,234 / 124,432 / 699,343 / 3,395,544. The spectral
witnesses add 9,418 / 56,224 / 224,050 integer entries for H6/H8/H10 and
95,108 / 554,403 / 2,190,186 bytes. These extra factors are part of the proof,
not free compression. No determinant upper or source factors enter discovery.

## What changed

1. **Measured numerical losses.** Raw primal/dual/Gram arrays are saved before
   exact export. Diagnostics separately measure the intended residual,
   equality error, negative Gram spectrum, PSD clipping, factor rounding,
   multiplier rounding, and final exact residual. Row/column conditioning
   significantly improved H6 and H8 in the recorded trials.
2. **Exact fixed-number quotient.** Linear factors in the mixed cubic
   dictionary can be represented by cubic factors modulo the number ideal.
   Creation and annihilation have different sector denominators. This is a
   proved equivalence with the body-two multiplier, not a locality heuristic.
   It reduces H6 Gram entries from 217,268 to 209,780; the quotient H6 trial
   had a worse bound than the selected nonquotient conditioned trial.
3. **Stronger exact residual bounds.** For a k-body coefficient matrix on
   exterior k-tuples, V >= ell I implies R_k >= ell binomial(N,k). Exact
   Gershgorin bounds already make conditioned H8 pass, at width 0.001481342.
   Rational Gram decompositions plus an exact residual Gershgorin check
   tighten that to 0.0007505851. Neither floating eigenvalues nor an arbitrary
   supplied Gram matrix are trusted.
4. **Sparse multiplier repair.** With H, b and all SOS factors fixed, a sparse
   LP redistributes the residual through the number ideal. Exact replay
   chooses the improvement. H6 wedge width improves from 0.00009657247 to
   0.00008440851. Unconditioned H8 improves from about 0.01140778 to
   0.009918200, still failing. Large cancelling multipliers can worsen the
   original coefficient-L1 bound, so these results require the named wedge
   replayer.
5. **Replay reuse and exact Gram symmetry.** The original verifier now exposes
   its already validated exact residual to extensions; its public result and
   refusal gates are preserved. Exact Gram products use symmetry and
   contiguous columns. On the same H8 artifact, the complete spectral interval
   replay decreased from 43.85 to 36.46 seconds locally with identical exact
   endpoints. This is one observed comparison, not a broad speedup claim.

See `cubic_number_frame.md`, `wedge_residual_bound.md`, and
`wedge_spectral_bound.md` for the identities and their scope.

## Negative evidence and remaining bottleneck

The unconditioned H8 SCS run produced an original interval of about 0.026285
Ha. Conditioning reduced PSD-clipping coefficient loss from 0.024790 to
0.002506 Ha at the same 600-second solver budget. Spectral residual bounds
recover more accuracy but cannot repair an underconverged original solution.

H10 constructed 112,701 coefficient rows and 2,347,238 map nonzeros. It spent
154.24 seconds building maps, 673.24 seconds in the solve, and 388.13 seconds
in the original export. It achieved only 175 solver iterations. Its numerical
candidate lower was -12.3917353, compared with a reference upper of
-12.3682995; the exact original width was 0.03221293. The stronger replay
still gives 0.02334631. This does not establish the optimum of the cubic
relaxation, so it cannot distinguish a cone limitation from poor convergence.

Both initial H6 Clarabel primal trials and the quotient explicit-dual trial
hit their 360-second external limits without certificates. Those failures
are retained. Successful SCS runs reported `optimal_inaccurate`; that status
is acceptable only because the exported rational proof is checked exactly.

Post-discovery H6 pruning has a separate validated artifact: the unconditioned
source reduces to 238 rows and 33,400 nonzeros while its original exact interval
is 0.001319972 Ha. The conditioned source retains 252 rows and 34,216 nonzeros
at width 0.0003966475. All original solve and pruning/replay costs remain
charged; neither demonstrates direct sparse discovery.

Conditioned H8 also compresses from 2,574 to 576 original SOS rows and
191,936 nonzeros, at width 0.0008204347 Ha. Its original certificate is
2,125,203 bytes; the residual proof adds 696 factors and 559,994 bytes.
Three thresholds were tested; all trial costs and the final independent
`python -S` interval replay are recorded under `factor_structure_h8`.
Replay was not faster in this measurement, and the original 721.56-second
full-dictionary discovery remains charged.

## Upper witnesses and validation cost

FCI is validation-only. Integrals are reconstructed from the frozen rational
Hamiltonian, without rerunning SCF or importing a different orbital gauge.
Upper amplitudes are rationalized and their Rayleigh quotient is replayed
exactly against that same H. The H8 full selected reference has 2,468 nonzero
integer states; H10 has 31,752, selected from 63,504 FCI coefficient pairs.
Zeros after rounding are recorded separately from top-k truncation.

The bounded streaming upper backend performs O(K times Hamiltonian terms)
integer CAR actions without a K-by-K retained matrix or action cache. It
replayed H10's 227,026,800 word/state checks in about 16 seconds. It does not
remove exponentially growing FCI discovery or witness support costs. Old
small-witness oracle limits remain unchanged. These are separate upper
artifacts, not Hilbert-space-free upper discovery.

## GPU gate and retained infrastructure

The two existing A10s remain warm intentionally and are controlled by SSH/CLI.
No new instance was launched for this campaign and neither was torn down.
The existing combined listed rate is $2.58/hour; this is not a final invoice.

Warmed float64 Torch eigensolves on the A10 were slower than one CPU thread
for batches 8x368 and 16x112. An H10-sized 8x725 batch was 1.52 times faster
on GPU (0.229 vs 0.350 seconds); 8x225 remained slower. These isolated kernel
measurements exclude transfer, assembly, and the rest of SCS. They support
investigating selective large-block acceleration, not an end-to-end claim.

## Verification and reproducibility

The research tests cover exact number-frame relations, raw physical map
reconstruction, conditioning, primal/dual agreement, wedge signs and sector
lifts, malformed spectral factors, component minima, streaming CAR signs,
large-reference refusal boundaries, and fixed-factor multiplier repair.
The focused suite passes 20 tests; the symbolic/orbit suite passes 17.
Full exact interval replay runs with `python -S`. Remote/download hashes
match for all 42 retained result files. Initial source snapshots, final
source snapshot and manifests live in `results/lambda_runs/cubic_precision`.

Example full replay:

```sh
python -S research/certificate_scaling/cubic_interval_replay.py \
  --certificate results/lambda_runs/cubic_precision/downloaded/h8_scs_quotient_conditioned/certificate.json \
  --reference results/certificate_scaling/cubic_precision/reference_h8_full_checked/h8/upper.json \
  --method spectral \
  --proof results/certificate_scaling/cubic_precision/wedge_spectral/h8_conditioned/witness.json \
  --out /tmp/spectra-h8-interval.json
```

An early agent spectral draft used incorrect tuple indexing and is explicitly
marked INVALID. It was replaced, independently audited, and covered by new
k=2/k=3 regression controls; only the corrected subdirectory witnesses count.

## Next falsifiable work

First separate H10 convergence from the relaxation's actual limit. Persist
and resume canonical SCS x/y/s states under exact problem-identity checks;
CVXPY's default cache must not be assumed to retain inaccurate iterates.
Compare certified width and full cost after additional iterations. In
parallel, test direct factor selection without building every pricing entry.
A small post-hoc factor set is evidence for a compression experiment, not the
required input-only structural condition.

A second mathematical route is positive finite-range amplitudes and exact
local-energy dynamic programming for stoquastic fixed-N chains. Existing
spin-chain parent controls are already known in this repo. Extending them to
local fermionic hopping and larger amplitude families remains work, and
positivity/local range alone does not guarantee accuracy. See
`positive_amplitude_structural_route.md` for the derivation and prior art.

The goal remains ACTIVE: no structurally guaranteed accurate discovery
algorithm for realistic growing molecular active spaces has been established.
