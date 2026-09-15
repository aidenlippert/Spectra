# Library integration and the H10 accuracy milestone

The original H10 rational Hamiltonian now has a locally rechecked complete
fixed-N interval of **0.479085203 mHa**. Both the orbital-rotation/MPS
upper and the full singlet-plus-nonsinglet lower were replayed with standard
Python integer/rational arithmetic. No full fixed-N determinant space enters
this construction or acceptance. The earlier 0.479085203 mHa GPU-discovered
certificate and all failed antecedent attempts remain preserved.

This is an adaptive follow-up to the frozen transfer campaign. The local-basis
MPS, transported one-/two-/three-body moments and changed numerical schedule
are explicit dependencies. It is not a success of the original timed-out
canonical random-state H10 rule, and it is not a new smaller SOS family.

## Measured library effects

SuiteSparseQR replaces a dense QR factorization of ideal-multiplier coefficient
columns. It preserves the numerical span checks, energy normalization and the
original exact acceptance requirements. This factors a constraint-coordinate
matrix, not a many-electron Hamiltonian matrix. cuDSS replaces only the shifted
normal preconditioner; the unshifted constrained solve is unchanged.

The following runs use the same host, FP64 precision, prepared arrays, zero
Gram start, 200 iterations and fixed mu=2. The setup costs are included. These
are fixed-work comparisons, not times to achieve a given certified accuracy.

| Case | Backend | Process seconds | Setup seconds | QR seconds | Iteration-kernel seconds | Peak host GB |
|---|---|---:|---:|---:|---:|---:|
| H8 | dense_200 | 18.945 | 7.489 | 7.146 | 8.824 | 0.598 |
| H8 | sparse_cpu_200 | 11.989 | 0.533 | 0.181 | 8.816 | 0.191 |
| H8 | sparse_hybrid_200 | 10.723 | 1.047 | 0.181 | 6.817 | 0.731 |
| H8 | sparse_cudss_200 | 10.835 | 1.322 | 0.182 | 6.535 | 0.864 |
| H10 | dense_200 | 216.205 | 169.564 | 168.384 | 39.589 | 3.982 |
| H10 | sparse_cpu_200 | 48.450 | 2.433 | 1.279 | 39.067 | 0.491 |
| H10 | sparse_hybrid_200 | 26.729 | 2.915 | 1.285 | 16.698 | 0.981 |
| H10 | sparse_cudss_200 | 26.333 | 3.190 | 1.314 | 15.925 | 1.110 |

The machine-readable receipt includes final objective/width differences,
numerical ranks and normal-map consistency errors for assessing equivalence.
Neither an optimizer status nor a floating-point bound prediction accepts an
energy certificate.

cuDSS and CHOLMOD were also compared with SciPy SuperLU on the actual H8/H10
normal matrices. Repeated-solve timings include RHS/output transfers. The first
probe accidentally charged CUDA context initialization to its first CPU setup;
that original result is preserved and the corrected probe records initialization
separately. Its repeated-solve measurements were unaffected by that setup issue.

FLINT rational arithmetic was tested separately from compiled integer Gram
products. Every comparison below requires exactly equal upper and lower
endpoints to the original standard-library receipt.

| H8 arithmetic | Lower replay seconds | Complete upper/lower replay seconds |
|---|---:|---:|
| h8_stdlib | 94.337 | 154.371 |
| h8_flint | 53.145 | 111.979 |
| h8_flint_gram | 91.377 | 151.374 |

These are paired measurements, not a broad timing distribution. FLINT rational
arithmetic is selected for the integrated experiment. The Gram-only variant
did not establish a useful improvement. Its focused tests include arbitrary
large integers, empty factors, invalid words, mixed charges and noninteger
factor rejection. The original accepting source files are unchanged.

The first integrated FLINT replay failed because the upper check had warmed a
CAR cache with Python Fraction objects, then the lower used FLINT rationals.
The retry loads the lower dependencies before selecting arithmetic and clears
that cache. A focused regression reproduces the warm-cache transition and checks
identical exact spin polynomials. The failed process and its cost are preserved.
The public replay caller also enforces the even-electron assumption required by
the integer-spin argument; all accepted fixtures satisfy it.

cuTensorNet was measured on an actual H10 MPS transfer with tensor shapes
(32,32), (32,2,63), (32,2,63), (2,2). It did not beat the existing sparse Numba
CPU transfer: about 42.7 microseconds for the CPU implementation, 120.2
microseconds for cached resident cuTensorNet, and 399.9 microseconds including
transfers. All compared values agreed to about 1e-16 relative error. This is
one contraction test, not a complete DMRG benchmark or a claim about larger
tensor networks. Its initial missing-CVXPY import failure is retained; adding
that existing transitive dependency enabled the retry.

## Complete-certificate runs and dependency cost

| Proposal | Proposal seconds | Proof seconds | Charged mixed-host component sum |
|---|---:|---:|---:|
| hybrid_adaptive | 612.610 | 403.310 | 1732.453 |
| sparse_cudss_adaptive | 424.700 | 269.340 | 1410.573 |

The shared local preparation contribution is **716.533 seconds**:
integrals, exact orbital rotation, fresh local MPS, upper acceptance, nonsinglet
construction, coefficient-map preparation and transported moments. The first
proof column is the standard-library lower replay; the second rechecks the
upper and uses FLINT rationals for the lower. The latter sum conservatively
charges upper verification again. The sums combine local Mac preparation and
remote A100-host computation; they are not measurements of one fresh complete
run on a single host. Earlier failed searches, benchmarks, installation and
independent replay are additional campaign costs, not omitted dependencies.

The fresh integrated H10 proposal starts again from zero Gram/ideal coordinates;
it uses no prior winning checkpoint or FCI vector. Exact local replay is the
authority for the reported final interval. Numerical FCI was used separately
as a reference in the parent campaign, not as a construction input.

## Cost, preservation and limits

The owned A100 instance is terminated. Launch through confirmed absence was
69.47 minutes at the observed $1.99/hour rate, for an estimated
**$2.30**, including boot and idle
time. The provider invoice is not available. Pre-existing instances and the
separate earlier GPU task were left untouched.

All remote job outcomes, setup failures, transfer receipts, source hashes,
installed versions and Linux process-time/memory records are retained. Linux
maximum resident sizes in the time files are KiB; local budget receipts use
bytes. Inner process durations are not added again to their enclosing remote
wall times. Earlier local research and failed H10 starts remain in the parent
all-attempt ledger.

This establishes a new H10 model certificate and working numerical-library
integrations. Broad transfer, favorable large-system scaling, competitive
complete cost and experimentally useful predictions still need evidence.

- [Exact local H10 result](/Users/aidenlippert/Documents/Spectra/results/transfer_solver_20260915/adaptive/h10_correlated_guide/hybrid_adaptive/exact_local_stdlib/interval.json)
- [Machine-readable measurements](/Users/aidenlippert/Documents/Spectra/results/nvidia_followup_20260915/summary.json)
- [Additional chemistry/library shortlist](/Users/aidenlippert/Documents/Spectra/research/nvidia_followup_20260915/LIBRARIES.md)
- [Broader transfer and physical-model results](/Users/aidenlippert/Documents/Spectra/research/transfer_followup_20260915/RESULTS.md)
