# Constructive compression investigation — 16 September 2026

**The requested major breakthrough is not achieved.** The new contribution is
a working exact tensor spectral-witness checker and a precise diagnosis of a
failed compact lower-bound construction. No new competitive solver, scaling
theorem, physical prediction, or general many-body solution is claimed.

The experiment stayed on the fully interacting open 2×4 Hubbard ladder with
U/t=8 and N_up=N_down=4. This is an eight-site lattice model, not molecular H8.
Its declared sector has 4,900 configurations. The target is a total energy
interval of 0.001t. Earlier accepted results remain intact.

## What is now checked

A cold charge-MPS calculation, built from local Hamiltonian terms and a charge
counter seed, produced an exactly replayed upper energy

\[
 U=-3.0256171860355323\ldots\,t.
\]

The stored state has 3,878 nonzero integer tensor entries, maximum bond 64.
The exact upper replay used 1,290,192 integer multiplications, at most 400
transfer entries, 0.3993 seconds as a fresh process and 20,447,232 bytes peak
RSS. It did not enumerate configurations or construct a many-body matrix.
This modest state compression is not a complete compressed proof.

The polynomial checker also completed an accepted, deliberately loose
lower-bound validation on this same model:

\[
 L=-9.9038193275\,t,
 \qquad U-L=6.878202141464468\ldots\,t.
\]

Its exact acceptance score was about 0.10913605 against the strict threshold
1. All eight recurrence steps were replayed. That bound is much weaker than
previous retained lower bounds and is solely an acceptance-path validation.
It is not the campaign accuracy result. The recurrence stores 26,208 tensor
entries across nine vectors; maximum replay RSS was 67,518,464 bytes and the
fresh process took 9.1565 seconds. No size or speed advantage follows.

Twenty-four focused tests pass. They include explicit CAR/particle–hole
binding on the small validation sector, tensor contractions, charge flow,
exact square-root enclosures, an adversarial erased-ground-component example,
acceptance and refusal on the actual eight-site proposal, and altered seed,
Hamiltonian hash, charge-flow and sector checks. The fresh exact replay imports
no NumPy, SciPy, quimb or PySCF. It is separate from the numerical constructor,
not a third-party replication or a wholly separate Hamiltonian implementation.

## The target failure is measured exactly

For the proposed approximately 0.0009t lower-to-upper separation, the bond-64
recurrence has exact reconstructed residual norm upper bounds:

| Step | Residual norm upper |
|---|---:|
| 1 | 2.863998815882951e-9 |
| 2 | 4.556568455882371e-9 |
| 3 | 0.01756020121410984 |

At step 3, the geometric-envelope residual allowance alone is
**2.1361500569160237**, against available overlap **1**. The contribution is
computed from exact squared norms and rational outward square-root bounds.
No claimed SVD discarded weight enters acceptance.

Extending this same prefix cannot make this particular geometric-envelope
gate pass at the same scaling. This is a narrow exact diagnosis. It is not a
proof that every bond-64 construction, a sharper polynomial gate, another
seed or the broader tensor method must fail.

The state energy itself is already adequate for the finite-model target:
combining its new upper with the preserved *enumerated* lower from the
positive-cone campaign gives width 0.00030609706446771515t. That comparison is
not a new enumeration-free complete certificate. It shows why improving the
upper alone does not remove the outstanding proof dependency.

## Numerical screens, including unsuccessful work

These are numerical diagnostics, not accepted lower bounds. The entries and
energy below refer to the end of cold state discovery. The timings include
the run's subsequent failed filter or preparation attempt.

| Bond cap | Cold numerical energy / t | State nonzeros | Fresh process seconds | Peak RSS bytes | Target gate |
|---|---:|---:|---:|---:|---|
| 16 | −2.90448028346310 | 564 | 3.6087 | 112,263,168 | Failed |
| 32 | −3.02127438239366 | 1,643 | 7.1059 | 254,296,064 | Failed |
| 64 | −3.02561718603606 | 3,878 | 34.6748 | 500,908,032 | Failed |
| 96 | −3.02591051325414 | 6,303 | 37.0211 | 598,638,592 | Failed |

The bond-64 preparation lost its numerical lower-overlap guarantee before
finishing the requested 300 power steps. The bond-96 preparation retained a
positive diagnostic overlap but the subsequent filter still failed. Its cold
state alone uses more entries than the full 4,900-amplitude reference. These
measurements do not justify another blind increase in tensor size.

The bond-96 warmup produced a worse Rayleigh energy than its cold state; the
screen records both. No bond-96 upper was promoted to an exact accepted result.

## Mathematical work across the routes

The reviewed [derivation](DERIVATION.md) proves the following sufficient rule.
For a seed with guaranteed ground overlap gamma, arbitrary stored recurrence
vectors and rigorous residual bounds eta_j,

\[
2z^k\|v_k\|+
\frac{2}{1-z^2}\sum_{j=1}^kz^j\eta_j<\gamma
\]

certifies the stated lower energy. The Hubbard seed has gamma=1. This replaces
the previous dense matrix-positivity check by exact tensor contractions.
The missing construction is one that makes the weighted residual small enough
while retaining manageable representation and complete cost.

Other routes were examined but did not close that gap:

* Relative resolvent preconditioning has a standard convergent-series
  criterion. The elementary local norm estimate is inconclusive here.
* Linked-cluster charge elimination needs derived, uniform coefficient and
  tail bounds. No such usable U/t=8 bound was established; an unjustified
  draft criterion was withdrawn rather than reported as a theorem.
* Positive tensor factorizations certify positivity of the proposed matrix,
  but still need a compact positivity certificate for its Hamiltonian residual.
* An identity seed in Hilbert–Schmidt operator space extends the polynomial
  theorem to arbitrary finite Hermitian models. It does not establish small
  MPO bonds or inexpensive residual certification, and it is not implemented
  by this Hubbard-specific checker.

The [primary-source check](PRIOR_ART.md) identifies established positivity,
filtering, inexact-recurrence and renormalization ingredients. No novelty claim
is based on restating them.

## Accounting and reproduction

Host: Apple M1, 8 GiB RAM, macOS 26.4.1. One bounded numerical process ran at a
time, with BLAS threads limited by the existing process wrapper. No paid
compute, library installation, Git write or publication occurred.

The 18 metered subprocesses total **110.9275 seconds**; largest process RSS was
**598,638,592 bytes**. This sum includes unsuccessful screens, an unavailable
optional dependency probe, an initially failed unit-test run, replay,
validation and repeats. It is **not total research wall time**: code authoring,
agent reasoning, literature lookup and small unmetered inspections/tests are
not converted into an invented end-to-end runtime. The failed optional import
was gmpy2; it was not required or installed.

The numerical scripts evolved during investigation. This is not a frozen
held-out campaign. The manifest captures the final files; each run ledger
records its command and measurements. The state-discovery dependency used by
the filter is the newly constructed upper in this directory, not a previous
full-state solution. The old strong lower appears only in the reported
comparison.

From the repository root, the accepted upper replays with:

```sh
.venv-correlated/bin/python -B -m research.constructive_compression_20260916.upper_checker results/constructive_compression_20260916/upper64_exact
```

The target refusal and the loose acceptance replay with:

```sh
.venv-correlated/bin/python -B -m research.constructive_compression_20260916.filter_checker results/constructive_compression_20260916/filter64_exact/proof.json
.venv-correlated/bin/python -B -m research.constructive_compression_20260916.filter_checker results/constructive_compression_20260916/filter64_exact/validation_loose_proof.json
```

The target result is in
`results/constructive_compression_20260916/filter64_exact/proof_checked.json`;
its `accepted` field is false. The loose validation is separately named.
Run ledgers, test output, exact fractions and all unsuccessful proposals are
retained under `results/constructive_compression_20260916/`.

The broad objective remains open. This pass supplies a working verifier and
an exact failure measurement for one new construction, not a completed major
breakthrough.
