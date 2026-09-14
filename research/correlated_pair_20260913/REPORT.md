# Correlated states and two-sided certificates

Completed local pass, September 14, 2026. This report supersedes the preliminary
branch reports in this campaign. Existing campaigns were preserved.

**The correlated-upper milestone is achieved on H6 and corrected H8. The
inexpensive correlated lower-proof milestone is not.** Direct Hamiltonian
optimization produced charge-preserving rational MPS witnesses within the
requested 0.500 mHa upper budget on both cases. Their complete intervals use
the inherited full-cubic lower proofs.

| Exact rational model | Lower, Ha | New upper, Ha | Full width, mHa | Upper deterioration from accepted reference, mHa |
|---|---:|---:|---:|---:|
| Original H6 | −6.333113620626 | −6.332955060480 | **0.158560** | **0.103566** |
| Corrected H8 | −9.256055014194 | −9.254878154398 | **1.176860** | **0.426275** |

Decimals are summaries; the linked receipts contain exact rational endpoints.
The obsolete H8 reference and target were not used.

- [H6 complete replay](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/mps/h6_b48_real/interval.json)
- [H8 complete replay](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/mps/h8_spatial_warm144/interval.json)

## What this does and does not establish

The accepting upper calculation no longer needs an FCI vector. The MPS is
optimized from the Hamiltonian, rationalized, and checked by a separate
standard-library contraction implementation. Neither its discovery nor its
upper verification builds a full determinant Hamiltonian or exports the full
wavefunction. The inherited lower already avoided full fixed-N enumeration;
its expensive Gram discovery remains a dependency.

**This is not yet evidence that the difficult correlations have been compressed
into a smaller representation than the reference calculation.** The retained
MPS has more nonzero tensor entries than the enumerated reference has nonzero
amplitudes on these small fixtures. The method supplies a usable factored
upper, not a demonstrated complexity advantage.

| Retained upper | H6 | H8 |
|---|---:|---:|
| Maximum exported spin-site bond | 48 | 144 |
| Nonzero rational tensor entries | 1,012 | 10,093 |
| Reference witness nonzero amplitudes | 200 | 2,468 |
| Equal-spin determinant space dimension | 400 | 4,900 |
| State artifact bytes | 20,623 | 206,857 |
| Largest sparse transfer environment, entries | 296 | 2,098 |
| Integer multiplications in upper replay | 2,034,689 | 143,247,748 |
| Exact upper replay, seconds | 0.592 | 32.951 |
| Inherited lower replay in the same run, seconds | 6.669 | 35.184 |

The H8 successful discovery path used a bond-96 spatial MPS, then its own
bond-128 and bond-144 continuations. Those three recorded runs cost **414.36
process-wall seconds**, including their construction and export. The last
94.46-second internal solve alone is not its discovery cost. The winning H6
run took 8.55 internal seconds; unsuccessful and alternative experiments are
charged separately in the campaign ledger.

## Upper construction and accepting identity

The proposer uses established [Quimb MPO/DMRG machinery](https://quimb.readthedocs.io/en/latest/operator/operator-basics.html).
Quimb 1.15.0 was installed in a separate `.venv-correlated` environment; the
existing environments were not changed. Numerical versions are recorded in
the discovery receipts. Real CAR/Jordan–Wigner MPOs, numerical number
penalties, and two-site sweeps optimize the state. H8 benefits from grouping
the two spin sites of each spatial orbital and compressing the numerical MPO.
That compression is only part of proposal generation.

The numerical optimizer is penalty constrained, rather than a fully U(1)
symmetric DMRG implementation. Export establishes exact charge conservation:
left-canonical bond charge matrices propose a gauge, forbidden entries are
zeroed, remaining entries are rounded to integers over 10^9, and unreachable
bond vertices are removed. Every retained tensor edge must advance the
declared alpha/beta occupation labels correctly. Scalar endpoint charges
enforce the full particle count and spin projection.

The accepting path contracts the **actual rounded state**, including its norm,
against the original exact rational Hamiltonian:

\[
U=\frac{\langle\phi|H|\phi\rangle}{\langle\phi|\phi\rangle}.
\]

It uses integer transfer matrices and shared prefixes of local CAR-word
products. Site denominators cancel in the quotient. Numerical optimizer
energies, MPO truncation estimates, gauge eigenvalues, and discarded tensor
sizes are not accepted as error bounds. H6/H8 MPS witnesses certify particle
number and spin projection; no pure-total-spin claim is made for them.

Source: [mps_exact.py](/Users/aidenlippert/Documents/Spectra/research/correlated_pair_20260913/mps_exact.py), [mps_round.py](/Users/aidenlippert/Documents/Spectra/research/correlated_pair_20260913/mps_round.py),
[mps_direct.py](/Users/aidenlippert/Documents/Spectra/research/correlated_pair_20260913/mps_direct.py), [mps_spatial.py](/Users/aidenlippert/Documents/Spectra/research/correlated_pair_20260913/mps_spatial.py).

## Structured lower: implemented, but not competitive

The main experiment implements a fixed-guide adaptation of
[Hastings’s 2024 paired self-consistent SOS construction](https://arxiv.org/abs/2412.03564).
It is not a reproduction of the paper’s generalized rotated-overlap algorithm.
The guide is the HF particle-hole transform with a rational chemical potential
inside the positive-weight interval. The weighted positive factors are

\[
\sum_i w_i\big[(d_i+\tau_i)^\dagger(d_i+\tau_i)
                       +\tau_i\tau_i^\dagger\big],\qquad w_i>0.
\]

The creator-distribution map exactly reproduces the nonscalar perturbation in
the cross terms before rounding. Cubic anticommutators cancel degree six by
CAR algebra. Five damped coefficient updates retain all induced terms, track
the scalar separately, and round only proposal coefficients. Acceptance
reconstructs the actual squares plus sector chemical shift, checks charge and
positivity, and bounds the entire remaining polynomial by its coefficient L1
norm. No expectation is substituted for an operator identity.

On H6 the lower improved from −60.151593 Ha for the guide control to
**−6.984202 Ha** after five updates. It uses **24 weighted positive factors and
1,584 correction coefficients**, with a 0.591563 Ha residual allowance. Paired
with the accurate MPS its interval is **651.246 mHa**. Discovery took 90.52
recorded seconds; the selected lower replay took 9.19 seconds. Fewer factors
did not give either the required accuracy or cheaper replay than the retained
cubic certificate. The H8 two-update attempt timed out at 210 seconds and
produced no accepted final structured certificate.

The correlated state then directly taught the factors. Exact MPS moments
proposed each correction’s scale by minimizing its positive-factor
expectation. Four blends were accepted and ranked by the actual full bound:

| Blend toward state-optimal scales | Positive-factor expectation, Ha | Full interval, mHa |
|---|---:|---:|
| 0 | 0.067122 | **651.246** |
| 1/4 | 0.063721 | 1,265.765 |
| 1/2 | 0.061292 | 2,183.080 |
| 1 | 0.059349 | 4,095.005 |

Better annihilation of the correlated state made Hamiltonian reconstruction
worse. Exact selection rejected those changes. Each row also checked
\(G\le U-L\le G+2\varepsilon\), where G is the summed positive-factor
expectation. This is an actual state/proof coupling experiment, with a negative
result for this particular rescaling rule. It does not rule out richer
correlated corrections.

Source and evidence: [fixed_guide.py](/Users/aidenlippert/Documents/Spectra/research/correlated_pair_20260913/self_consistent/fixed_guide.py),
[H6 iteration history](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/self_consistent/h6_fixed/receipt.json),
[state-guided experiment](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/state_guided_h6/receipt.json).

## Density conditioning and the complement branch

The [June 2026 density-interaction extension](https://arxiv.org/abs/2606.31765)
motivates a separate exact test. The accepted implementation uses actual CAR
products, a paired density-conditioned correction, and a four-mode resonant
example. An allowed six-mode example leaves five sixth-degree terms; deleting
them would require a **23/147 Ha** coefficient-norm allowance. Quartic-only
closure is therefore not assumed.

The restricted attractive-density parent allocation fails on the molecular
PH guides. An exact max-flow/min-cut calculation yields the following simple
necessary subset inequalities, using all modes:

| Case | Attractive-edge allocation required, Ha | Available positive guide capacity, Ha | Exact deficit, Ha |
|---|---:|---:|---:|
| H6 | 9.786704277212 | 4.145619771922 | 5.641084505290 |
| H8 | 14.666468900234 | 5.584000342660 | 9.082468557574 |

These certify failure of the specified allocation family **after dropping
positive density interactions**. They do not show that the physical
Hamiltonian, the full density operator, or a correlated/cluster extension is
infeasible. The toy resonant construction is verified; no useful molecular
resonant-cluster lower was obtained.

For the retained H6 MPS, the exact residual calculation gives

\[
\|(H-U)\psi\|^2=0.0001639613684067571\ \mathrm{Ha}^2.
\]

Construction and contraction of H² cost 38.22 seconds and generated 49,405
operator terms. No determinant sector was built. Each adjacent two-spatial-
orbital reduced density matrix is exactly positive definite, checked in
charge blocks of dimension at most four. Thus no nonzero operator supported
on one of those clusters exactly annihilates this rounded MPS. This obstructs
that exact local-parent construction only. A bound on the **entire orthogonal
complement** above U was not established, so the residual is not converted
into a ground-energy lower bound.

Evidence: [density closure and cuts](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/density_conditioned/exact_diagnostic.json),
[correlated residual and local-parent checks](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/complement_h6.json).

## Response, transfer, and CH₂

At matched complete widths of 0.160560 mHa and 1.178860 mHa, a response-assisted
control reused the same cubic lower and inherited sector-gap proofs as the
direct route. Response reconstruction and replay added **2.15 seconds on H6
and 4.90 seconds on H8**, with no accuracy or discovery gain. A competitive
response-assisted structured terminal proof was not obtained. The earlier
operator-action savings are preserved; this comparison supplies no solver
speedup claim. See [matched response control](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/response_comparison.json).

The upper and direct structured-lower rules were frozen before generating
**H6 at 1.91 Å**, absent from the existing fixture inventory. Generation used
ordinary integrals/RHF and one HF determinant check, with FCI disabled. Both
endpoints were then discovered from that Hamiltonian without a cubic teacher:

- Lower: −6.642889888550 Ha.
- Upper: approximately −5.272650981672 Ha.
- Complete width: **1,370.238907 mHa**; the tight-interval target did not transfer.
- Source hashes remained unchanged across generation and final acceptance.

This is a full-interval transfer result, not just an upper improvement. Its
width exposes the weakness of the independent lower construction. See the
[freeze](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/transfer_freeze.json) and
[new-geometry interval](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/fresh_h6_1p91/interval.json).
H10 was not attempted because the H8 accuracy and cost results do not justify
that expansion in this pass.

For the declared frozen-core, fixed-geometry CH₂ rational model, new correlated
MPS states plus the inherited enumerated spin-sector lower certificates give

\[
E_S-E_T\in[40.6043047894,\ 40.6247594413]\ \mathrm{mHa}.
\]

The triplet is therefore lower within this model; the difference interval is
**0.02045465 mHa** wide. Exact S² expectations and the exact commutator
[H,S²]=0 certify a nonzero implicit projection into the requested total-spin
sector. A global lower bound controls the projection’s energy allowance.
The rounded unprojected tensors are not claimed to be exactly spin pure.
The projected upper allowances are below 1.6×10⁻⁹ mHa. Numerical state
construction took about 4.2 seconds per spin; complete exact model-gap replay
took 18.75 recorded seconds. Lower verification still enumerates 400/225
spin-projection labels. No physical-model error, zero-point correction,
experiment comparison, or experimental advantage is certified.

Evidence: [CH₂ spin interval](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/ch2/interval.json).

## Accounting, validation, and preservation

The ledger records **2,032.69 process-wall seconds** for completed solver,
installation, experiment, and test runs, including failures and timeouts,
before preservation scans. The largest logged single-job peak was
974,766,080 bytes; simultaneous aggregate memory was not measured. Jobs used
one BLAS/Numba thread. Historical cubic discovery costs of 259.309/721.563
seconds are separate inherited costs, not included in that new-work sum.

This is a lower bound on total campaign computation: two early failed
optimizer attempts have missing timing receipts, and two early successful
smoke summaries preserve only partial internal timing. Consequently the
prospective accounting objective was not fully satisfied; the initial
2,400-second budget is not claimed as a completely audited end-to-end total.
All final construction and acceptance runs have durable budget receipts.

Final campaign collection: **11 passed, 3 explicitly skipped prototype modules**, with ten mutation subtests passing. The eleven accepted tests cover: exhaustive tiny CAR words through degree
four against an independent occupation-bit oracle; repeated-mode signs;
energy/norm contractions; input and charge mutations; unequal-weight
creator identities and literal paired products; PH involution; exact density
cuts; retained higher-degree closure; positive-definiteness checks; and a
Quimb grouped-MPO/charge-export comparison against the independent exact
oracle. An initially wrong test mutation used a valid hole charge; its
correction and both test logs are retained. Exact H6/H8 and CH₂ accepting runs
used Python `-B -S` without numerical packages.

The preservation audit found **zero changes among 9,643 inherited files**
relative to this campaign’s starting inventory. One analysis-only JSON had
already differed from the parent campaign’s seal at the start; both hashes
are recorded. The parent manifest was not rewritten to conceal that mismatch.
See [audit.json](/Users/aidenlippert/Documents/Spectra/results/correlated_pair_20260913/audit.json).

Early incorrect implementations were rejected, retained as history, and made
to refuse execution. Their preliminary reports do not support any claim here;
see [REJECTED_PROTOTYPES.md](/Users/aidenlippert/Documents/Spectra/research/correlated_pair_20260913/REJECTED_PROTOTYPES.md).

## Replaying the main results

From the repository root, with a local Python interpreter:

```sh
python -B -S -m research.correlated_pair_20260913.intervals h6 results/correlated_pair_20260913/mps/h6_b48_real
python -B -S -m research.correlated_pair_20260913.intervals h8 results/correlated_pair_20260913/mps/h8_spatial_warm144
python -B -S -m research.correlated_pair_20260913.ch2
python -B -S -m research.correlated_pair_20260913.finish_transfer
```

These commands recompute and overwrite the corresponding receipt files; use
fresh output copies or the lower-level check functions when preserving the
sealed run’s timing receipts. Reference-vector replay in `intervals` is a
separately marked comparison and is unnecessary for validity of the new
upper-plus-inherited-lower interval.

The next mathematical bottleneck is a correction family that keeps both
state annihilation and Hamiltonian reconstruction accurate. This pass shows
why optimizing either one alone is insufficient. The sub-1.6 mHa milestone
survives with independently constructed correlated uppers; a small, inexpensive
state-and-proof pair remains unachieved.
