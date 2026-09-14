# Coupling transfers, but this compact rule misses most of the useful structure

**The first reusable contraction rule failed to recover most of the fitted
H6 proof's improvement.** It uses 108 generated combinations and gives a
12.740724 mHa interval, versus 13.010569 for the quadratic baseline. The
previous fitted 216-combination proof remains at **3.415421 mHa** and is
frozen. The new rule preserves only **2.81% of that fitted bound gain**, despite
capturing 56.75% of the proof's aggregate squared operator size.

The same frozen rule improves the bound on both new H6 geometries and the
larger H8 fixture. Those improvements are modest, and every interval remains
above 1.6 mHa. This supports the usefulness of coupling across cases; it does
not establish a compact solution of the difficult correlations.

All twelve new energy certificates passed independent exact replay. No
additional adaptive H6 enrichment was run in this pass.

![The fixed contraction rule compared with quadratic and separate controls](/Users/aidenlippert/Documents/Spectra/results/mechanism_transfer_20260913/transfer.png)

## First, diagnose the frozen 3.415 mHa interval

The preceding comparison and its full-family witness are complete. Its
independent physical audit distinguishes two exact limitations on the
entire completed 492-generator family:

| Quantity | Certified value, mHa |
|---|---:|
| Best found interval | 3.415421 |
| Actual error of its lower bound | between 3.360427 and 3.415421 |
| Minimum possible interval with the frozen upper, from the full-family witness | 1.513465 |
| Minimum lower-bound error relative to the true ground energy, from that witness | 1.458470 |
| Distance from the found lower to the accepted family ceiling | 1.901956 |

The family ceiling leaves the 1.6 mHa target unresolved. Its 1.901956 mHa gap
above the found lower cannot be assigned wholly to missed directions: the
ceiling itself may be loose, and no full-family optimum has been established.

The frozen upper's error is at most **0.054994 mHa**. The full collective
representation envelope is **0.108331 mHa**, although the explicit lower-tail
shift paid by this certificate is only 0.000003632 mHa. Its exact residual
penalty is **0.009140 mHa**. Even allowing the upper and representation terms
their full ranges, at least **3.242955 mHa** remains between the retained
ground energy and the certificate's scalar before residual correction.
The dominant remaining gap is therefore in the lower construction.

These allowances are not counted twice. In particular, the additional
0.006569 mHa loss relative to the reported floating objective is a different
view of export accuracy, not another penalty to add to the exact residual.
The [exact diagnosis](/Users/aidenlippert/Documents/Spectra/results/mechanism_transfer_20260913/interval_diagnosis.json)
records the algebra and its limits. The
[frozen comparison report](/Users/aidenlippert/Documents/Spectra/research/trace_pricing_20260913/REPORT.md)
records the accepted family witness and separate physical-reference replay.

## What structure the successful operators contain

The analysis used canonical operator coefficients and the uniform fixed-N
trace, avoiding the poor conditioning of the original frame coordinates.
The accepted proof contains 198 stored cubic/linear anticommutator factor
rows. Across its eight blocks, **30 singular directions account for 99%**
of their squared operator size, and **36 account for 99.9%**. These are
approximate factor-spectrum measurements, not a 30- or 36-direction energy
certificate.

Four spin-flipped block pairs have relative operator-Gram differences
between **0.10% and 0.39%**. Thus spin pairing is a recurring feature. The
measured differences are nonzero; no exact identity between fitted blocks
is inferred from that near agreement.

We tested the shared generation rule

    B[i,s,t,u; alpha] = sum[k,p] lambda[k]^alpha L[k][i,p] Q[k,s,t] a[p,u].

It contracts the pattern and annihilation indices through the Hamiltonian's
own spatial matrices, retaining all four density spin components. Linear
annihilation operators are included for joint optimization. Three individual
exponents and their three pairs were inspected on the training proof only.
The pair **alpha = 1/2 and 1** had the greatest trace projection coverage and
was frozen before new fixture generation or energy solves.

This gives **18s generated operators** for s spatial orbitals in the tested
symmetry decomposition: 108 for H6 and 144 for H8. No numerical duplicate
removal occurred in any tested case. Exact modular rank checks establish
their full independence in canonical operator coefficient space.

There is an exact identity behind part of this rule. The spin-summed alpha=1
contraction is the cubic term in [Hret,a_i], with known linear corrections.
Exact CAR subtraction verified this for all twelve training modes. The
alpha=1/2 weighting and independent spin components are candidate extensions.
The [proof](/Users/aidenlippert/Documents/Spectra/research/mechanism_transfer_20260913/PROOF.md)
distinguishes that identity, rounded numerical directions and accepted energy
certificates. No novelty claim is made for the equation-of-motion identity.

## Transfer under the frozen rule

| Case | Role | Generated directions | Quadratic interval | Coupled rule | Same rule, separate |
|---|---|---:|---:|---:|---:|
| H6, 1.4 Å | Training | 108 | 13.010569 | 12.740724 | 13.010536 |
| H6, 1.2 Å | New geometry | 108 | 10.756607 | 10.535842 | 10.756712 |
| H6, 1.8 Å | New geometry | 108 | 13.000525 | 11.685282 | 13.000524 |
| H8, 1.4 Å | Larger, held out from rule selection | 144 | 26.895633 | 26.676823 | 26.895812 |

All intervals are in mHa. Each row uses the same rational Hamiltonian,
fixed-particle sector, reference upper and tail across its three controls.
H8 already had project results; it was held out from selection of this rule,
not claimed to be an entirely untouched molecule. The two changed H6
geometries were newly generated after the rule was frozen.

The coupled gains over the baseline are **0.269845, 0.220765, 1.315243 and
0.218810 mHa**, respectively. A few separate-control results are slightly
weaker than the baseline because these are finite numerical solves followed
by exact correction; earlier stronger bounds remain retained. No optimum
claim follows from these controls.

The spatial-pattern count follows the fixed 2s-2 rule: ten for H6 and
fourteen for H8. Their independently checked tail widths are respectively
0.108331, 0.540096, 0.005720 and 0.187166 mHa in table order. The original
training Hamiltonian and tail were preserved.

## What the rule actually costs

| Coupled result | H6 training | H6, 1.2 Å | H6, 1.8 Å | H8 |
|---|---:|---:|---:|---:|
| Added block dimensions | four 6, four 21 | same | same | four 8, four 28 |
| Gram entries including baseline | 7,632 | 7,632 | 7,632 | 21,504 |
| Compact certificate, bytes | 103,457 | 104,120 | 99,115 | 252,404 |
| Expanded certificate, bytes | 895,899 | 896,765 | 849,727 | 2,771,398 |
| Stored factor rows, including baseline | 309 | 309 | 300 | 505 |
| Nonzero direction coefficients | 2,892 | 2,892 | 2,892 | 7,184 |
| Nonzero compact factor coefficients | 5,748 | 5,748 | 5,519 | 15,895 |
| Whole case discovery, all three controls, seconds | 21.540 | 22.028 | 21.977 | 131.477 |

Factors and directions use denominators 10^10. The largest factor integers
use 32--33 bits; the largest direction integers use 34 bits. Exact rank
records and per-solve construction, solve, export, residual and verification
costs are in the ledger. Factor rows and approximate trace-effective rank
are different quantities.

The implementation still prepares the complete 492-generator H6 or
912-generator H8 frame before inserting the compact rule. This costs about
4.2 seconds for H6 and 23.3 seconds for H8 and is fully charged. The linear
formula count therefore does not establish linear total complexity. For
H8, the coupled trial alone takes 52.129 seconds, including 18.809 solving
and **30.922 exact acceptance**; verification is already a substantial cost.

All cases completed inside their 240-second budgets. Total case discovery
cost was **197.022 seconds**, including the failed setup attempts. Training
structure analysis added 4.436 seconds. Including discovery of the fitted
training proof gives **743.196 seconds of cumulative causal discovery**.
Earlier full-family diagnostic and audit costs remain in the frozen branch's
ledger; this figure is not the entire project's elapsed time.

The new geometry inputs initially failed because the SDP Python lacked
PySCF. They were generated with the existing molecule environment, with
the failed attempts charged against the original budgets. No successful
scientific case was rerun and no rule was retuned. Fixture generation,
including those failures and the reference work, cost **6.539 seconds**.

Numerical FCI on the two new H6 cases used a 400-determinant spin sector and
took **0.024 and 0.021 seconds** for the eigensolves, after integral/RHF
setup. Their rounded 200-amplitude vectors give independently checked
rational upper bounds. FCI is much cheaper as a numerical energy comparator
on these small cases; its floating energies are not presented as rigorous
ground-energy lower bounds. H8 uses its previous 1,000-amplitude reference.
Upper-reference discovery remains many-body work, separate from the lower
mechanism.

## What this changes about the next question

Operator-size concentration and shared spin structure are real observations
in the fitted proof, but this experiment shows that they do not identify
the energy-critical cancellations well enough. A rule with 56.75% trace
coverage retained only 2.81% of the fitted lower-bound gain. Treating trace
coverage as an accuracy proxy would have produced a misleading success.

The useful next target is a structured factorization that preserves those
coefficient cancellations while sharing spin and orbital factors. It should
be evaluated by exact retained bound gain before being treated as a mechanism.
These tested geometries are now validation data, not untouched targets for
further rule selection. Another compact candidate needs a fresh transfer
test. The current results justify continuing mechanism discovery, not
claiming either a general solution or a proved failure of the full family.

The [direct source comparison](/Users/aidenlippert/Documents/Spectra/research/mechanism_transfer_20260913/RELATED_WORK.md)
maps the actual code to the requested spin-adapted SOS/v2RDM construction
and established column-generation methods. It identifies the restricted
cubic anticommutators and exact residual/tail handling without claiming
that direction selection or SOS itself is new.

## Verification and artifacts

Three new checker tests passed, including an independently known sixteen-mode
example, malformed input refusals and exact reproduction of the old H6 lower.
Two additional exact-rank tests cover rational independence, dependence and
refusal of noninvertible modular denominators.
The exact equation-of-motion identity passed for all twelve training modes.
**All twelve new intervals passed fresh standard-library replay in 80.457
seconds**, with no NumPy, SciPy, CVXPY or PySCF loaded on that accepting path.
The short identity check overlapped this replay after discovery. All **528
files** in the frozen trace-pricing manifest remain unchanged.

The [protocol](/Users/aidenlippert/Documents/Spectra/research/mechanism_transfer_20260913/PROTOCOL.md),
[frozen rule](/Users/aidenlippert/Documents/Spectra/results/mechanism_transfer_20260913/frozen_rule.json),
[training analysis](/Users/aidenlippert/Documents/Spectra/results/mechanism_transfer_20260913/training_analysis.json),
[fresh replays](/Users/aidenlippert/Documents/Spectra/results/mechanism_transfer_20260913/fresh_replay.json)
and [cost ledger](/Users/aidenlippert/Documents/Spectra/results/mechanism_transfer_20260913/cost_ledger.json)
preserve the measurements, failures and hashes. The PNG and SVG comparison
figures are both available in the results directory.

From the workspace root:

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m unittest \
  research.mechanism_transfer_20260913.test_core -v
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m \
  research.mechanism_transfer_20260913.replay
```

Discovery entry points refuse existing output directories. The fixture
builder uses `.venv-molecule/bin/python`; SDP discovery uses the numerical
environment recorded in `environment.json`. The rule and completed results
must remain frozen when reusing this branch for future investigations.
