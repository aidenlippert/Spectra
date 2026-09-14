# Compact response: a complete H6 control and the remaining compression problem

The response mechanism now participates in an exactly replayed **1.000 mHa
full H6 ground-energy interval**, using the frozen upper. This crosses the
three requested numerical thresholds. The response program itself is compact
and has an operator error certificate that enumerates no states. Completing
the lower proof still expands all 924 six-electron determinants into conserved
blocks. **This is a successful finite control, not a compact complete solver.**

On the separate H6 1.6 Å fixture, the unchanged response rule lowers a
single-determinant Hartree–Fock upper by **20.933188 mHa**, verified exactly.
That is transfer of a useful response component; a tight full interval has
not transferred. The physical-validation track now has a sourced CH2 pilot,
with model and solver error treated separately.

## The frozen full-bound comparison

| Quantity | Width, mHa | Outcome |
|---|---:|---|
| Previously retained restricted-family bound | 3.4154210655 | Improved by 70.721% |
| Original H6 target | 1.6000000000 | Reached |
| Exact old 492-generator family floor, with this upper | 1.5134646590 | Crossed strictly |
| New response plus expanded retained proof | **1.0000000000** | Exact acceptance |
| Previously available enumerated physical reference | At most 0.0549943928 | Still tighter than this control |

The new endpoints, shown approximately, are

    −6.334058626233 ≤ E0(H) ≤ −6.333058626233 hartree.

Their difference is exactly 1/1000 hartree. The upper is unchanged; no upper
improvement is credited to this comparison. These energies concern the frozen
rational Hamiltonian in its entire N=6 sector. They are not laboratory energy
predictions or bounds on physical-model error.

Crossing the old-family floor proves that this accepted lower lies beyond
what that restricted family can deliver at the same upper. It does not prove
a computational advantage over complete enumeration. Earlier project work
already used Schur responses and obtained tighter enumerated bounds. The
present experiment isolates the analytic molecular response construction and
charges its remaining solve. [Exact comparison and scope](/Users/aidenlippert/Documents/Spectra/results/compact_response_20260913/comparison.json)

## What was compressed

Let Q fix double occupancy of the last spatial orbital, P=I−Q, and write
H−b as blocks A, B, B*, D. The certificate proves δI ≤ D ≤ MI directly from
the molecular density factors, one-body terms and certified tail.

The stronger lower endpoint retains the positive contribution from hopping
out of the fixed-occupancy sector. For a spatial density matrix
L=[[A,v],[v*,ell]], exact CAR algebra gives

    Q Q(L)^2 Q = Q [(Q(A)+2ell)^2 + 2||v||² − Q(vv*)] Q.

This improves the old rule that discarded the entire positive density square.
The H6 lower endpoint at b=U−0.001 Ha rises from about 0.190674 to 0.564391 Ha
before outward rounding. It is a bound on this restricted sector, not a
whole-molecule excitation gap. On H8 and the different H6 geometry, the old
rule gives a negative endpoint; this strengthened rule gives a positive one.

| Case | Rounded δ, Ha | Rounded M, Ha | Polynomial degree | Certified response penalty, mHa | Program bytes |
|---|---:|---:|---:|---:|---:|
| Frozen H6 | 0.5643 | 14.6918 | 25 | 0.000520765 | 501 |
| H8 development control | 0.2426 | 39.4789 | 70 | 0.000737819 | 501 |
| H6 1.6 Å response-rule holdout | 0.0779 | 13.6083 | 72 | 0.000837319 | 464 |

The same Chebyshev stopping rule was used throughout: an operator residual
penalty at most 10^-6 Ha, order cap 256, and outward endpoint rounding on a
10^-4 grid. H8 helped develop the endpoint formula, so it is not labeled a
holdout. The 1.6 Å fixture existed in previous work but was not used to fit
this response rule.

The response X=p(D)B* is specified by a polynomial recurrence; no inverse or
expanded operator powers are needed to describe it. Its degree counts D
applications per right-hand side. The recurrence certifies

    E = B* − DX,       ||E||²/δ ≤ η.

It remains necessary to prove

    K = A − BX − X*B* + X*DX ≥ η I.

The program certifies the response error uniformly over retained vectors;
it does **not** certify this retained inequality. A short program also does
not make each action of D cheap. [Derivation](/Users/aidenlippert/Documents/Spectra/research/compact_response_20260913/DERIVATION.md)

## What the complete control still expands

The H6 proof checks every conserved block, verifying spin-up-number and
orbital-parity conservation from the Hamiltonian. P has dimension 714 and Q
dimension 210. The largest full block has 200 states; the largest retained
block has 152 and largest eliminated block 52.

The recurrence is applied to every coupling column. After rounding X to
rational entries, replay recomputes its actual residual exactly; it does not
assume rounding preserved the ideal polynomial error. Rational triangular
factors prove K−ηI positive through an exact diagonally dominant residual.
All 12 blocks pass, including uncoupled sectors.

| Charged object/work, per complete construction | Amount |
|---|---:|
| Physical basis labels generated | 924 |
| Hamiltonian word/state checks | 848,232 |
| Full conserved-block matrix entries | 132,004 |
| Explicit response entries | 23,910 |
| Retained matrix entries | 76,646 |
| Stored retained triangular-factor entries | 38,680 |
| Selected response D actions on individual columns | 17,800 |
| All six degree diagnostics, including the selected degree | 36,312 column actions |
| Expanded closure certificate | 654,983 bytes |
| Shared H6 fixture, tail and frozen upper files | 80,419 bytes |

An exact physical-vector witness also shows that the simpler sufficient test
A−BB*/δ ≥ 0 fails at this target: its expectation is −0.0454303254 Ha. This
refutes that particular scalar bound, not every possible scalar approximation.
Numerical degree diagnostics fail at orders 1, 2, 4 and 8, and pass at 16 and
26; the values are not monotone. Only order 26 was exported as the complete
exact response proof. [Full replay and scalar obstruction](/Users/aidenlippert/Documents/Spectra/results/compact_response_20260913/closure_replay.json)

## Transfer without a new FCI calculation

The fixed sparse trial procedure takes the retained part of a supplied upper
vector, constructs one response direction, optimizes a two-dimensional
Rayleigh quotient, and replays the rounded integer amplitudes exactly.

| Trial input | Final explicit amplitudes | Upper improvement, mHa | Interpretation |
|---|---:|---:|---|
| Frozen H6, inherited 200-amplitude FCI upper | 200 | 0 accepted | Candidate is worse by 0.000000378 mHa; retain incumbent |
| H6 1.6 Å, one Hartree–Fock determinant | 53 | **20.933188** | Exact variational improvement without FCI discovery |

The fresh trial uses 72 D actions, acts on 3,744 source-state occurrences,
and performs 3,436,992 Hamiltonian word/state checks. It stores no Hamiltonian
action cache or sector matrix. Nevertheless, generated destinations reach
all 200 states of the relevant symmetry block, including projected-out
destinations. Its 53 final amplitudes do not represent all the work.

On this fresh fixture, the initial HF-based program target b=U_HF−0.001 Ha
is actually ruled out as a full lower: the accepted trial lies 19.933188 mHa
below it. The positive Q-sector interval and response certificate remain
valid. This is direct evidence that a good response certificate alone cannot
stand in for the retained proof. [Exact trial results](/Users/aidenlippert/Documents/Spectra/results/compact_response_20260913/trial_replay.json)

## Verification and cost

**15 focused tests pass.** They check the compressed-square identity, a
noncommuting exact Chebyshev example, the block congruence, independent CAR
actions, complete block coverage, and rejection of changed responses,
factors, bounds and bindings. The H8 full-sector expansion cap is tested.
All accepting replays run with Python -S and load none of NumPy, SciPy,
CVXPY or PySCF. This is independent fresh-process replay by the supplied
checker, not verification by an outside research group.

A fresh-process H6 benchmark, starting from the frozen input files, took
**2.474 seconds**: 1.417 seconds for response/retained-proof construction and
serialization, then 1.056 seconds for independent exact lower and upper
replay. It includes all six degree diagnostics. Its output certificate
matches the archived one. Measured peak resident memory was 53.36 MiB for
construction and 32.42 MiB for replay. This is one local run, not a scaling
or comparative performance result. [Cold benchmark](/Users/aidenlippert/Documents/Spectra/results/compact_response_20260913/cold_benchmark.json)

The cost ledger records **21.693 seconds of stage times**, including probes,
replays, tests, the benchmark and preservation checks. Nested timings are
not added twice. The full H6 block construction is repeated 13 times across
those paths: 12,012 generated basis labels and 11,027,016 word/state checks.
Sparse trial work and 2,500,632 upper-replay word/state checks are additional.
The stage sum is not calendar research time; literature work, editing and
exploratory reads are outside it. Historical integral, factor and FCI
discovery are explicitly inherited, not eliminated. No new molecular fixture,
SDP or PySCF FCI call was made; the expanded H6 work above remains many-body
computation. [Cost ledger](/Users/aidenlippert/Documents/Spectra/results/compact_response_20260913/cost_ledger.json)

All 744 files in the preceding manifest are unchanged. For a read-only full
replay, run from the project root:

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.compact_response_20260913.verify
```

## What comes next

The original question is still whether a small set of interference patterns
can carry the hard correlations while the remainder is certified collectively.
This pass establishes a compact response description and controls its error.
The 714-dimensional retained proof is the exposed obstacle.

The next mathematical test should replace the expanded retained factors by
a small set of response-coupled directions plus an independently proved
positive remainder, keeping the same H6 target and upper. Success means
preserving the accepted bound while reducing total construction, action and
replay work. Further polynomial tuning without reducing that retained work
would not answer the main compression question.

Alongside it, the [physical-validation pilot](/Users/aidenlippert/Documents/Spectra/research/compact_response_20260913/PHYSICAL_VALIDATION_PILOT.md)
starts with CH2's singlet–triplet gap and a matched coupled-cluster baseline.
It separates total-spin certification, solver intervals and physical-model
corrections. The reference audit corrected a superseded CCl2 measurement;
no chemical prediction or experimental advantage is claimed yet.
