# Direct molecular-control construction: results and unresolved objective

The general compression problem is **not solved**. In particular, this pass
did not construct an enumeration-free certificate for the original H8
four-phase, amplitude-0.5 control task.

It did implement a different, precisely scoped capability: exact observable
certificates constructed directly from Hamiltonian coefficients and a supplied
compact state, with no full-sector trajectory or matrix. This works in a
**much stronger control regime** and transferred to two other molecular inputs.
An exact obstruction now explains why this particular norm-bound mechanism
cannot solve the original amplitude regime by simply improving its norm
estimates or changing its pulse duration.

## Accepted new component

The constructor retains every term of the interacting Hamiltonian. It bounds
its influence collectively through two exact commutator bounds while solving
the two-orbital control algebra exactly. It needs scalar MPS moments, not a
time-dependent many-electron wavefunction.

The target is a dimensionless orbital-population contrast `D <= -0.6`. All rows
include pointwise error 0.001 Ha in each of the two controls, initial trace
distance 0.0005, and integrated phase-flip dephasing rate 0.001.

| Input | Spatial orbital pair, zero-based | Constant v, Ha | Duration, atomic units | Robust D interval, rounded outward |
|---|---:|---:|---:|---:|
| Original H8 | 3, 5 | 64 | 25/1024 | [-1.984548, -0.795337] |
| Asymmetric H6 | 0, 5 | 32 | 25/512 | [-2, -0.992444] |
| Asymmetric water | 0, 6 | 64 | 25/1024 | [-2, -0.696358] |

All rows use `u=0`. H8 retains the source observable and initial state. For the
two transfer inputs, the rule was frozen before execution: select the most and
least occupied spatial orbitals using exact initial MPS contractions, break
ties by smallest index, fix angle 25/8, and double v from 1/2 until the target
passes or a declared cap is reached. Selection is replayed. These are selected
orbital contrasts, not a common measured chemical observable.

The H8 amplitude is **128 times the original 0.5-Ha limit**. The mechanism makes
the pulse sufficiently brief that a coarse but rigorous bound on the
interacting drift remains useful. This is a real finite-model certificate,
but it does not show that Spectra has learned an economical representation of
the original weak-control dynamics. It supplies no calibrated laboratory
actuator, preparation procedure, or physical-model guarantee.

The input state is the specified normalized rational MPS, not an assumed exact
ground state. The constructor/checker uses Python integers and Fractions.
Saved endpoints, success flags, and approximate commutator norms cannot replace
fresh arithmetic. A separate numerical H8 calculation gave D approximately
-1.38927074536, inside the new interval. That reference did enumerate the
4,900-state sector and is charged separately; acceptance does not use it.

## What the obstruction proves

For constant W-only control with `u=0` and `|v|<=0.5`, the global
commutator-norm-sum envelope cannot certify the negative target at **any
duration**, even if its two operator norms are known exactly. Applying the
commutators to one balanced-spin determinant gives exact norm lower witnesses.
Those force the best upper endpoint of this envelope to remain nonnegative.
At angle 25/8, its endpoint is at least 1.4308116314.

This is a proof-family limitation. It is not a lower bound on the physical D,
not a proof of unreachability, and not a limit on the original four-phase
protocol. It identifies the missing information: a successful weak-control
certificate needs to retain useful state-dependent response or correlation
information instead of bounding all drift through these two global norms.

The complete inequalities, exact arithmetic construction, robustness argument,
and a correctly signed prospective adjoint identity are in DERIVATION.md.

## Tensor construction and diagnostic results

An exact operator constructor now represents the supplied 2,912 CAR terms as
a sparse matrix product operator with 3,916 symbolic local connections and
maximum bond 154. It re-expands exactly to the original coefficient dictionary.
This takes 0.4499 seconds internally and never constructs a configuration
basis. It is a reusable operator component, not a completed dynamics solver.

The direct tensor-state path did not pass the control target. Its first
bond-64 Hamiltonian-action error bound was about 1.34 million. Canonicalizing
the operator and state reduced the bound to 192.9873, still unusable. A loose
sufficient bound does not establish that the best approximation is bad.

The supplied H8 MPS already has maximum bond 144 and 10,093 nonzero entries.
Checked Schmidt-tail diagnostics constrain how aggressively that particular
state can be compressed:

| Maximum MPS bond | Canonical-basis state-error lower bound | Declared local-basis state-error lower bound |
|---:|---:|---:|
| 32 | 0.11217077 | 0.00899726 |
| 64 | 0.05660582 | 0.00522857 |
| 96 | 0.02734849 | 0.00320666 |
| 128 | 0.01112959 | 0.00174357 |

These rounded-down lower bounds concern vector norm and the declared orbital
basis/order. They are not lower bounds on an observable error. For the
diagnostic tolerance 0.005, they rule out bond 128 in the original order and
bond 64 in the tested local basis. They do not prove the next bond will suffice
for dynamics.

The local rotation used 28 exactly orthogonal rational Givens rotations, with
bounded numerical application error about 1.8851e-9. Its one-particle proposal
was inherited from prior work. A fresh PySCF proposal attempt timed out at 90
seconds; the cause was not established. No claim of fresh orbital discovery is
made. Numerical tensor bounds assume the explicit IEEE binary64/BLAS rounding
model. The exact strong-control accepting path does not depend on them.

## Preserved reference and input identity

A fresh replay of the uploaded archive reproduced the accepted 24-coordinate
short and 64-coordinate long protocols, and refused the insufficient
eight-coordinate proof. Its 20 original tests passed. It took 47.8699 seconds
and peaked at 173,621,248 bytes. This remains an enumerated replay: 4,900
balanced-spin configurations and 903,620 Hamiltonian entries.

The uploaded Hamiltonian is in canonical RHF orbitals. It is exactly
phase-equivalent to the older local H8 fixture via spatial signs
`[1,-1,1,1,1,-1,1,1]`. This leaves D unchanged and flips W. A complete driven
equivalence would additionally require the corresponding state and waveform
transformation. No source Hamiltonian or MPS was edited.

## Measurement and verification

These timings are local component measurements on the same M1 Mac with 8 GB
memory. They are not matched-protocol speedups:

| Work | Wall time | Peak process memory |
|---|---:|---:|
| Original enumerated archive replay, including original tests | 47.8699 s | 173.6 MB |
| Exact sparse operator construction and symbolic replay | 0.5373 s | 32.2 MB |
| New H8 strong-control construction process | 1.7202 s | 41.7 MB |
| Fresh H8 accepted arithmetic in final replay | 1.7104 s | Included below |
| Final H8 replay with low-amplitude refusal, obstruction and four corruptions | 9.9969 s | 40.9 MB |
| H6 and water construction, selection and fresh replay together | 4.0629 s | 36.6 MB |
| Separate enumerated numerical check of the new H8 pulse | 6.7665 s | 291.0 MB |

The final focused suite passed 17 tests. They cover independent small-system
fermion actions, exact operator reconstruction, control identities,
trigonometric enclosures, generalized input/refusal behavior, tensor residual
enclosures, orbital gates and Schmidt-tail bounds. Fresh replay rejects changed
endpoints, amplitudes, uncertainty budgets and an invalid zero-quadrature
premise. There is no formal theorem-prover verification or claim of bug-free
software.

ACCOUNTING.json in the result directory lists every metered subprocess,
including unsuccessful attempts and repeated checks. The cumulative subprocess
time is not a cold complete solver cost. It excludes inherited integrals/MPS
discovery, prior orbital-proposal discovery, editing, reasoning, and some early
unmetered subagent diagnostics. Peak memory is a maximum across measured
processes, not a sum. The largest recorded tensor diagnostic used about
607.1 MB. No paid compute or library installation was used; nothing was pushed
to GitHub.

Some early subagent outputs contained encoding, coefficient-scaling or
square-root-scale mistakes. Those implementations were replaced or removed,
the affected calculations rerun, and their outputs explicitly excluded in
artifact_disposition.json. The corrected results are the only ones above.

## Reproduction and remaining objective

The standalone exact component bundle includes the literal Hamiltonians and
MPS inputs, current checker sources, integrity hashes and fresh-replay script.
It includes no original reduced embeddings or teacher trajectories. From its
root, run:

```text
python3 -B -S replay.py --out ../fresh-direct-control-replay
```

The output directory must be new and outside the bundle. This reconstructs
the three strong-control certificates, the weak-control envelope obstruction,
the exact operator representation and refusal tests using only the standard
library. Repository-only tensor diagnostics additionally need the existing
numeric environment and retain their separate arithmetic assumptions.

The outstanding scientific result is still a useful certificate for the
original weak-control dynamics whose discovery and acceptance both avoid
full-state dependence at manageable cost. Nothing here closes general many-body
scaling, physical actuator mapping, state preparation, or synthesis. The
strong-control component and proof-family obstruction are retained advances;
they must not be substituted for that unresolved result.
