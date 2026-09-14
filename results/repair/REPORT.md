# A constructed representation supplies the missing prediction input

Spectra now constructs an executable predictive state from observable records of
the existing thermal workload. It requests additional observations when the old
records do not identify the future response, derives the state update without
receiving the thermal equations, and supplies the resulting predictions to the
unchanged policy bridge. A previously unjustified operating policy becomes
certifiably feasible within the admitted exact linear model class.

This removes a supplied-state-and-update assumption for this bounded family.
The construction uses established realization theory. It is neither a new
scientific method nor evidence of compounding, hardware validation, or general
mastery of matter. The mission and constructive bridge remain unchanged.

## The workload and the two different missing distinctions

The existing equations in [mission_examples.py](../../experiments/mission_examples.py)
are unchanged:

\[
\begin{bmatrix}h_{t+1}\\c_{t+1}\end{bmatrix}
=\begin{bmatrix}1/2&1/4\\1/4&1/2\end{bmatrix}
\begin{bmatrix}h_t\\c_t\end{bmatrix}
+\begin{bmatrix}1/2\\0\end{bmatrix}u_t.
\]

The reset is zero, controls satisfy \(0\le u_t\le1\), and the original terminal
requirements are \(c\ge1/10\) and \(h\le3/10\). These are the workload's
mathematical quantities, not calibrated temperatures or energy measurements.
This pass explicitly extends the requested horizon to two preparation steps
followed by two operating steps. The old two-step records are preserved.

**Discarded information.** Preparations \(a=(1,1/4)\) and \(b=(0,3/4)\)
both give current hot output \(h=3/8\), but cold outputs \(1/8\) and \(0\).
Their next hot outputs under zero input differ by \(1/32\). A predictor using
only current hot output must therefore incur error at least \(1/64\) on one
of these histories. The cold output is already a declared observable port.
The hot-only representation is a diagnostic ablation introduced for this pass;
the original equation evaluator already tracked both temperatures. This is not
evidence of a previously hidden physical sensor or a newly discovered law.

**Insufficient evidence about the update.** All nine original two-step records
identify only the first two Markov vectors,
\(g_0=(1/2,0)^\top\) and \(g_1=(1/4,1/8)^\top\). Retaining both current
outputs does not make the longer input/output law identifiable from those
records. For example,

\[
A'=\begin{bmatrix}1/2&1/4\\1/4&1/4\end{bmatrix},\quad
B'=\begin{bmatrix}1/2\\0\end{bmatrix},\quad C'=I
\]

reproduces every original record, yet disagrees on the consequential policy
below. This witness was constructed for the independent analysis; the learner
was not given a menu containing either model. The learner's automatic diagnostic
is a rank test for identifying Markov vectors, not an autonomous search for
policy-specific countermodels or a minimal experiment.

## What the executable learner actually constructs

The learner receives only a JSON file containing the declared model prior,
control records, reset declarations, and exact intervals for two output ports.
It receives neither physical state matrices nor an intended rank. The adapter
extracts the ports from the old task-output interface and discards the evaluator's
internal trajectory. A separate learner process imports no thermal evaluator.

The prior is one-input, two-output, zero-reset, strictly proper, deterministic
discrete-time LTI behavior of order at most \(D=4\). These are supplied
assumptions. They are not inferred by an apparent rank plateau.

On the old records the learner returns `needs_evidence`. It requests an input
impulse followed by seven zeros, observing both ports after every step. After
those exact equation evaluations, it constructs a rank-two Hankel realization
and selects these predictive coordinates:

\[
z_1=h_t,\qquad z_2=\text{hot output after one zero-input step}.
\]

The discovered matrices are

\[
z_{t+1}=
\begin{bmatrix}0&1\\-3/16&1\end{bmatrix}z_t
+\begin{bmatrix}1/2\\1/4\end{bmatrix}u_t,
\qquad
\begin{bmatrix}h_t\\c_t\end{bmatrix}
=\begin{bmatrix}1&0\\-2&4\end{bmatrix}z_t.
\]

Thus the additional predictive statistic encodes the cold-side distinction.
It was selected algebraically from observable response data, not supplied with
its physical name and update law. The method and ordered candidate test family
are supplied conventional mathematics; the selected coordinates and matrices
are data-derived outputs. This is a reconstruction of two required coordinates,
not a reduction below the physical system's dimension.

Initialization is executable: start at \(z=0\) after a confirmed zero reset and
update through the recorded preparation controls. The current implementation
refuses unknown reset/history and uncertain measurement intervals. It does not
provide observation-based reinitialization or a noisy-state estimator.

## Why the guarantee extends past the acquisition records

Let \(G\) contain \(g_0,\ldots,g_{2D-1}\). For every reset record the
convolution equation is \(XG=Y\). The independent checker requires
\(\operatorname{rank}X=2D\), verifies this equation against the raw observable
records, and checks all candidate identities

\[
C_R M^k b=g_k,\qquad 0\le k<2D.
\]

The difference between the true and candidate impulse responses has realization
order at most \(n+r\le2D\). Cayley–Hamilton makes every later difference zero
once these first \(2D\) vectors vanish. Therefore every finite reset input
history has the same outputs, conditional on the declared order and LTI prior.
This is an exact theorem, not a statistical extrapolation from the word sweep.
See the [construction and proof](../../research/repair/realization_proof.md).

Independently, the supplied-equation checker derives

\[
\Pi=\begin{bmatrix}1&0\\1/2&1/4\end{bmatrix}
\]

and verifies \(\Pi A=M\Pi\), \(\Pi B=b\), \(C_R\Pi=I\), and reset
compatibility. Physical matrices occur only in this separate checker. These
identities also establish the predictive interpretation of the selected tests.

The result covers causal, well-defined feedback using preserved ports and
producing controls in \([0,1]\), while the component dynamics remain unchanged.
Output equality implies identical feedback inputs, so the guarantee composes
along those histories. A fixed 16-step feedback replay also agrees exactly.
There is no supplied reaction–thermal interaction law: isolated reaction and
thermal records do not authorize inventing one. That coupling remains an open
input, documented in the [obstruction note](../../research/repair/evidence_obstruction.md).

## The policy whose status changes

The fixed operating search examines two controls from
\(\{0,1/4,1/2,3/4,1\}\), sorted by total input then lexicographically.
After preparation \((1,1/4)\), its first candidate, \((0,0)\), succeeds.
This is a small feasible-policy construction, not a difficult search or a global
preparation optimization result.

| Case | Final hot | Final cold | Original constraints |
|---|---:|---:|---|
| Repaired prediction for \((1,1/4,0,0)\) | \(19/128\) | \(17/128\) | Pass |
| Independently evaluated original equations | \(19/128\) | \(17/128\) | Pass |
| Initially compatible alternative \(A'\) | \(9/64\) | \(11/128\) | Cold constraint fails |
| Same hold after other preparation \((0,3/4)\) | \(15/128\) | \(3/32\) | Cold constraint fails |

The passing route has cold margin \(21/640\), hot margin \(97/640\), and
total preparation-plus-operation input \(5/4\). It is absent from the original
and acquisition trajectories. Before acquisition, compatible models disagree
on feasibility; afterward, the finite-order certificate establishes the response
for every compatible system in the declared class.

The unchanged bridge consumes the **learned prediction**, with learned response
Lipschitz bounds, and retains the feasible point. Its overall domain-optimization
status remains `unresolved`: the submitted one-point cover information does not
certify a global comparison. No stronger optimality claim is attached.

## Resources, checks, and retained assumptions

The nine original records supply 18 scalar outputs. The additional experiment
would require one reset, eight applied inputs and 16 scalar observations if
performed as a streaming physical experiment. This run instead made eight exact
prefix-evaluator calls, totaling 36 simulated transitions. No laboratory
measurements or hardware operations occurred.

The combined dataset has 17 records, 34 scalar outputs, design rank eight, and
an \(8\times4\) Hankel matrix. Exact construction/checking is finite rational
linear algebra. General costs depend on record count, order bound, port count
and rational bit length. Each deployed state/output update takes
\(O(r^2+pr)\) arithmetic operations and \(O(r)\) mutable state; exact integer
bit lengths can grow with history length. The replay helper caps stored histories
at 256 inputs; the input parser caps \(D\le8\), ports at four, records at 128,
and input rational numerators/denominators at 256 bits. It is a bounded research
implementation, not a universal noise-tolerant inference service.

The recorded wall times were approximately 0.040 s for initial inference,
0.00027 s for added equation calls, 0.046 s for acquisition/construction in the
learner subprocess, and 0.168 s for checking, policy work and fixed evaluations.
These are stage timings from one replay, excluding untimed setup/serialization
and the separate full test suite; they establish no performance advantage.
The 1,093 fixed input words and 16 feedback steps are additional mathematical
regression work, charged to checking/evaluation rather than physical acquisition.

| Source of uncertainty | Status in this result |
|---|---|
| Physical-model applicability | Unbounded; exact LTI/order assumptions retained |
| Representation error | Zero under the checked model/class and resets |
| Numerical error | Zero for accepted exact rational calculations |
| State-estimation error | Zero only with confirmed reset and complete input history |
| Measurement uncertainty | Exact equation data only; uncertain intervals refused |
| New physical coupling | Not identified or certified |

Independent replay checks the evidence binding, incompatible-policy witness,
learned response and bridge receipt. Eight focused tests include rejection of a
self-consistent fabricated model that fails the raw evidence, an independent
three-state example, uncertainty refusal, and malformed records. The full suite
passes 209 tests. Historical source/data/result and certificate-archive hashes
are checked without rewriting the older receipts. Test count is implementation
evidence, not scientific capability.

The precise discharged assumption is: **given the admitted exact finite-order
input/output class and accessible reset experiments, the state coordinates,
initialization procedure and update/prediction map no longer need to be supplied.**
The procedure constructs and checks them. Model completeness, physical evidence,
nonlinear extension, unknown coupling and global alternatives remain unresolved.

## Established foundations and artifacts

The realization is a conventional [Ho–Kalman construction](https://doi.org/10.1524/auto.1966.14.112.545).
The known-equation invariant-row-space check also yields dimension two; this is
the linear specialization of the constrained-lumping principle underlying
[CLUE](https://arxiv.org/abs/2004.11961). No superiority over CLUE, ARX, or other
realization methods was tested or claimed. The distinction between identification
and decision-specific evidence follows the established
[data-informativity framework](https://arxiv.org/abs/1908.00468).
The [primary-source map](../../research/repair/primary_sources.md) separates these
foundations from the still-open research questions.

- [Frozen protocol](../../research/repair/PROTOCOL.md)
- [Observable-data learner](../../experiments/repair_realization.py)
- [Acquisition and policy integration](../../experiments/repair_run.py)
- [Result and resource ledger](result.json)
- [Independent replay receipt](verification_receipt.json)
- [Implementation audit](../../research/repair/implementation_audit.md)

From the project root, reproduce with the configured Python runtime:

```text
python3 -m experiments.repair_run
python3 -m unittest discover -s tests -v > results/repair/tests.log 2>&1
python3 -m experiments.repair_verify
```
