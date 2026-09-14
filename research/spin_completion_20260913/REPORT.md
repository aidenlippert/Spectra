# Spin-changing interference improves the compact H6 bound

**Sixty-four coupled combinations narrow H6's certified interval from
13.010569 to 7.414219 mHa**, a 43.01% reduction. The same combinations used
separately give 9.486940 mHa: coupling improves the found bound by another
**2.072721 mHa**. The new lower exceeds the entire previous diagonal-spin
family's exact ceiling by **2.566609 mHa**.

This is a useful missing correlation mechanism, but the **1.6 mHa target
remains unmet**. Both complete-family controls timed out without a lower
certificate. This pass therefore establishes neither the new family's
optimum nor an obstruction to reaching the target within it.

The original question remains the organizing test: can the difficult part
of the electron problem be represented by a small collection of interference
patterns, with the remainder certified collectively? Here the ten spatial
Hamiltonian patterns and their **0.108331264 mHa** collective remainder are
unchanged. The new work improves the joint correlation constraints attached
to those patterns. It does not yet supply a small sufficient description.

## What was added and why

The previous family admitted Q_k^(up,up) a_i and Q_k^(down,down) a_i.
The new one also admits Q_k^(up,down) a_i and Q_k^(down,up) a_i, using
exactly the same spatial matrices. These operators probe correlations that
exchange spin between orbital components. Their positivity conditions apply
to every physical state; they are not instructions to manipulate a state.

The complete list contains 492 generators, partitioned into four blocks of
dimension 93 and four of dimension 30. Every previous 252-frame generator is
included. No numerical basis truncation was used.

Replaying the previous full-frame dual and evaluating the new positivity
conditions produced **sixteen exactly verified negative directions**. The
strongest normalized expectation is about -0.01784948. These directions
exclude the prior nonphysical moment assignment. Their value as energy
constraints was then tested by separate optimization and exact lower replay.

## Small coupled sets

The fresh adaptive search started from the same full quadratic baseline.
Its first round used the exact counterexample's separators; later rounds
used the current numerical dual as a proposal. At each of four rounds, both
one and two new independent directions per symmetry class were optimized
and exactly replayed. Continuation used accepted energy gain per measured
trial cost plus pricing. Every round selected the larger, sixteen-direction
bundle, so this does **not** demonstrate an advantage over a fixed sixteen-
direction schedule.

| Selected correlation combinations | Certified interval, mHa |
|---|---:|
| 0: fresh quadratic baseline | 13.010569 |
| 16 | 12.280137 |
| 32 | 10.465265 |
| 48 | 8.850775 |
| **64, coupled** | **7.414219** |
| Same 64, separate contributions | 9.486940 |

Four additional trial certificates, at 8, 24, 40 and 56 combinations, are
retained and charged. The search stopped at its four-round limit after
**76.177 seconds**, before exhausting its 240-second budget. A forty-
combination result already exceeded the previous full-family ceiling.

At the same count of 32 combinations, the previous diagonal-spin search
gave 11.429079 mHa, versus 10.465265 here. This is an observed comparison of
different discovered directions and different discovery ancestry, not a
matched policy benchmark. The exact ceiling comparison is the stronger
evidence that the added family contributes new information.

The 64-combination molecular lower is -6.340472845570542... Ha and the
separately replayed rational upper is -6.333058626233001... Ha. Although the
best numerical proposal was `optimal_inaccurate`, its exact coefficient-L1
penalty is included: 0.001460231 mHa. Numerical status does not certify the
result; the reconstructed rational SOS and tail proofs do.

## Representation and complete costs

| Quantity | Selected 64 | Full 492 |
|---|---:|---:|
| Added joint blocks | 8 × 8 | 4 × 93 and 4 × 30 |
| Gram entries, including quadratic baseline | 6,236 | 43,920 |
| Numerical map nonzeros, conditioned | 251,268 | 12,129,902 |
| Numerical map nonzeros, full control without conditioning | — | 2,509,230 |
| Compact accepted lower certificate | 94,565 B | none exported |
| Expanded accepted SOS certificate | 645,205 B | none exported |

The selected directions contain **3,936 nonzero frame coefficients**, and
the compact certificate contains another **4,736 nonzero factor coefficients**
including the baseline. These counts are distinct from the ten Hamiltonian
patterns. The sizes above exclude the frozen Hamiltonian, its 6,683-byte
retained model, and its 3,667-byte tail certificate.

The last selected trial took 0.517 seconds to construct, 1.726 seconds to
solve, 0.031 seconds to export, and 5.921 seconds for exact acceptance.
The 76.177-second adaptive total also includes imports, shared construction,
the baseline, all candidate trials, pricing, exports and writes.

The exact separator diagnostic took **45.850 seconds**. Including it makes
the current small-search cost **122.027 seconds**. Its source counterexample
also required a previous full-frame solve and dual construction totaling
100.484 seconds; including that ancestry gives **222.511 seconds**. The
separate-contribution control took another **13.309 seconds** and inherits
the selected directions' entire discovery cost.

The original full-family control hit its outer limit at **240.278 seconds**
without an exported certificate. Its whitening transformation produced a
dense coefficient map. A separately declared full control omitted that
transformation while retaining every operator and PSD entry. Map nonzeros
fell by **79.31%**, but this control also timed out, at **240.339 seconds**.
Both attempts remain failed results. Map sparsity alone did not deliver a
completed full solve under the chosen budget.

Total new diagnostic, adaptive, ablation and full-control work was
**615.970 seconds**, including both failed controls. This excludes fresh
validation, development work and the 100.484-second prior counterexample
ancestry. The old rational upper and frozen molecular representation have
their own previous discovery costs. These are single-run observations;
BLAS thread environment variables were set to one, while Clarabel retained
its automatic thread setting. No general speedup or scaling claim follows.

## What this means for the next pass

There is now certified evidence that spin-changing **joint** patterns remove
part of the missing correlation information. The old family's obstruction
does not apply to the enlarged family: an accepted proof already exceeds it.

The next useful experiment is to continue enriching the **small coupled
subspace**, using its current dual's remaining violations and a fixed total
budget. Its final optimization took about 1.7 seconds, whereas neither dense
full-family control completed. The present four-round result leaves the
accuracy question open. More directions are justified only by actual
certified gain and total cost. Switching families again should be motivated
by new evidence; no exact obstruction for this completed spin family has
been established.

This remains one frozen finite-basis H6 test. It establishes no general
many-body compression theorem, experimental accuracy, or physical synthesis
capability. No many-body sector was enumerated in the new discovery or
acceptance. The upper witness's prior determinant-amplitude discovery is
outside the compact lower method.

## Evidence and reproduction

The [proof](/Users/aidenlippert/Documents/Spectra/research/spin_completion_20260913/PROOF.md)
describes the exact operator and residual acceptance. The original
[protocol](/Users/aidenlippert/Documents/Spectra/research/spin_completion_20260913/PROTOCOL.md)
and the separately declared
[full-control follow-up](/Users/aidenlippert/Documents/Spectra/research/spin_completion_20260913/RAW_CONTROL_PROTOCOL.md)
record the limits and preserve both failures.

Nine focused tests cover diagonal-spin containment, correct spin-changing
CAR and adjoint order, exact factor composition, sixth-degree cancellation,
tail transfer, invalid-input refusal, exact versus numerical pricing maps,
complete-frame inclusion, accepted/rejected separator conditions, and refusal
to overwrite existing full-control results or logs. The seven exact tests
run with `python -S`; the two map tests use the numerical
environment. The previous pass's **132 manifest files remain unchanged**.

All **ten molecular intervals and sixteen separators** passed fresh exact
replay in **92.605 seconds**, together with the complete previous dual ceiling
and rational reference upper. No NumPy, SciPy, CVXPY or PySCF was loaded on
that accepting path. This is recorded in
[fresh_replay.json](/Users/aidenlippert/Documents/Spectra/results/spin_completion_20260913/fresh_replay.json).
All costs, certificates, proposals, timeout metadata, source snapshots and
hashes are in
[results/spin_completion_20260913](/Users/aidenlippert/Documents/Spectra/results/spin_completion_20260913).

From the workspace root with the configured Python environment:

```sh
python -S -m unittest research.spin_completion_20260913.test_core \
  research.spin_completion_20260913.test_diagnostic -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest \
  research.spin_completion_20260913.test_discovery -v
python -S -m research.spin_completion_20260913.replay \
  --campaign results/spin_completion_20260913/campaign \
  --out /tmp/spin-completion-replay.json
```

`diagnostic.py` repeats the separator computation; its recorded run used an
external 180-second watchdog. `campaign.py` and `raw_control.py` launch their
bounded stages into new output directories. Lower discovery refuses existing
output directories. All accepting paths use the standard library and the
existing CAR, collective-tail, dual and reference-upper verifiers.
