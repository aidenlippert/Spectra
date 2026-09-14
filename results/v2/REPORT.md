# Complementary-learning result receipt

The finite model theorem passes its declared exact and acquisition gates. The broader scientific-intelligence claim remains open. Read [the proof](../../research/v2/PROOF.md) for assumptions and derivations.

## Exact result

| Retained response laws | Optimal one-shot error |
|---|---:|
| none | 7/16 = 0.4375 |
| A | 3/8 = 0.3750 |
| B | 3/8 = 0.3750 |
| AB | 1/4 = 0.2500 |

The complementary saving is **1/16**, or 6.25 percentage points. The second law reduces error by 1/16 before the first law and 1/8 after it. This is a statement about acquired information, not superiority of the inference algorithm.

Exact policy verification checks 520 finite Bellman nodes. The chain verifier checks 27 combinations of factor count and retained laws through six factors, over each complete product-action menu. Its adaptive lower-bound cross-check includes horizons zero through five.

## Acquisition and transfer

Each world used **2,784 acquisition shots**: 174 direct A shots and 2,610 B reference shots compiled using A. The exact failure union bound is 7.458653775e-05 per world. The unconditional AB error upper bound is 0.250018647; the conservative complementary-saving lower bound is 0.062481353. These are coverage guarantees conditional on the declared source, model class and shot independence.

| Hidden-world seed | Both masks recovered | Unused control pairs checked exactly | Empirical frozen error | Empirical AB error |
|---|---|---:|---:|---:|
| 17 | True | 3249 | 0.448730 | 0.254272 |
| 29 | True | 3249 | 0.438232 | 0.253906 |
| 43 | True | 3249 | 0.437744 | 0.248901 |

Each row has 8,192 sampled target trials. Conditional expected errors match the exact theorem in every world. The JSON contains simultaneous Hoeffding intervals across all 21 method/world results. Sampled outcomes corroborate the analytic result; they do not establish its universal mathematical part.

Exact-answer lookup has zero held-out hits and matches the frozen baseline. Structured retrieval and independent full-feature GF(2) inference match AB actions and errors exactly. The A-only and B-only results are information-erasure ablations after acquisition; they are not competitors denied shared data without disclosure.

## Resource comparison

The optimal reset target-only four-shot error is **103/256**, above the chosen threshold 3/10. Any fixed-horizon reset target-only strategy meeting the threshold therefore needs at least five shots. The acquisition-plus-one-shot strategy meets the threshold with budget `2784+N`, crossing the conservative `5N` bound at **697 tasks**. For 8,192 tasks, budgets are **10,976 versus at least 40,960**.

This comparison is limited to fixed-horizon target-shot budgets. It is neither an expected-stopping-time bound nor an advantage over a globally learning conventional agent. Primitive preparations have different qubit counts. The JSON separates qubit preparations, retained records, random draws, arithmetic operations and measured runtime; no apparatus energy or hardware time result is claimed.

## Representation revision: a separate result

The sampled-only learner detects interaction `(0,1)`, then `(1,2)`, and predicts all unused higher-weight controls in the two tested quadratic systems. A cubic term agrees on every fitting context but triggers a held-out mismatch and abstention. A same-data full-feature conventional solve agrees with both fitted models.

Warm and reset acquisition both cost **682 shots** across the two systems. Held-out diagnostics add 270 shots; the cubic negative control adds 476. The measured warm saving is **zero**. Interaction discovery is demonstrated within the finite admitted vocabulary; sample compounding from representation revision is not.

## Verification and artifacts

57 focused and existing regression tests pass. They include independent dense quantum checks, exact likelihood and tail arithmetic, corrupted-certificate rejection, rank and sample refusals, wrong-prefix channel failure, unused-context transfer and nonlinear misspecification.

- [Complete observations and policy certificates](compounding_results.json)
- [Compact result ledger](summary.json)
- [Chain arithmetic checks](chain_certificates.json)
- [Revision observations and mismatch evidence](revision_results.json)
- [Protocol fixed before the integrated run](../../research/v2/protocol.json)
- [Statistical audit](../../research/v2/statistical_audit.md)
- [Robustness and composition inequalities](../../research/v2/compositional_transfer_bound.md)

No physical experiment, new physical primitive, universal solver or general scientific intelligence has been established. The next open gate is reusable strategy or representation discovery on broader coupled physics with an advantage that survives matched stateful baselines.
