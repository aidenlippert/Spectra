# Reproducing the interacting-certificate campaign

All commands run from the Spectra repository root. The measured host uses the
existing `.venv-correlated` numerical environment, `.venv-molecule` chemistry
environment, and isolated `.venv-interacting-libs` overlay. The numerical runner
pins one BLAS/Numba thread and one heavy process at a time. Its interpreter paths
describe this local host; they are not a portable environment installer.

## Replay the accepted original-H8 certificate

The accepting path needs Python's standard library and the repository's checker
modules. It does not need NumPy, SciPy, PySCF, quimb or the optimization overlay.
Use a new output directory, outside the sealed campaign:

```bash
python -B -S -m research.interacting_scaling_20260915.replay_bundle \
  results/interacting_scaling_20260915/models/h8_matched_direct \
  results/interacting_scaling_20260915/cases/h8_matched_direct_level0 \
  /tmp/spectra-h8-interacting-replay
```

This copies the proposed witnesses, rederives the exact orbital rotation,
rechecks the MPS upper, both spin-domain lower proofs and the complete remainder,
then compares both exact original-H endpoints with the retained claim. It writes
no files into the sealed source case. A previous interval or success flag alone
cannot satisfy it. The output directory must not already exist.

For a retained canonical reference, use the existing `strict_replay` caller with
its `--rotated` argument, naming that run's `_base` local state directory. The
reference's copied local MPS is a moment guide, not a canonical MPS witness.

## Fresh discovery

Before sealing, a new unique development run can be started with:

```bash
python -B -S -m research.interacting_scaling_20260915.cold unique_h8_run \
  results/transfer_solver_20260915/cases/h8_cold/specification.json \
  --solve-seconds 450
```

Add `--reference` and use `--solve-seconds 840` for the extended retained canonical
reference, with a newly discovered local state and transported moments. Use
another unique name. Both
paths regenerate integrals; they must be compared on equal original inputs.
Every stage and failed adaptation level has a bounded receipt. The full parent
clock includes all stages and intervening orchestration; stage sums are separate.

Here, fresh/cold refers to new molecular integrals, state discovery and proof
discovery, not a reboot or reinstall of the numerical software. OS caches and
reusable compiled library kernels may already exist. Initial environment work
is charged to the campaign ledger rather than repeated in every cold clock.

The preserved MPS routine has a generic `initializer` description mentioning a
random MPS even when an explicit seed is supplied. Its actual `initial_MPS_path`,
the calling command and the retained seed are authoritative for these runs.
The cold wrapper supplies a new globally charge-correct product seed; coupling
continuations are explicitly recorded as such. The historical routine and its
original output strings were not silently rewritten after the protocol froze.

After sealing, create a new campaign namespace before doing new discovery.
The existing runner refuses to append to a sealed campaign. Do not delete or
rename a manifest to bypass that rule. Earlier state/proof inputs are allowed
only in explicitly labeled conditional or continuation tests.

## Frozen held-out tests and coupling

`results/interacting_scaling_20260915/heldout_protocol.json` binds the constructor,
models, nested levels and resource caps before held-out results. The supported
test names are `h12_heldout`, `water_heldout`, and `h4_631g_heldout`. The latter
changes the orbital-space model, not just molecular geometry.

The separately timestamped `capacity_protocol.json` was also declared before
H12 inputs were generated. If the initial H12 attempt refuses only the matrix
or row envelope, it permits a fresh retry with up to 4,000,000 Gram entries and
400,000 coefficient rows. It changes no operator supports, ideal equations,
state algorithm or solve limits. Each pipeline has a 3600-second cap, and both
attempts are charged together. Peak RAM is measured, not hard capped; the
permission was based on the preceding H10 run's approximately 1.01 GB peak.

The adaptive coupling driver uses the same frozen nested proof rule at lambda
0, 1/4, 1/2 and 1. It permits all global-N charge distributions. Only the state
is continued between couplings; optimization starts independently at each one
and can warm-start between exactly nested levels of that same Hamiltonian.
Only lambda=1 is transferred to the original molecular input. Intermediate
models are diagnostic support scalings.

`fragment_control.py` is a separate zero-coupling comparison. It enumerates each
small fragment over every charge and spin projection, then combines local bounds
over the total particle number. It refuses nonzero coupling. Its local losses
are retained when summed. It is not an input to the interacting state or proof.

## Audit artifacts

* `runs/`: every bounded child attempt, including failures and installation costs.
* `executions/`: declared stage order, stage limits and completion status.
* `cold/`: fresh parent clocks, input protocols and every attempted level.
* `cases/*/prepared/`: directly generated maps and representation counts.
* `cases/*/exact/`: independently checked local or canonical certificate bundles.
* `cases/*/original_interval.json`: exact two-sided original-model transfer.
* `duals/`: exported family, affine/nullspace repair and independent exact check.
* `coupling/`: coupling definitions, continuation dependencies and all outcomes.
* `physical_model_ladder.json`: exact solver subtraction and separate numerical
  chemistry controls; no certified physical-error interval is claimed.
* `posthoc_h12_result.json` and `residual_h12_result.json`: separately declared
  optimizer continuations, with unchanged source maps and explicit earlier costs.
  These are not substituted for the frozen H12 outcome.
* `h12_recovery_result.json`: the main H12 attempt after correcting the caller's
  certificate handoff. Its driver uses the existing `--nonsinglet-source` argument;
  the earlier failed initialization and the frozen constructor remain preserved.
* `h12_main_refinement_result.json`: one additionally declared main-proof residual
  phase, gated on the preceding complete replay and charged beyond the initial
  diagnostic envelope. It uses identical stored maps and state, with fresh exact
  replay of any resulting candidate.

Final file integrity is checked by `seal.py check`. Hash integrity is not itself
a proof of the mathematical claims; the exact replay tools provide that check.
