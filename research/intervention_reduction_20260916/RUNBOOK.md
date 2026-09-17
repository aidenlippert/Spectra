# Reproduce the intervention-reduction pass

Run from `/Users/aidenlippert/Documents/Spectra`. The exact accepting modules
use the Python standard library and the repository CAR algebra. They do not
import NumPy, SciPy, PySCF, quimb, or numba. Numerical proposals and references
use the existing `.venv-correlated` environment. No installation, paid compute,
Git push, or external service is needed.

The fixture and MPS inputs are under
`/Users/aidenlippert/Documents/Spectra/results/transfer_solver_20260915/cases`.
The exact-model and inherited-state digests are in ACCOUNTING.json. Preserve
these inputs. The commands below create new outputs and refuse to overwrite
existing ones; select fresh names when replaying again. The bounded runner
sets numerical thread counts to one and records every instrumented process.

## Focused tests

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -B -S -m unittest \
  research.intervention_reduction_20260916.test_exact \
  research.intervention_reduction_20260916.test_trajectory \
  research.intervention_reduction_20260916.test_robust \
  research.intervention_reduction_20260916.test_kernel -v

.venv-correlated/bin/python -B -m unittest \
  research.intervention_reduction_20260916.test_predict -v
```

Run numerical proposer tests separately; the accepting checker deliberately
refuses a process that has imported numerical packages.

## Fresh exact reconstruction and three reusable H6 queries

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -B -S \
  -m research.intervention_reduction_20260916.budget \
  --name reproduce_h6_bundle --seconds 180 -- \
  /opt/homebrew/Caskroom/miniconda/base/bin/python -B -S \
  -m research.intervention_reduction_20260916.bundle \
  results/transfer_solver_20260915/cases/h6_asymmetric/fixture.json \
  results/intervention_reduction_20260916/h6_snapshot/proposal_r16.json \
  results/intervention_reduction_20260916/reproduced_h6_bundle \
  results/intervention_reduction_20260916/h6_snapshot/switch_r16.json \
  results/intervention_reduction_20260916/h6_snapshot/interior_r16.json \
  results/intervention_reduction_20260916/h6_snapshot/three_r16.json
```

This reconstructs the action kernel from the original Hamiltonian and integer
columns, then drops the configuration data before query acceptance. It does not
trust a previously printed receipt or an imported Gram matrix.

## Replay the designed pulse and its 5% functional robustness

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -B -S \
  -m research.intervention_reduction_20260916.budget \
  --name reproduce_h6_policy --seconds 180 -- \
  /opt/homebrew/Caskroom/miniconda/base/bin/python -B -S \
  -m research.intervention_reduction_20260916.bundle \
  results/transfer_solver_20260915/cases/h6_asymmetric/fixture.json \
  results/intervention_reduction_20260916/h6_snapshot/proposal_r24.json \
  results/intervention_reduction_20260916/reproduced_h6_policy \
  results/intervention_reduction_20260916/h6_design/candidate_1_r24.json \
  --max-deviation-Ha 1/2000
```

The expected functional lower endpoint is greater than 1/40. The printed
`target_met` flag is the separate state-error target and is false at this larger
perturbation allowance. CONTROL_POLICY.json records these two criteria explicitly.

## Discovery and refinement

`discover.py CASE OUTPUT --budget 512 --ranks 4 8 16 --selection snapshot`
rebuilds the original snapshot proposals from the inherited MPS and fixture.
It charges its selected-configuration eigensystems. It is not a cold molecular
calculation from geometry and does not avoid determinant labels.

`enrich.py FIXTURE PARENT CHILD --added 8` retains all parent columns exactly,
adds directions from driven residuals on the training corners, and charges all
reached labels. `predict.py FIXTURE PROPOSAL OUTPUT --policy heldout_switch`
creates a candidate trajectory. `--stable` selects stable Taylor propagation;
`--chebyshev` selects the interpolation proposer used in the final water batch.
Use `--schedule-from TRAJECTORY` to preserve a previously selected pulse while
changing the reduced model.

`design.py FIXTURE PROPOSAL OUTPUT_DIRECTORY` runs the finite 729-policy search.
Its three exported candidates are numerical proposals; each still needs exact
acceptance. `trajectory.py` performs direct physical-action replay; `robust.py`
adds the control-neighborhood guarantee; `bundle.py` compiles once for many
queries. `compare.py` compares all mathematical fields of direct and compiled
receipts. Every entry point supports `--help`.

## Independent reference and audit

`oracle.py FIXTURE PROPOSAL TRAJECTORY RECEIPT OUTPUT` builds a full magnetic-
sector matrix with an independent ladder routine and compares sparse numerical
propagation to the certificate. It explicitly enumerates its reference sector.
It is a diagnostic, not the proof checker and not an undisclosed teacher.

`python -B -S -m research.intervention_reduction_20260916.audit` checks bindings,
inherited-column preservation, final claims, and every available process ledger.
It writes ACCOUNTING.json and CONTROL_POLICY.json. It does not replace exact
replay. Re-running it after new experiments changes the ledger totals; the
report's quoted totals refer to the completed original pass.

The accepting boundary imports `experiments/marginal_symbolic.py` for CAR
polynomial parsing and Hermiticity. The process runner reuses
`research/correlated_pair_20260913/budget.py`. A final SHA-256 manifest records
the new source, artifacts, and these dependencies. The manifest detects later
changes; it is not a mathematical certificate or a novelty claim.
