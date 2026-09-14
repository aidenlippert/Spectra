# Reproduce the local-operator pass

Run from `/Users/aidenlippert/Documents/Spectra`. The recorded environment is in
`results/molecular_identity_20260913/environment.json`. The numerical selector
uses the existing NumPy/SciPy/CVXPY environment; accepting replay uses only the
standard library and exact CAR routines.

## Recheck the exact structural obstruction

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.molecular_identity_20260913.replay \
  --scan results/molecular_identity_20260913/initial_scan.json \
  --witness results/certificate_scaling/commutator_dictionary/h4_exact_dual/witness.json \
  --out results/molecular_identity_20260913/obstruction_rechecked.json
```

## Recheck an energy interval

The same command applies to each saved certificate; choose the matching H4/H6
reference. Replay reconstructs every rational square, the global-number ideal,
the complete residual, and the independent rational upper expectation.

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.certificate_scaling.cubic_interval_replay \
  --certificate results/molecular_identity_20260913/h4_coupled/round_3/certificate.json \
  --reference results/certificate_scaling/active_space_ladder_references_aligned/h4/upper.json \
  --out results/molecular_identity_20260913/h4_interval_rechecked.json
```

## Repeat discovery in a new output directory

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/homebrew/Caskroom/miniconda/base/bin/python \
  -m research.molecular_identity_20260913.coupled \
  --fixture results/certificate_scaling/active_space_ladder/h4/fixture.json \
  --out results/molecular_identity_reproduction/h4
```

Replace both `h4` components with `h6` for transfer. The module writes its fixed
policy, starts from the creator-channel cone, prices all four-orbital subsets,
adds up to eight negative cubic directions per round and their adjoints, and
performs three enrichment rounds. Each solve requests a 60-second solver limit;
construction, canonicalization, export, and pricing are additional recorded costs.
Floating solver output and near-degenerate eigenvector choices may differ across
versions or machines; every resulting certificate still needs exact replay.

The separate-local-block control uses
`research.molecular_identity_20260913.adaptive` with the same fixture/output CLI.
It starts with all three-orbital blocks and adds eight four-orbital supports in
each of two rounds. `local_blocks` produces its initial all-three-orbital control.
The saved H4 separate-block run reused its already measured initial control;
its seed cost remains in the ledger.

For the H4 coupling ablation, decode
`h4_coupled/round_3/additional_groups.json`, remove each `merge_into` field, and
pass those exact polynomial groups to the existing commutator solver with the
same creator-channel base. This keeps the learned operators and removes cross
terms with the base; its dictionary-discovery cost is the preceding coupled run.

## Focused controls

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.molecular_identity_20260913.test_local
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/homebrew/Caskroom/miniconda/base/bin/python \
  -m unittest research.molecular_identity_20260913.test_discovery \
  research.certificate_scaling.test_commutator_dictionary
```

`audit` collects completed experiments, replays missing interval receipts, and
checks hashes of certificates/references that were already replayed. It refuses
unfinished experiments or changed frozen inputs. The final manifest records
the source snapshot and artifact hashes; the previous commutator solver source
and the narrow extension diff are preserved alongside results.

H6 initially exposed an unhandled symmetry channel before round-two construction.
The exception log and failure receipt are retained. The fix opens a new PSD block
when no base block has the requested charge/symmetry, and the saved run resumes
without repeating selection or completed solves. H4 selected no such channels.
This changed implementation coverage, not the numerical selection policy. The
recovered pre-failure wall time is explicitly labeled an estimate from file
timestamps; the time spent diagnosing the implementation bug is not solver time.
