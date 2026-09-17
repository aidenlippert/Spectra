# Reproduce the acceptance and selected-channel continuation

Run commands from `/Users/aidenlippert/Documents/Spectra`. Every result is
conditional on the shared H12 inputs recorded in `results/acceptance_channels_20260915/protocol.json`.
The original campaigns are preserved. This directory does not initiate any GitHub operation.

The numerical interpreter is `.venv-correlated/bin/python`; the installed
SciPy version is 1.18.1 and SCS is 3.2.11. SuiteSparseQR is loaded from the
already installed `.venv-interacting-libs` overlay. The numeric runtime cache,
attempted separately, timed out while copying; its partial copy was never added
to an import path. All numerical runs used the existing installed packages.
Pin BLAS/OMP threads to one. The campaign runner does this and records each
subprocess, timeout, wall/CPU time, peak child RSS and failure.

Focused algebra, conic-formulation, and refusal-path tests:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv-correlated/bin/python -B -m unittest discover -s research/acceptance_channels_20260915 -p 'test_*.py'
```

Replay an exported four-direction certificate in a new folder, checking its
actual upper and lower witnesses and transfer to the original Hamiltonian:

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -B -S -m research.acceptance_channels_20260915.replay_candidate dense_t2 --name dense_t2_fresh
```

The checker requires the exact original model in
`results/interacting_scaling_20260915/models/h12_heldout`, the fixed inputs under
`results/acceptance_channels_20260915/h12`, and the actual proposed certificate.
No numerical optimization output or supplied success flag is accepted as proof.
Use a new `--name` for every invocation; existing receipts are never overwritten.
`matched_continuation` selects the unchanged-family comparison certificate.

Each executed construction command is saved verbatim in
`results/acceptance_channels_20260915/runs/*.json`; the companion `.log` contains
its output. These are the authoritative run commands and accounting, including
unsuccessful experiments and storage-related retries. Repeating a numerical
command requires a fresh proposal tag and a fresh run name. A saved checkpoint
is the last iterate; the exported `raw.npz` is the best accepting-score iterate.
Do not silently interchange them in a matched comparison.

The main entry points are:

| Module | Role |
|---|---|
| `accepted_fit` | Restricted L1 fits with fixed Gram matrices or nonnegative whole-block rescaling |
| `scs_control` | Separate numerical algorithm for the unchanged equality SDP |
| `scs_accepted` | Full accepted-L1 conic formulation with every original Gram variable free |
| `channel`, `channel_membership` | Exact supplied-channel algebra and constructive membership in the current cones |
| `select_channels`, `dense_t2` | Selected spin-change-three directions, quartic maps and nested small Gram extensions |
| `mixed_t2` | Separate spin-change-one diagnostic and joint solve preserving the earlier block |
| `map_consistency` | Direct polynomial evaluation compared with the molecular coefficient map |
| `trace_diagnostic` | Numerical trace/nullspace diagnostic; never an accepted obstruction |
| `replay_candidate` | Fresh complete original-model replay with the original standard-library checker |
| `paired_exact`, `paired_replay` | Experimental exact paired-contraction path, requiring comparison with the original full molecular receipt |
| `iteration_comparison` | Read-only comparison of common logged iterations, separating equal-time outcomes from a claim about the added correlations |
| `residual_diagnosis` | Exact arithmetic on the accepted receipt; limits only residual-scalar replacement for that fixed identity |
| `singlet_null_channels` | Exact coordinates for 24 physical singlet-null directions in two actual H12 dictionaries; no dual repair or energy claim |

The exact external H8 source dual described in the attachment is absent. The
supplied polynomial and its current-cone membership can be checked locally;
the external negative dual evaluation cannot be reproduced without that bundle.
Current-cone membership uses the supplied polynomial in current orbital labels;
it does not establish a basis match or transport of the external physical operator.

The local input cache under `~/.cache/spectra/acceptance_channels_20260915`
avoids offloaded Documents files during numerical work. `input_cache.json`
records each copied input and its hash, including comparison with the sealed
prepared-file hashes. Accepted replay folders copy their needed inputs and do
not rely on an ephemeral runtime cache. Cache population and failed reads are
included in the new metered stages, not presented as an algorithmic speedup.

The attempted optional SciPy cache did not complete within 300 seconds. Its
timeout is retained in the ledger; it is not a working runtime or a speedup result.

The original dense-candidate replay first timed out after 900 seconds. The
separate `extended_replay` caller grants 1,800 seconds to the unchanged exact
checker and records that explicitly; it changes no energy-acceptance condition.
The original failed run remains part of the cost. A paired-contraction replay
uses a separate output name and must reproduce the original exact lower and
remainder before a verifier speedup is reported.

The final report combines exact replay receipts, timed run records and the two
read-only diagnostics. `summarize` refuses to finalize while a metered run is
still active. Its cost total is a continuation ledger including failed attempts,
not a cold calculation from integrals. Do not attribute the equal-time interval
difference to the selected directions without reading `iteration_comparison.json`:
the four-direction run completed more iterations and was worse at common ones.

SCS API and packing references: [semidefinite cones](https://www.cvxgrp.org/scs/api/cones.html),
[Python interface](https://www.cvxgrp.org/scs/api/python.html), and
[solver settings](https://www.cvxgrp.org/scs/api/settings.html).
