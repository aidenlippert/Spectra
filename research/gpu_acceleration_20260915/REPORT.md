# Measured acceleration: 2.8× faster complete H6 calculation

The selected fresh H6 run completed in **146.792 seconds (2.45 minutes)**,
including new molecular integrals, state discovery, both lower-proof
constructions, optimization, and exact verification. The original local
reference took **411.123 seconds (6.85 minutes)**. The measured combined
improvement is **2.80×**; it is not a GPU-only speedup.

The new exact interval is **0.054185855 mHa** for the newly generated rational
Hamiltonian. An exact orbital-sign comparison and coefficient-norm allowance
transfer it to the original asymmetric H6 fixture at **0.054185895 mHa**.
Both are below the 1.6 mHa requirement. The original, tighter approximately
0.021517 mHa result remains intact. This is a comparison at the requested
accuracy threshold, not at identical interval widths.

## What improved

1. **Earlier numerical refinement.** Switching mu from 2 to .03 once the
   predicted interval is below 10 mHa and the primal residual below 1e-4 avoids
   spending the entire original 300-second first-stage budget. The existing
   numerical stopping tolerances and final exact acceptance remain.
2. **Selective CUDA work.** Blocks smaller than 128 stay on CPU; larger blocks
   use FP64 CUDA eigensolves. Blindly moving every block to the GPU was worse.
3. **Existing faster exact projection.** The project already had a canonical
   CAR-word projection routine. Using it during preparation, with Python 3.12,
   reduced the H6 projection measurement from about 127 seconds to 25 seconds.
   Its entries and selected rows matched the frozen original exactly. This
   combines a runtime change and an existing algebraic implementation; neither
   is credited to GPU hardware. The accepting projector remains unchanged.
4. **Concurrent independent construction.** Nonsinglet and retained-map
   preparation run in separate processes, with disjoint outputs and JIT caches.
   This reduced the complete run from 179.212 to 146.792 seconds. All 38 checked
   fixture/state/proof/array files were byte-identical to the serial run, and
   the exact upper bounds were equal.

## Separate the GPU contribution

These component comparisons used the same A100 host, identical prepared
inputs, the same refinement schedule, and one CPU BLAS thread. They include
optimizer setup, transfers, recording, and export; exact replay is separate.

| Component | CPU, improved schedule | GPU for all block groups | Selected hybrid | Hybrid runtime reduction vs CPU |
|---|---:|---:|---:|---:|
| H6, 2,200 iterations | 39.061 s | 43.591 s | **33.478 s** | **14.3%** |
| H8, 2,800 iterations | 138.549 s | 123.867 s | **111.629 s** | **19.4%** |

The selected hybrid outputs passed the original exact checker at
**0.054172531 mHa on the frozen H6 input** and **0.642880121 mHa on the frozen
H8 input**. CPU/GPU checkpoint differences were below 3.3e-12 relative to the
specified normalization. Both cases informed development; this is not a
held-out reliability study. H8's complete fresh pipeline was not benchmarked
in this pass.

Against the old recorded two-stage optimization times, the selected component
is 9.38× faster on H6 and 2.97× on H8. Those larger ratios combine the changed
schedule, implementation, and hardware. The old runs used the original local
machine; the paired CPU/GPU table above is the hardware comparison.

## Complete fresh H6 accounting

| Step | Elapsed seconds |
|---|---:|
| New integrals | 2.318 |
| New MPS discovery | 37.576 |
| First exact upper check | 1.867 |
| Nonsinglet proof construction | 33.923, overlaps the next row |
| Retained preparation | 43.038, overlaps the preceding row |
| Hybrid optimization and export | 34.570 |
| Complete exact replay, including upper recheck | 27.413 |
| **Total critical-path time** | **146.792** |

The procedure and source hashes were frozen before this run. No earlier MPS,
coefficient maps, optimization checkpoint, or proof supplied its inputs.
Software installation, provisioning, downloads, development failures, and an
additional local audit are separate from this per-solve time. CUDA software
caches were already warm. The final cloud and original local environments
differ, so this full-run ratio is not an isolated hardware comparison.

After downloading, **505 files matched their remote SHA-256 hashes**. A second
local replay with the frozen original standard-library checkers reproduced
both rational endpoints exactly in **16.636 seconds**, importing no numerical
libraries. That additional replay is not included in the 146.792-second cloud
solve. The original-fixture comparison charges a norm allowance of
1/50,000,000,000 hartree to each endpoint.

## Unsuccessful and exploratory routes

- Ordinary batched GPU eigensolves were slower; more CPU threads and independent
  block threads did not establish a useful improvement. One early in-process
  BLAS thread-count probe crashed and was retained as a failed probe.
- CUDA initially lacked toolkit headers; explicit toolkit dependencies fixed it.
- Five fresh-run attempts stopped during environment setup: an older quimb API,
  a virtual-environment path error, missing networkx, missing CVXPY, and an SCS
  option incompatibility. Their logs and costs remain recorded. Dependency and
  solver-option checks passed before the successful complete runs.
- The first complete configuration using older CPU software took about 294
  seconds, despite its faster optimizer. It exposed preparation as the next
  bottleneck; the later complete runs measured the fixes.
- GMP rational probes reproduced exact endpoints but gave variable timing
  improvements. GMP integer transfers did not establish a gain. This remains
  experimental and is not part of the selected accepting path.

The selected implementation is an isolated research component, not a claim of
globally optimal performance, general scaling, or new chemical predictive
accuracy. Further gains will depend on the actual block sizes and remaining
CPU construction costs.

## Reproduction and evidence

Use the commands in [README.md](README.md), with the final cold-run options
`--fast-prepare --parallel-builds`. Final software: Python 3.12.14, NumPy 2.2.6,
SciPy 1.18.1, CuPy 14.2.0, quimb 1.15.0, PySCF 2.14.0, CVXPY 1.9.2, and
SCS 3.2.11. Four projection contract tests passed on CUDA and the CPU fallback.

Evidence is under `results/gpu_acceleration_20260915/`:

- `cases/cold_h6_hybrid_v8/`: frozen protocol, per-stage logs, complete execution,
  exported certificate, exact interval, local replay, and original-fixture bridge.
- `measured_summary.json` and `backend_differential.json`: paired component results.
- `parallel_build_differential.json`: serial/concurrent artifact equality.
- `remote_logs/`: kernel trials, dependency logs, environment and remote hashes.
- `download_verification.json`: verified transfer of all retained remote artifacts.
- `termination.json` and `final_instance_status.json`: compute lifecycle.

One newly launched A100 SXM4 40 GB instance was used at $1.99/hour. At shutdown
request, its approximately 49.44-minute lifetime cost **about $1.64**, excluding
taxes or provider billing adjustments. This includes setup and unsuccessful
experiments. Only this instance was requested for termination; the two existing
A10 instances were left active. No main-thread source files were edited.
