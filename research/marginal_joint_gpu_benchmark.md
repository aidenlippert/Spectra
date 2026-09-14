# Joint projector CPU/GPU benchmark

The bounded harness is
`results/marginal_graded_hubbard8/discovery/joint_gpu_benchmark.py`.
It prepares the actual94 local reflection/spin blocks, all32 joint Gram
blocks and1,280 Gram columns. Float64 matrices use physical reflection
norms and symmetric weighted Gram normalization. Profile batches retain
per-sector derivative tensors in memory; ragged block shapes are handled
separately. Nonzero half and charged penalties enter the density objective
with the correct sign and factor1/5.

Reproduce the final local CPU receipt from the Spectra root:

```sh
OPENBLAS_NUM_THREADS=1 python -m results.marginal_graded_hubbard8.discovery.joint_gpu_benchmark --evaluations 6 --batch 3 --ratio 0.5
```

Observed timings in `results/marginal_graded_hubbard8/gpu_benchmark/receipt.json`:

* Matrix/Gram preparation:11.460854 seconds.
* Six batched profile evaluations:0.113940 seconds.
* Scalar SciPy smallest-eigenvalue reference:0.105665 seconds.
* Gram largest-eigenvalue solves:0.005810 seconds.
* Batched versus scalar maximum discrepancy:2.54e-14 or smaller.

No GPU ran: CuPy is absent locally. There is no measured GPU speedup.
The initial incomplete harness receipt is retained separately and should
not be used for a GPU comparison. In this small run, preparation dominates
end-to-end numerical time; much larger resident batches must be measured
before drawing acceleration conclusions. The exact rational proof replays
remain separate CPU acceptance gates.

The GPU path transfers resident matrices once, warms a batched eigensolve,
times profile and Gram kernels separately, synchronizes CUDA timing points,
and records copy time, resident bytes and matrix-batch size. Vendor solver
workspace is explicitly unmeasured. It refuses CPU/GPU eigenvalue mismatches
beyond1e-9 absolute plus1e-10 relative tolerance. Computation failures are
not relabeled as missing CuPy. Limits are64 evaluations and16 per batch.
The GPU execution path itself remains untested until a CUDA device is used.

On an existing Lambda instance, activate a CUDA-compatible CuPy environment,
confirm a visible CUDA device, copy the source tree, and run the same command.
Use `--evaluations 64 --batch 16` for the bounded larger comparison. Keep the
JSON receipt and compare end-to-end timings as well as resident kernel time.
The equivalent exact certificates must still pass the standard-library
CPU replay before any scientific result is accepted.

Lambda dashboard/key discovery is recorded in `migration/GPU_ACCESS.md`.
The registered aiden-mac key exists, but no instance was running at that
inspection. No paid resource was launched by this task. An instance launch
scope is being coordinated in the old task; the available credit balance
is not a blanket spending authorization.
