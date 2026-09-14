# Verified Lambda GPU handoff

The user authorized responsible use up to $2,600 total and explicitly approved Lambda's final license/billing agreement. They requested exactly one final handoff to the continuation task, then archive the old task and stop all messages in both directions. Do not reply to the archived task or route subagent reports there.

## Running instance

- Name: spectra-gpu-benchmark-20260912
- Type: NVIDIA A100-SXM4-40GB, 30 vCPUs, displayed $1.99/hour
- Region: us-west-2 (Arizona)
- SSH: `ssh -i /Users/aidenlippert/.ssh/id_ed25519 ubuntu@129.146.98.55`
- Registered key: aiden-mac. Private/public key match verified; SSH login successful.
- Dashboard: https://cloud.lambda.ai/workspace/4c315b664b1444fa8f10eb579e0c0482/instances
- Python: `/home/ubuntu/spectra-gpu-venv/bin/python` (3.10.12)
- Remote source/data: `/home/ubuntu/spectra-gpu-20260912`
- Environment: `CUDA_PATH=/usr OPENBLAS_NUM_THREADS=1`.
- Pinned packages: NumPy 2.2.6, SciPy 1.14.1, CuPy CUDA12x 14.2.0, cuda-pathfinder 1.8.1; driver 580.105.08; CUDA 12.8.
- No API key created. Browser account is already signed in. Instance termination must use Lambda; guest shutdown is not assumed to end billing.

## Measured result

Use the grouped resident evaluator, `resident_batch_benchmark.py:evaluate_groups`, when batching independent numerical candidates. It groups all 94 actual physical sectors into 12 matrix sizes, assembles on GPU with float64, and returns minima. Each candidate includes all 94 sectors.

| Candidates | Batch | CPU scalar reference | GPU median including result copy | Speedup |
|---|---|---|---|---|
|256|128|3.4976 s|0.9058 s|3.86x|
|1024|256|14.0788 s|3.5225 s|4.00x|

CPU is single-threaded SciPy smallest-eigenvalue reference on the same instance; these are not comparisons against an optimized multicore CPU search. CPU grouped implementation was slower than scalar and both agree. All three GPU repetitions passed a 1e-9 absolute plus 1e-10 relative tolerance; maximum observed error was 7.98e-13. Floating-point agreement does not certify positivity or a mathematical bound.

The original ungrouped 64-candidate port was slower on GPU (1.7917s versus 1.1298s CPU), and individual Gram eigensolves were slower too. Do not route small or sequential solves blindly to the GPU.

Startup matters: first grouped warmup/compilation was 14.75s. The subsequent 1024-candidate process had 1.14s warmup and 4.84s total load/transfer/warmup/first search, excluding generation of physical matrices and the CPU comparison. Reuse matrices and warmed processes for sustained search. Mac physical matrix preparation alone took 10.43s in the original harness; a separate preparation with two additional direct-assembly checks took 20.20s. Resident matrix data is about 24 MB; peak constructed matrix batch about 369 MB, excluding solver workspace.

## Reproduce

From the Mac, run `run_remote.sh` in this directory. Or on the instance:

```sh
cd /home/ubuntu/spectra-gpu-20260912
CUDA_PATH=/usr OPENBLAS_NUM_THREADS=1 timeout 300 /home/ubuntu/spectra-gpu-venv/bin/python resident_batch_benchmark.py --data resident.npz --out resident_recheck.json --evaluations 1024 --batch 256
```

The existing data is frozen around the earlier refined certificate, not the latest scientific optimum. Rebuild the resident tensors from current profiles/projectors when integrating into current research. Preserve orbit normalization, both-sided weighted Gram normalization, and objective direction. Existing exact integer/rational CPU replay remains the acceptance gate. Current timings justify this batched numerical path only; they make no general chemistry/scalability claim.

## Receipts and ownership

Local directory: `/Users/aidenlippert/Documents/Spectra/results/lambda_runs/20260912T092845Z`.

- `launch_receipt.json`: authorization, instance, source archive hashes.
- `a100_baseline.json`, `mac_baseline/receipt.json`: original same-source benchmarks.
- `resident_a100_256.json`, `resident_a100_1024.json`: repeated GPU measurements and complete output comparisons.
- `resident_preparation.json`: direct physical matrix reconstruction discrepancy 1.43e-14 or less.
- `resident_batch_benchmark.py`, `prepare_resident.py`, `resident.npz`, `run_remote.sh`: reusable numerical artifacts.

All downloaded receipt data/source hashes matched the local frozen files. Grouped CPU evaluator additionally matched independent scalar assembly at batch sizes 1, 3, and 7, including incomplete batches, for all 94 physical sectors.

The continuation task owns this running resource after the one final handoff. Start useful bounded work promptly, record elapsed spend at $1.99/hour, copy results back, and terminate the instance in Lambda when no useful job remains. $2,600 is a ceiling, not a target. Do not provision additional GPUs without evidence they improve the useful work per dollar. The old task will be archived and will not monitor, maintain, or respond about this instance.
