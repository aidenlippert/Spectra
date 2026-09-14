# Warm Lambda research pool

The user explicitly requested API/CLI control, multiple A10 instances, and
keeping the instances running across experiment batches. This supersedes the
earlier per-batch teardown policy. **Both instances are intentionally running.**

| Worker | Instance ID | SSH address | Work completed |
|---|---|---|---|
| 1 | fab09ca8165b41d7bfc5649f413dda3f | ubuntu@129.159.32.200 | 30 LP experiments, four concurrent processes |
| 2 | b55288de9761445f8c589c21829f853c | ubuntu@146.235.200.232 | Nine SDP comparisons, three concurrent processes |

Both are `gpu_1x_a10`, California/us-west-1, 30 vCPUs and 200 GiB RAM each.
Current API rate is $1.29/hour each: **$2.58/hour combined, $61.92 per 24 hours**.
The existing $2,600 spending ceiling remains a ceiling, not a target. There is
no automatic teardown configured. The completed batches are idle; warm capacity
does not mean experiments continue when no task is running.

The first instance appeared active while the UI setup was in progress; it
matched the configured A10, region, and aiden-mac key, and SSH confirmed no
research workload before use. The second launch has an API request and response
receipt here. Both are explicitly registered in `pool.json`.

## Control without browser automation

```sh
python3 research/certificate_scaling/lambda_pool.py status
python3 research/certificate_scaling/lambda_pool.py types
python3 research/certificate_scaling/lambda_pool.py launch --name spectra-next-worker --region us-west-1
```

Launch is limited to one A10 per request and checks the current price and region
capacity. A single coordinator owns resource mutations; experiment workers use
SSH. Termination is explicit and refuses unregistered IDs. Do not terminate
these after each batch under the current user instruction.

The API credential is outside the repository at
`~/.config/spectra/lambda-api-key`, mode 0600. It is read by the local CLI,
never sent to experiment workers, printed, or placed in source archives.
The downloaded temporary credential file was removed after secure storage.
Official API reference: https://docs.lambda.ai/api/cloud

SSH uses the existing `~/.ssh/id_ed25519`. Both hosts have
`/home/ubuntu/spectra-venv/bin/python` with NumPy 2.2.6, SciPy 1.14.1,
CVXPY 1.6.5, and SCS 3.2.8. Source lives in `/home/ubuntu/spectra-direct`;
outputs in `/home/ubuntu/spectra-direct-results`. The current jobs use CPU
solvers; their presence on GPU hosts is not evidence of GPU speedup.

## Reproducibility

The LP lane executed `source_v1.tar.gz`; the comparison lane executed
`source_v2.tar.gz`. Their manifests record exact source hashes. All 91 LP output
file hashes and all 28 comparison output file hashes were checked after download.
All 39 returned lower certificates were independently replayed with exact
rational arithmetic and paired with separately validated upper witnesses.
`results/certificate_scaling/direct_intervals/summary.json` records those intervals.

The LP batch took 15.05 seconds after setup; the comparison batch took 26.55
seconds. These are batch observations, not controlled scaling or GPU benchmarks.
The two lanes ran at different times while the second host booted and installed
dependencies. Future batches can dispatch simultaneously without that setup cost.
