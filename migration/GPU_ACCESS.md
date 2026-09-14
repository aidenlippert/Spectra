# Current warm A10 pool — supersedes historical termination notes

The user explicitly requested multiple A10 hosts, API/CLI access, and keeping
instances running between batches. Both are intentionally ACTIVE:
- worker 1: fab09ca8165b41d7bfc5649f413dda3f, ubuntu@129.159.32.200
- worker 2: b55288de9761445f8c589c21829f853c, ubuntu@146.235.200.232
Each $1.29/h; combined $2.58/h. Current cap remains $2,600 total.
Do not tear down after each batch. Do not message archived tasks.
CLI: `/opt/homebrew/Caskroom/miniconda/base/bin/python research/certificate_scaling/lambda_pool.py status`.
API key stored outside repo at ~/.config/spectra/lambda-api-key (0600).
SSH key unchanged. Remote venv /home/ubuntu/spectra-venv; source
/home/ubuntu/spectra-direct. The later adaptive campaign used spectra-adaptive
and spectra-adaptive-v2; all32 remote jobs finished, downloaded and exactly
replayed. No research process remains live; both hosts intentionally warm.
Current result: research/certificate_scaling/STRUCTURAL_ATTACK.md.
The later structural campaign ran24jobs (23certificates,1timeout), all downloaded;
no researchprocess remains live. Active structuralgoal continues inthis task.
Full ownership, commands, source manifests, and lifecycle:
results/lambda_runs/direct_discovery_pool/README.md.

---

# Latest GPU campaign — 2026-09-12 certificate scaling

The current task launched one A100-SXM4-40GB at $1.99/hour under the recorded
$2,600 total authorization, ran ten bounded experiment jobs, retrieved results,
verified69 remote file hashes, and independently replayed37 returned certificates.
The instance spectra-certificate-scaling-20260912 (158.101.19.250) has been
TERMINATED through Lambda. Provider UI verified: No running instances.
Observed running interval estimate about$0.29, not an invoice; excludes the
unobserved interval before the first running-state observation. All receipts,
used source and results archives: results/lambda_runs/certificate_scaling/.
The old IP is not active. SSH key remains the existing aiden-mac key. No API key
created. No agents or tasks should contact the archived historical task.
See research/certificate_scaling/RESULTS.md for authoritative findings.

## Historical status follows

# GPU status — latest continuation

The user approved responsible Lambda spending up to $2,600 total. The one
A100 instance handed to this task was used for785 current-model scientific
candidates, results were retrieved and hashes compared, then it was
terminated through Lambda. Dashboard verified: No running instances.
Observed-lifetime cost estimate including the earlier benchmark: about$0.72
at$1.99/hour; not a final invoice. See
results/lambda_runs/diagonal_science/lifecycle.json and
research/marginal_gpu_diagonal_search.md.

The SSH key remains available for a future authorized bounded run; the old
IP129.146.98.55 is no longer an active instance. The old task is archived and
MUST NOT receive messages or reports. Future compute belongs to this task.
Detailed now-historical instance/environment instructions are in
results/lambda_runs/20260912T092845Z/GPU_HANDOFF.md.

## Historical access discovery

# Existing Lambda access — forwarded 2026-09-12

The old task inspected the signed-in Chrome Lambda dashboard at cloud.lambda.ai.
Workspace: 4c315b664b1444fa8f10eb579e0c0482. Registered key aiden-mac exactly
matches /Users/aidenlippert/.ssh/id_ed25519.pub; private key path is
/Users/aidenlippert/.ssh/id_ed25519. No credential content is stored here.
Fingerprint: SHA256:0ZfVHYg8EKcjuIWfV0Mby6n/E5ubYUgbOuq/TqBcOJ0.

The dashboard showed no running instances; no launch occurred. Historical
ubuntu@68.209.75.199 timed out before authentication and is not known active.
No API key was found or created. User has about $2,600 of credit and requests
GPU acceleration, but a bounded paid instance launch scope is still pending
in the coordinating old task. Do not interpret the balance as a spend budget.
The benchmark work remains local until an active instance is available.
