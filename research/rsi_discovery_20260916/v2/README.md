# Spectra discovery v2

This extends the earlier discovery harness with program-level CAR/SOS discovery,
operational retained methods, typed scientific states, dependency-aware replay,
and independent acquisition/discovery/transfer campaigns. See `RESULTS.md` for
outcomes and `DERIVATIONS.md` for exact statements and limitations.

The reference is the existing Spectra algebra and solver, including its existing
spin-twirl caching and paired anticommutator optimization. Neither existing
optimization is counted as a new discovery. A small local ladder oracle checks
CAR semantics independently; it is not a molecular eigensolver.

## Reproduction

Run from the isolated repository root. Use the pinned `.venv-discovery-local`
environment (`requirements.lock`) and the current native subscription client.
The local link points to `/Users/aidenlippert/.local/share/spectra-rsi/runtime-20260916`.
It replaced execution through the original `.venv-discovery` after macOS
offloaded dependency files in Documents. Both environments are preserved.
Set `PYTHONPYCACHEPREFIX=/private/tmp/spectra-rsi-pycache` when running here to
avoid attempts to load offloaded project bytecode caches.
The model backend uses the app-bundled Codex executable if present. The older
Homebrew client could not use GPT-6-Astra; its failed attempts remain recorded.
No model is allowed to use tools during a scientific proposal.

```sh
.venv-discovery-local/bin/python -m pytest -q research/rsi_discovery_20260916/test_discovery.py research/rsi_discovery_20260916/v2
.venv-discovery-local/bin/python -m research.rsi_discovery_20260916.v2.acquire /absolute/new/acquisition
.venv-discovery-local/bin/python -m research.rsi_discovery_20260916.v2.campaign /absolute/new/campaign --policy-archive results/rsi_discovery_v2_20260916/policy_archive.json
.venv-discovery-local/bin/python -m research.rsi_discovery_20260916.v2.quotient /absolute/new/quotient
.venv-discovery-local/bin/python -m research.rsi_discovery_20260916.v2.domain_transfer results/rsi_discovery_v2_20260916/confirmatory001 /absolute/new/conformance
.venv-discovery-local/bin/python -m research.rsi_discovery_20260916.v2.compound_campaign results/rsi_discovery_v2_20260916/confirmatory001 results/rsi_discovery_v2_20260916/later_process_conformance001 /absolute/new/compound-study
.venv-discovery-local/bin/python -m research.rsi_discovery_20260916.v2.frontier_search campaign results/rsi_discovery_v2_20260916/frontier003/inputs /absolute/new/frontier-search
.venv-discovery-local/bin/python -m research.rsi_discovery_20260916.v2.correlation_search campaign /absolute/new/correlation-search results/rsi_discovery_v2_20260916/pilot001/proposal.json results/rsi_discovery_v2_20260916/acquisition003/proposal.json
.venv-discovery-local/bin/python -m research.rsi_discovery_20260916.v2.summarize results/rsi_discovery_v2_20260916
```

Outputs must use new directories. Prior campaigns are never overwritten. The
four-arm run uses five independent acquisitions, two proposal calls per arm,
150-second per-call deadlines, serial scientific evaluations and frozen transfer
methods. Native inference may run concurrently; each call has an independent
ledger/process for CPU accounting. All transfer inference has finished before
serial final timing starts. Candidate execution has wall/CPU/file limits;
macOS does not provide the Linux address-space limit used by the runner. RSS
receipts are process-family high-water marks, not additive per-attempt peaks.

The corrected polynomial-map follow-up has six fresh discovery campaigns
conditional on one shared acquired method, two proposals per arm and 180-second
proposal/evaluation deadlines. Every arm has the existing ordinary monomial-map
interface. Acquired arms use its learned backend. A second transfer reruns every
frozen algorithm with the same ordinary backend, so inherited primitive speed
cannot by itself count as improved discovery. The first failed study and the
invalidated follow-up remain in the evidence and accounting. The discarded
design had no development evaluations or transfer outcomes.

`workflow.py` runs fresh integrals -> product-initialized charge-MPS discovery ->
exact upper check -> direct collective/local coefficient maps -> numerical SOS
proposal -> original exact full-N lower acceptance. It accepts an optional frozen
molecular specification and fixed bond/solver budgets. Its generated-map adapter
checks every actual solver map before use. Independent accepting processes use
Python `-S` and refuse numerical-package imports.

`correlation_search.py` proposes executable choices of overlapping local and
coupled collective supports directly from original Hamiltonian coefficients.
It selects supports before any cubic map exists, with the same support-count
budget as the window baseline. The retained map compiler then builds the actual
solver maps. The unchanged program transfers from fresh H4 to fresh H6; original
full-N acceptance reports actual intervals and all workflow-stage costs. This
small physical exploration does not supply independent RSI significance.

`frontier.py` reads only a rational fixture and MPS from a specified directory.
It computes an algebraic observable reduction and exact dynamical defect bounds
for the frozen original four-phase control task. A loose enclosure remains a
conditional proof obligation. It never imports the archived reduced embeddings,
training snapshots or full-state trajectories located beside those inputs.

`frontier_search.py` lets native proposals generate different sparse observable
representations. Candidate construction, independent exact checking and fresh
certificate replay run in distinct standard-library-only processes. The checker
fixes identity, target observable, Hamiltonian, MPS, controls, uncertainty and
tolerance. It refuses non-Hermitian, inexact or out-of-domain proposals. Failed
and conditional programs are retained as high-risk branches, unavailable as
accepted callable capabilities. This exploratory branch is separate from the
controlled RSI experiment and is run after its timing work finishes.

## Module boundaries

`language.py` defines pure candidate programs and their explicit capabilities.
`model.py` records native proposals, failures and usage. `records.py` owns typed
states/actions, method manifests, all context-read dependencies and environment
fingerprints. `algebra.py` defines fixed mathematical targets. `evaluate.py`
owns acceptance and timing. `campaign.py` owns matched selection and transfer.
`ideal.py`/`quotient.py` check exact declared-ideal witnesses. `certified.py`
adds separately checked template operations. `workflow.py` adapts a retained
method to the existing domain solver. `frontier.py` preserves a conditional
control branch; it is not a replacement energy evaluator or a fake inverse-
design/experimental adapter.

`library.py` reopens acquired methods in later processes, verifies manifest and
source hashes, honors historical horizons and installs only the recorded
dependency closure. `compound.py` defines the independently checked polynomial
operator problem; `compound_campaign.py` owns its four-arm follow-up.

Executed study source files are archived beside their protocols. The dependency
snapshot records every copied existing Spectra module. The final source audit,
cost ledgers, raw proposal transcripts, exact certificates and frozen programs
are in `results/rsi_discovery_v2_20260916`. No universal-program, whole-pipeline
speed, physical-accuracy, or lifetime-payback claim follows from a kernel timing.

Method packages record source, input/output contracts, assumptions, domain,
evidence, dependencies, cost regime, failure cases and usage. Applicable methods
are automatically exposed by domain and actual acquisition state. An unverified
high-risk method is not callable as an accepted capability. Missing replay
branches are unknown. Changed environments require live evaluation; dependency
checks prevent later observations from appearing before their acquisition.

Correct outputs, research-discovery improvement, physical chemical accuracy,
whole-pipeline speed and cumulative payback are distinct gates. Passing one does
not silently count as passing another.
