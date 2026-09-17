# Spectra discovery loops

This is a working, bounded discovery system in the isolated `codex/rsi-discovery`
worktree. It proposes and executes certificate constructors, checks conjectures
and proofs, retains methods, learns exploration policies from recorded outcomes,
and reopens those methods in a separate future campaign. It does not establish
a general many-body solver, a new theorem, scientific novelty, or a complete
runtime advantage.

The first executable campaign uses small molecular models and the existing
enumerated reference certificate checker. It is an integration and research
control, not an integration of every existing Spectra research branch. Hydrogen
chains are two input families; the system also runs non-chain water models.

See [RESULTS.md](RESULTS.md) for measured outcomes, failures, costs and limits.

## Run

From `/Users/aidenlippert/Documents/Spectra-rsi-discovery`:

```sh
.venv-discovery/bin/python -m pytest -q research/rsi_discovery_20260916/test_discovery.py

# Default: local executable mutations, no model call.
.venv-discovery/bin/python -B -m research.rsi_discovery_20260916.campaign \
  results/rsi_discovery_20260916/new_campaign

# One proposal through the supported Codex CLI, using its existing login.
.venv-discovery/bin/python -B -m research.rsi_discovery_20260916.campaign \
  results/rsi_discovery_20260916/new_model_campaign --provider codex

# Separate fresh campaign reopening an acquired method library.
.venv-discovery/bin/python -B -m research.rsi_discovery_20260916.campaign \
  results/rsi_discovery_20260916/new_future_campaign \
  --reuse-from results/rsi_discovery_20260916/campaign_002
```

Output directories must be new. Old witnesses and failed attempts are preserved.
The default run has four predeclared molecular tasks, one worker, a maximum of
one optional model proposal, and fixed stage limits. The separate reuse run
uses two additional frozen geometries. Nothing schedules unattended runs,
pushes to GitHub, or purchases compute.

The clean `.venv-discovery` environment passes `pip check`; exact versions are
in [requirements.lock](requirements.lock). Python 3.12 was used. To recreate it,
create a standard isolated virtual environment and install that lock file.
The earlier `.venv-rsi` is retained but is not the campaign runtime. The proof
run used the installed Lean 4.30.0-rc2 toolchain, with no Mathlib dependency.
If Lean is unavailable or a proof fails, the statement remains a conjecture
and the acquired proof capability is not enabled.

## Three connected loops

1. **Algorithms.** Restricted Python programs construct their own numerical
   eigenvector and rounded Cholesky witnesses. Mutations change eigensolver
   strategy, factor precision and margins. The optional model supplies another
   executable program. Every candidate is run in a separate process; the
   independent checker reconstructs and validates the resulting witness.
2. **Lemmas.** Z3 searches for a counterexample to proposed cross-term bounds.
   Its rational witness is replayed exactly. CAR verifies a finite occupation
   projector identity. Lean checks two reusable, explicitly integer arithmetic
   statements. Failed proofs stay unproved; neither a Z3 answer nor a handful
   of examples becomes a general theorem. The checked cross-term lemma and
   benchmarked integer Gram component enable an adaptive remainder constructor.
3. **Research policy.** Recorded trees expose only reachable action metadata
   and already-observed outcomes. A finite search varies branching width,
   depth/continuation preference, stopping patience and family priority. Fixed,
   simple adaptive and selected policies are then executed on fresh histories
   with the same proposal tree and resource envelope. Replay cannot answer an
   unrecorded branch. New branch outcomes must actually be evaluated.

The retained library stores source and evidence by hash. Reuse reopens the
library, audits its history, checks evidence bindings, reruns the Lean templates,
and independently checks every new scientific certificate. Prior success grants
permission to try a method; it never grants acceptance on a new task.

## Independent acceptance

For every exact rational Hamiltonian, the evaluator reconstructs all blocks of
the entire fixed-number sector, with no omitted spin-projection sector. It
uses the existing `enumerated_baseline.blocks`, `shifted`, and `quotient`
functions and `compact_response.closure.check_factor`. Each lower witness
proves that its shifted matrix is a rational Gram matrix plus an exactly
diagonally dominant nonnegative remainder. A nonzero integer Rayleigh vector
provides the upper endpoint. The global interval must be at most **1.6 mHa**.

The independent check runs with Python's site packages disabled and rejects
numerical imports. It does not execute candidate code. It checks task identity,
complete blocks, exact integer factors, positive denominators, Hermiticity,
particle/spin-count conservation, positivity, Rayleigh quotients and the target.
Candidate output flags and reported bounds cannot replace this reconstruction.
The original checker remains the acceptance authority after acquisition.

The exact component task computes the integer Gram matrix of a triangular
factor. Its outputs must equal the independent Python-integer implementation
on actual molecular factors and zero, signed and large-integer stress cases.
The FLINT implementation is benchmarked using three executions per input.
This is finite program equivalence evidence, not a universal proof of arbitrary
code. New molecular certificates still receive the original independent check.

Lean checks `cross_term_bounds` and `compose_lower_errors` over **integers**.
Their axiom reports permit only Lean's standard foundational axioms
`propext`, `Quot.sound`, and `Classical.choice`. No `sorryAx` or custom axiom is
accepted. These familiar inequalities are useful proof components, not novelty
claims. No new PSD acceptance rule, quotient equivalence, or physical-model
guarantee is silently introduced.

## Evidence and costs

Each output directory contains:

- `protocol.json`: frozen tasks, split, resource limits and complete source map;
- `history.sqlite`: append-only, hash-chained events including started operations;
- `objects/`: content-addressed source, prompts, task inputs, certificates,
  checker receipts and provider event logs;
- `attempts/`: executable candidates, requests, results and independent-check logs;
- Lean source and axiom reports;
- `report.json`: attempts, measured costs, policy comparisons and reuse accounting.

An interrupted operation remains visible. History and object hashes are audited.
This detects accidental alteration; it is not a digital-signature trust service
against a person who controls the entire local filesystem. Candidate programs
cannot write the history or change evaluation rules.

All automated stages in a run are charged, including failures, input generation,
model proposals, exact checking and policy search. The total wall clock includes
controller overhead; summed stage clocks exclude it. Candidate body timings,
whole attempt clocks and total campaign clocks are separate quantities. Child
RSS is labeled as the process-family high-water mark, not falsely summed or
presented as a per-attempt peak. Model token usage is recorded when supplied by
the provider. Subscription marginal price is unknown, not zero.

Earlier setup, manual implementation time/tokens and tests outside each campaign
are explicitly outside its meter. Accordingly, this is not a fully measured
lifetime research bill. The results report includes the initial failed run and
does not erase its acquisition cost. A failed or unmatched arm cannot yield a
cost-advantage number. No reward is based on output-file size, an unchecked
numerical objective, or test-only success.

Policies reserve bounded worker time, stop according to their observed history,
and charge any verification overshoot. Stages have wall and CPU limits;
candidate output files are limited. On macOS there is no enforced hard address
space limit. The candidate language disallows imports, dynamic execution,
introspection and host I/O. It is a restricted language, not a sandbox for
arbitrary Python packages. The supported Codex process retains its normal
read-only tool sandbox and is exempt from the candidate output-file limit.

## Integration coverage

| Area | Actual status in this change |
|---|---|
| Orbital representation | Fresh RHF canonical orbitals via the existing generator; no orbital-representation search |
| Correlated-state/MPS construction | Not integrated into this campaign; no MPS improvement claim |
| Local/collective lower-bound dictionaries | Existing module is snapshotted as a future seam; no dictionary search or large SOS pipeline run |
| Fermionic, spin and number algebra | Existing exact CAR and complete fixed-N/alpha-count blocks are exercised; finite projector identity checked |
| Direct coefficient-map construction | Exact Hamiltonian-to-block construction is used; no new scalable SOS coefficient map |
| Compact response/recursive elimination | Existing remainder checker is reused; no new response or elimination method |
| Rigorous remainders | Exact signed-row remainder checks and checked integer arithmetic lemmas; acceptance unchanged |
| Numerical optimization | Executable eigensolver/Cholesky/precision variants and exact acceptance; no new variational family |
| Exact dual obstructions | Not implemented; physical nullspace versus truncated ideal equivalence remains unresolved |
| Dynamics/observable closure | No evaluator or discovery campaign implemented |
| Robust control | No control evaluator; actuator limits are not relaxed or replaced |
| Hardware kernels/scheduling | Exact CPU integer Gram kernel and serial bounded scheduling; no GPU benchmark |
| Research-policy selection | Causal historical replay, finite policy evolution and real fresh policy runs |
| Chemical-model selection, other physical queries, inverse design, experimental planning | Separate scientific evaluators are required and are not implemented or validated |

The molecular guarantee covers the exported rational finite Hamiltonian.
Floating-point integral evaluation, physical-model error, experimental accuracy,
and larger-system scalability are not certified. Enumeration is explicit and
limited to 2,000 sector labels and at most 16 spin orbitals; this is not an
enumeration-free solution to the existing H8/H10/H12 research problem.

## Provenance and related implementations

[dependency_snapshot.json](dependency_snapshot.json) records 45 copied existing
modules, absolute source paths, byte counts and SHA-256 hashes. Their originals
were only read. Source hashes are checked at campaign start and before final
reporting. The worktree is based on commit
`695767960e47bc69607234071d35d53ed940cab7`; ongoing main-checkout research is not
modified. Git's initial full checkout and later broad status reads stalled in
local object traversal, so this task restored its index and deliberately
materialized the required source subset. Full repository checkout is not claimed.

The supplied Dream-RSI paper was read from
`/Users/aidenlippert/Downloads/2609.14858v1.pdf`. Its
[official repository](https://github.com/zhengkid/Dream-RSI) still marked its full
implementation as forthcoming when checked. This implementation adopts the
recorded-history replay idea; it does not purport to reproduce unreleased code
or the paper's reported results. Historical improvement is not a fresh-task
guarantee.

[ShinkaEvolve](https://github.com/SakanaAI/ShinkaEvolve) 0.0.7 was assessed using
its installed headless Codex route parser and actual full-patch engine. The
smoke check confirmed program replacement while preserving text outside the
evolution region. Its default headless route uses a separate Node wrapper;
the current runtime instead calls the already-authenticated native Codex CLI
directly. No full Shinka evolution campaign was run, and OpenEvolve was not
installed. No authentication secrets were extracted or treated as API keys.
