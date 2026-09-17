# Direct compact construction for the supplied H8 control task

The active objective is to remove full-sector trajectory discovery and checking
from the supplied finite-control calculation. A short final reduced evolution
does not establish success if either construction or acceptance first expands
the balanced sector. No universal many-body solution is presumed.

The input is the exact supplied canonical H8 Hamiltonian and rational charge MPS from
`Spectra_certified_control1.zip`, imported without changing its bytes. It has
16 spin orbitals, 8 electrons, spin counts (4,4), and an initial MPS with maximum
bond dimension 144. The operator controls and target D <= -0.6 are taken as
task data from the archive. The short T=2 schedule is the first primary test;
the longer T=21 constant schedule is a later stress test if measured costs allow.
The original reduced embeddings and full-model trajectories may be inspected
as explicitly charged references; they are not construction inputs for a claimed
independent compact method.

The primary routes are direct tensor-state dynamics, task-specific observable
closure, and compact exact contraction/certification. Construction must use the
Hamiltonian coefficients and the supplied compact initial state, not a supplied
list of all 4,900 amplitudes. Every operator term is represented or charged by a
rigorous remainder allowance. Local small-system oracles and a replay of the
imported enumerated reference are separate validation costs.

Compare the requested observable target and uncertainty budgets from the source
receipts. Do not silently replace them with a weaker task. Optional shorter-time
or relaxed-tolerance tests are diagnostics only. Source control coefficients are
mathematical actuators; no laboratory preparation or physical-model guarantee
is assumed.

Record source hashes, operator and state bond dimensions, tensor entries,
construction, contractions, failed attempts, peak process memory, and replay.
Use one heavy process at a time on the local 8 GB machine. Numerical libraries
can propose representations; acceptance must include their actual errors. No
paid compute, GitHub operations, or external messages are needed for this pass.

The observable-only strong-control experiment is a separately labeled regime.
It retains the exact H8 Hamiltonian, initial state, observable, and nonzero
uncertainty budgets, but allows a pulse amplitude of 64 Ha. It does not replace
or satisfy the primary frozen amplitude-.5, four-phase task. Its positive result
and the original task's unresolved status must both appear in the final report.
