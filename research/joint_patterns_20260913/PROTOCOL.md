# Joint learned-pattern constraints: frozen bounded pass

Keep the previous exact molecular models unchanged: H4 has six learned
density-square patterns and H6 has ten. Every H6 trial solves this same
ten-pattern Hamiltonian, including trials that use only a subset of its
patterns to construct additional constraints. Replay the same collective
remainder certificate and transfer the new retained-model lower bound.

The baseline is the full symmetry-partitioned quadratic fermionic SOS cone
(linear, particle-particle, hole-hole, and particle-hole operators), with a
one-body multiplier of the fixed-particle-number identity. Added odd
generators are Q_k a_p and a_p. A joint block imposes positivity of
{B^dagger,B} for every linear combination B within that symmetry block.
All cross-pattern terms are kept. CAR cancellation of every degree-six
term is checked exactly; none may be truncated. Separate-pattern blocks
form an ablation with the same generators but no cross-pattern products.

Predetermined cases: H4 baseline and joint six-pattern control; H6 baseline,
separate ten-pattern ablation, and joint two-, four-, and ten-pattern cases.
The factor order is the existing decreasing coefficient-eigenvalue order;
neither ground states nor old proof factors are selection inputs.

Each discovery process gets a 90-second end-to-end wall budget. Solver time
is capped at 60 seconds, reduced after construction to reserve 20 seconds
for rational export and exact replay. A parent-process watchdog enforces
the wall budget, including imports and file writes. A timeout or solver
failure remains a recorded failure. No certificate counts as within-budget
unless its exact replay finishes inside that budget. Development, pilot
runs, fresh independent replays, and pre-existing input construction are
reported separately. Single-thread BLAS is used.

Acceptance uses rational factors, the existing CAR SOS verifier, and the
previous exact tail verifier. The two anticommutator halves share one
rounded factor and are related by exact adjunction. Save compact pattern
references and count their fully expanded proof size and replay cost.
Floating objectives and solver statuses are diagnostics, not energy proofs.

Compare accepted lower bounds using separately replayed existing rational
reference upper witnesses. Their prior wavefunction-discovery cost is not
removed by this method and is excluded from lower-bound discovery. The
1.6 mHa interval is the success target. Also report gain, factor dimensions,
coefficient-map work, construction, solve, export, replay, and total cost.
A weak or unfinished solve does not prove this constraint span insufficient;
an exact feasible dual witness would be needed for such a claim.

This is a Hamiltonian-derived compression test within established partial
three-particle positivity ideas. See Mazziotti, [Phys. Rev. A 71, 062503
(2005)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.71.062503)
and [Phys. Rev. Lett. 117, 153001 (2016)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.117.153001).
No general many-body solution, asymptotic theorem, or priority claim follows
from these two molecular fixtures.
