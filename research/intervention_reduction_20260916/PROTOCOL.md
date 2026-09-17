# Molecular intervention reduction — frozen initial experiment

The task is to construct a reusable, checked reduced dynamical model on existing
interacting molecular Hamiltonians. This is a test of the larger many-body
objective; no universal algorithm or major breakthrough is presumed.

The initial comparison uses the unchanged H6 asymmetric, H8 cold, and asymmetric
water inputs from `results/transfer_solver_20260915/cases`. H4 was planned as an
integration control but was not run in this pass; exact two-level controls were
used for checker tests instead. The two controls are a spin-summed population difference and hopping
between the highest initially occupied spatial orbital and its next neighbor.
They do not commute. Both control amplitudes obey |u| <= 1/100 hartree.
The initial requested horizon is 10 atomic time units. The target is a uniform
state-vector error <= 1/200, hence population expectation error <= 1/100.

Discovery uses an existing MPS only to select configurations by a deterministic
prefix beam. The initial candidate configuration budgets are 64, 256, and 512;
reduced dimensions are 4, 8, and 16 where available. A selected configuration
space may exhaust a small physical sector; this must be reported. MPS discovery
is an inherited dependency, and no cold end-to-end speedup is claimed.

The accepting calculation reconstructs Hamiltonian and control actions on every
selected configuration from the exact CAR input. All reached external labels
are counted. It verifies the positive reduced metric, projected generators, and
the omitted action at the control-box vertices. No omitted action may be dropped
because its coefficient is small. Integer candidate vectors are defining data,
not claimed exact eigenvectors. The initial state is that explicitly specified
normalized vector, not an unknown ground state or an experimentally prepared
state.

Compare a conventional uniform residual bound with a componentwise transition
envelope that removes real diagonal phase evolution. Count construction,
unsuccessful candidates, replay, selected/reached configurations, storage, and
memory. This method is a conventional projected-subspace/residual construction;
a new implementation or a favorable result does not establish field novelty.

Numerical propagation is an independent diagnostic, never the accepting proof.
Exact tests include singular metrics, omitted leakage, noncommuting control
cross terms, nonorthogonal coordinates, malformed sectors, and phase-only
motion. No Git operations, paid compute, or experimental apparatus is required.

## Recorded adaptation after the first H6 failures

The initially selected low-energy eigenvectors failed the uniform 10-au target.
Response-weighted eigenvector selection also failed, including rank 64. This is
a bounded-search failure, not a theorem ruling out small reduced models.

The next representation is ordinary snapshot/POD reduction on the selected
configuration matrix. Training uses constant controls at the four amplitude-box
corners, at times 2.5, 5, 7.5, and 10 au. The first basis column remains the
selected-space ground eigenvector. Additional real columns come from real and
imaginary snapshot deviations. The three test protocols in `predict.POLICIES`
were not training trajectories. A direct full-space polynomial residual checker
charges numerical propagation, rational rounding, all segment jumps, and
normalization. This preserves cancellations discarded by the uniform envelope.

The H6 rank-16 two-stage test met the original 0.005 state-error target. Freeze
this snapshot rule, a 512-configuration beam, rank 16, order 10, and 0.25-au
segments for the first water and H8 transfer attempts. Their outcomes remain
unseen at this freeze. The all-waveform uniform claim is not substituted for
the narrower protocol-specific acceptance. Snapshot construction can require
the complete selected-sector eigensystem; record this dependency explicitly.

## Subsequent adaptations and new tests

The frozen rank-16 rule failed on water and H8. Residual/POD enrichment then
added eight directions per round using the same four training corners and times
0, 2.5, 5, 7.5, 10. Existing columns, especially the initial state, were preserved
exactly. Water reached the switch target at rank 32; H8 did not. At water rank 32,
the original fixed-step Taylor proposer became unstable and the checker returned
a huge bound. Its receipt is retained. Order-18 Taylor with a generator-dependent
step repaired time propagation. This is an adaptation, not a pass of the frozen
original time rule.

The existing H6 rank-16 model passed the two additional predeclared test pulses.
Exact Bernstein norm enclosures extended endpoint certificates to the whole
trajectory. Exact control-norm bounds extended each accepted pulse to all
measurable perturbations within specified integrated-deviation budgets and the
original amplitude box. The first standard budget corresponds to a pointwise
deviation of at most 1/400000 Ha per control over 10 au.

Before pulse search, the declared H6 inverse objective was a population increase
of at least 1/40 in the selected spatial orbital, with the same 10-au horizon,
0.01-Ha amplitude limits, and 1/400000-Ha perturbation allowance. The finite menu
contains 729 three-stage schedules, of durations 3, 4, 3 au, with each control
at its negative maximum, zero, or positive maximum. The top numerical candidate
was checked exactly; neither a global control optimum nor laboratory
reachability is claimed. It met the functional target at rank 16 but missed
the stricter 0.005 state-error target. Rank-24 enrichment preserved the chosen
pulse and initial state and met both. A later, explicitly larger 1/2000-Ha
perturbation allowance also preserved the functional target, while failing the
stricter state-error target. These two success criteria are not interchangeable.

Compilation of the exact action Gram matrix was introduced after profiling
repeated trajectory replays. Query-only timings exclude its charged original
CAR reconstruction. The compiled query path discards the configuration records;
the construction does not avoid them. Direct and compiled arithmetic must agree
on every mathematical receipt field, not merely on the pass/fail flag.

After water rank 32 passed the switch pulse, two other predeclared pulses were
tested without changing that basis. The interior pulse passed; the three-stage
pulse missed. A bounded rank-40 residual enrichment is allowed as another
adaptation, using unchanged training corners rather than the failed test pulse.
Chebyshev interpolation with exact integer conversion to power coefficients is
also tested as a numerical proposer to reduce the many time segments required
by stable Taylor propagation. It changes neither the operator model nor the
accepting residual theorem. All old candidates and unsuccessful costs remain.
