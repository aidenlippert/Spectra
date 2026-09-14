# V9 construction headroom protocol

Same 36 V7 development cases, physical assumptions, horizons, tolerance .001,
term/cache caps512, degreecap24. No reserved width5/7 ladder/X/Y outcomes.

Compare complete fresh calls: fraction-free adaptive Taylor; trajectory-wide
coordinate Galerkin; matched Krylov projected Taylor; direct projected ODE
collocation. Krylov arms try dimension8 then16, reusing each basis within its
order attempts. Projected Taylor tries degrees4,6,8,10,12,16,20,24; collocation
tries4,6,8,10,12. These are explicitly bounded conventional search schedules,
not optimal adaptive implementations or learned programs. Basis construction,
all failed attempts, linear solves, quantization, power/Bernstein conversion,
norm search and independent checker work count. Matrix solve dimension is capped
at384; no hidden Hilbert-space dense diagonalization is used.

For each polynomial, build Bernstein residual records once and attempt l1,
firstfit, weighted norm partitions until passing; exact checker always recomputes
accepted candidates with a fresh original V7 Generator. Failed norm bounds need
not be independently checked because they authorize no guarantee. Record them
and their costs. Structural budget failures close the dependent search branch;
other failures and all successes are retained.

The coordinate Galerkin method constructs a fixed-space trajectory, accumulates
all omitted derivative coefficients (including terminal), and grows the space
by largest accumulated weighted omission. Increase Taylor degree until either
an l1 certificate meets tolerance or the in-space tail is <=eps/4, then attempt
ordinary norm certification and grow/restart if needed. Default additions are
max(4, floor(basis_size/2)), capped by remaining support; 16 expansions maximum.
These choices are heuristics; exhaustion is a refusal, not impossibility.

Run one diagnostic pass first. Repeat paired timings only if a candidate
actually shows full-cost headroom at equal accuracy or recovers feasibility.
Any useful conventional method strengthens the baseline and is not relabeled
as acquired m1. No learner is authorized by a smaller polynomial alone.

Before promoting the one Galerkin recovery over full Taylor, add a stronger
BFS comparison through the same incremental projected recurrence and witness
code. Its frontier grows by complete BFS layers, retaining the same fixed
support/cache budget and the same tail-based degree choice; every discovery
column counts. Replay all five arms. This isolates basis selection from the
previous coarse BFS degree schedule and redundant witness reconstruction.
