# V9: trajectory-wide construction and direct ODE collocation

**The scientific-method acquisition goal remains unachieved.** V9 implements
and checks two materially different approximation-construction methods, then
rejects their present schedules on complete cost against stronger conventional
methods. No m1/m2 learner or heldout evaluation was run. V1–V8 remain preserved.

[Protocol](../../research/v9/HEADROOM_PROTOCOL.md) ·
[180 full calculation rows](headroom.json) ·
[mathematical constructions](../../research/v9/CONSTRUCTIONS.md).

## Results under fixed physical requirements

All arms use the same 36 development generators, central-Z observables,
gamma0/.2/2, horizons .2/.5, tolerance .001, live support/cache512 and degree24.
The collocation stage solve has an additional explicit dimension cap384.
Construction, failed rounds, norm proposals, conversions, solves, accepted
certificate checking and refusal work are counted. This is one diagnostic
pass, not a statistical confidence claim about timing or generalization.

| Conventional arm | Certified | Refused | Total seconds including failures |
|---|---:|---:|---:|
| Fraction-free adaptive Taylor | 26/36 | 10/36 | 2.176 |
| Trajectory-wide residual enrichment | 27/36 | 9/36 | 7.769 |
| BFS through the same adaptive recurrence | 27/36 | 9/36 | 4.240 |
| Matched projected Krylov/Taylor | 24/36 | 12/36 | 8.369 |
| Direct projected ODE collocation | 24/36 | 12/36 | 6.982 |

Residual enrichment is slower on every one of the 26 cases shared with full
Taylor. It recovers n6 XXZ, gamma2, T=.2, but that case was already solvable by
BFS. Through the strengthened common implementation, BFS takes about .147s
versus enrichment's .281s, despite keeping213 rather than181 basis labels.
The larger retained basis is cheaper to discover and certify here.

Collocation often reaches tolerance at lower polynomial degree than the
projected Taylor arm, but does not beat full integer Taylor on any common
successful case. Its projected solves, failed degree attempts and exact
residual conversion cost more than the shorter polynomial saves. The bounded
schedules are not claimed to be optimal implementations of every Galerkin,
Krylov or collocation method. They do not pass the acquisition headroom gate.

## What changed mathematically and computationally

The previous residual heuristic chose labels from individual Taylor steps.
V9 constructs a fixed-space trajectory, includes every omitted derivative
coefficient—including the terminal one—and chooses its next expansion using
accumulated residual contributions. Degree advancement and basis expansion are
separate decisions. Every accepted result is replayed with a fresh original
checker; proposer caches are not trusted.

The polynomial method solves derivative-collocation equations in a numerical
Krylov space. It does not interpolate already-computed trajectories. Its
floating stages become a rational polynomial whose initial error and full
residual are checked exactly. This avoids an unnecessary proposed requirement
for separate node-value certificates or algebraic Chebyshev nodes.

Both are supplied conventional methods. Improving how the implementation
constructs a calculation does not show that the system acquired that method.

## Incremental reuse diagnostic

A separate exact replay inspects four development trajectories. Reusing an
unchanged prefix would avoid approximately6–19% of generator multiply-adds,
while retaining up to3296 raw derivative entries. Difference propagation
reduces more generator work but adds substantial coefficient subtraction and
map merging:

| Case | Original generator multiply-adds | Difference-action multiply-adds | Additional merge visits |
|---|---:|---:|---:|
| n3 XXZ, gamma0, T=.2 | 590 | 397 | 472 |
| n3 XXZ, gamma2, T=.5 | 4974 | 4082 | 3434 |
| n4 mixed, gamma2, T=.5 | 128477 | 98952 | 52172 |
| n6 XXZ, gamma2, T=.2 | 22216 | 15007 | 13953 |

These counters are different operations, not interchangeable runtime units.
There is no demonstrated full-cost win. Fresh certificate checking remains
necessary. [Exact diagnostic](reuse_probe.json) and
[prefix/difference lemma](../../research/v9/incremental_basis_lemma.md).

## Verification and corrections

All128 accepted evolution certificates were exported, reparsed and replayed
through the unchanged V7 checker. All169 tests pass. Historical source/data
hashes and both V7/V8 certificate archives remain unchanged.
[Verification receipt](verification_receipt.json).

Root review replaced a prototype that repeatedly reconstructed each candidate
prefix, omitted terminal residual information, and had tests that permitted
either success or failure. It also corrected a Bernstein-subdivision proof:
the interval-maximum arithmetic did not give the claimed improvement; the
integrated Bernstein envelope does. Primary-source review corrected an author
attribution and the dimension-dependent cost of the collocation solve.

These corrections, tests and negative results are retained as evidence. They
are not counted as learned research operations or physical validation.

## Next decision

The present schedules do not justify autonomous acquisition. The next
mathematical preflight investigates a different representation: evolve an
explicitly solvable reference generator exactly and certify corrections due
to the remaining interaction. It must pay for representation construction,
mode growth, coefficient/readout enclosures and independent checking, preserve
the same physical calculation, and confront conventional exponential-integrator
baselines. Simple prefix reuse by itself does not close the measured gap.
The larger m1→m2 causal experiment and two untouched incremental gains remain
required before declaring compounding capability.
