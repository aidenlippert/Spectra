# Next basis attack: what residual expansion actually supplies

`residual_basis` is not a trajectory-wide adaptive Galerkin method. At step
`j` it computes one raw `G c_j`, adds at most `batch` labels outside the
current basis according to that step's coefficient magnitude, and projects
the next coefficient onto the enlarged set. A label omitted at step `j` is
represented only by an l1 residual witness; it is not available to generate
later coefficients. The `work` mode in `pruned_taylor` similarly ranks one
step's candidates by magnitude divided by a one-column work estimate. It is a
per-coefficient drop rule, not a basis-rank or accumulated-residual rule.

## Standard exact adaptive Galerkin alternative

For a finite Pauli basis B, form the exact projected matrix
`A_B=P_B G|_B` by applying every generator column in B. Propagate the exact
projected Taylor recurrence `u_0=P_B o`,
`u_(k+1)=A_B u_k/(k+1)`. At every order retain the full omitted residual
`r_k=(I-P_B)G u_k`, including collisions and signs. A residual-driven
expansion scores labels by their accumulated integration contribution
`|r_k(p)| T^(k+1)/(k+1)` (or by a predeclared interval norm), adds a declared
batch, and restarts the projected recurrence. Stop only when the sum of all
retained residual weights plus the Taylor tail has an independently checked
bound. Exporting the ordinary coefficient maps then passes the unchanged
power or Bernstein residual checker; no checker change is needed.

This is a conventional adaptive Galerkin/residual-enrichment method. It is
algebraically different from `residual_basis` because the selection state is
the accumulated trajectory residual and every column of the current basis is
part of one explicit projected operator, rather than one-step coefficient
rankings.

## Finite rank counterexample

Consider the exact linear operator on basis `{e0,e1,e2,e3}`:
`G e0 = δ e1 + e2`, `G e1 = K e3`, and `G e2=G e3=0`, with initial `e0`.
For `0<δ<1`, a one-label magnitude expansion with batch one chooses `e2`
and omits `e1`. Its projected recurrence terminates after the first term and
cannot represent the true second-order coefficient `δ K e3/2`. A rank-two
enrichment that chooses e1 (or a trajectory-residual score that accounts for
future column work) captures the path `e0→e1→e3`. This proves that basis rank,
not merely retained first-step coefficient mass, can control the needed
trajectory. It is a mathematical counterexample only; it is not a workload
instance or evidence of a headroom win.

## Cost and proof obligations

An adaptive Galerkin candidate must charge: every generator column used to
build `A_B`; all exact coefficient multiply-adds and collision merges; sorting
or scoring every residual label; basis insertion and duplicate checks; every
restart and failed batch; rational conversions/bit growth; support and
polynomial storage; witness construction for power or Bernstein integration;
and a fresh independent checker replay. The checker must recompute
`G c_k`, verify every residual and jump, and not trust the projected matrix or
the residual transcript.

The smallest plausible constructive difference is therefore one accumulated
residual score plus exact projected-column reuse, compared against the
strengthened Taylor (26/36 successes) and BFS (27/36 successes) baselines. It could improve
success only when a low-magnitude frontier has high downstream leverage; it
could reduce cost only if avoided columns and fewer restarts exceed scoring,
matrix-build, and checker costs. No such cost-positive measurement was run in
this bounded algebraic preflight, and no threshold or tolerance change is
admissible.

## Decision

The existing residual expansion should be described as a conventional
per-coefficient projection heuristic. A trajectory-wide exact Galerkin method
is well-defined and has a finite rank counterexample, but there is no
defensible total-cost or success advantage yet. Keep this as an untested
baseline candidate; do not call it acquired, and do not open held-out cases.

## Physical-applicability correction

The nilpotent high-gain example above is a generic linear-system example, not an admitted Lindblad generator. It does not justify amplification of a small operator-norm residual under positive unital evolution, which is contractive. Adding e1 alone also cannot represent e3; the intended enrichment adds both e1 and e3 in addition to the initial direction. The example distinguishes heuristic choices but supplies no physical headroom theorem. An actual candidate must exploit accumulated residual structure or representation cost in the admitted Pauli generator and pass that generator's exact checker.
