# V7 development headroom protocol

Status: prospective development diagnostic, before autonomous method acquisition.
No evaluation-distribution outcomes have been generated. This is not a
compounding experiment or a novelty claim for any supplied algorithm.

## Fixed scientific task

For a supplied rational spin Hamiltonian, uniform local depolarization
`gamma in {0, 1/5, 2}`, and one central Z observable, construct an explicit
piecewise polynomial at horizon `T in {1/5, 1/2}`. The required all-state
observable error is at most `1/1000`, conditional on the declared generator.
Damping and tolerances are identical in every arm and cannot be increased to
rescue a candidate. Two deterministic coupling families (XXZ chain and mixed
XYZ chain with transverse fields) and widths 3, 4, 6 form development problems.
Independent full matrix calculations are restricted to small systems. Sparse
support is capped at 512, polynomial order at 24, basis size at 512, and a
single attempted calculation has a 20-second soft deadline checked between
bounded operations. Refusals and prior unsuccessful attempts are charged.

## Conventional calculation arms

1. Sparse Taylor, complete generated support.
2. BFS coordinate projection with increasing depth.
3. Residual-coefficient-driven sparse polynomial expansion with increasing
   batch size (this is a modest conventional heuristic, not a new method).
4. Two-pass Arnoldi with a projected exponential polynomial.

All arms may use l1, lexicographic anticommuting first-fit, and weight-ordered
anticommuting first-fit. All may use exact monomial or Bernstein residual
integration. These identities and numerical methods are supplied knowledge.
A better implementation or normal form is incorporated in the baseline;
rediscovering it does not pass the autonomous acquisition gate.

For each arm, run the prospectively fixed sequence of settings to the first
accepted certificate. Report accumulated elapsed construction, solution,
witness construction and independent checking time; structural operation
counts; candidate and witness sizes; failures; and final bound. Timing includes
all attempted settings and group/normal-form searches. Counters are not a
complete rational bit-complexity model and are not converted into arbitrary
weighted cost units. Reference calculations and test/validation costs are
recorded separately from candidate execution.

Every accepted certificate binds the external requested horizon. The checker
recomputes every residual from the submitted exact rational polynomial and
supplied exact rational generator. It verifies norm partition coverage,
anticommutation, outward square-root bounds, initial discrepancy and jumps.
The candidate is the exported rational polynomial; there is no unverified
promise about an ideal unrounded numerical algorithm.

## Headroom decision

A prospective new mechanism must improve complete feasible calculation cost
over the strongest conventional arms with identical requirements, not merely
reduce retained observables or tighten a norm bound. Initial development
measurements locate potential headroom only. Any selected mechanism must be
replayed with randomized timing order, repeats and common cases before a
positive gate. A failed run counts as failure; success-only average cost is
insufficient. If conventional methods remove the putative advantage, retain
those methods and close that candidate, while keeping the larger goal active.

## Reserved evaluation distribution

Reserve widths 5 and 7, disordered coefficients, a ladder coupling family, and
central X/Y observables, using fresh seeds derived from the literal namespace
`Spectra-v7-method-transfer-heldout-2026-09-10`. Do not generate outcomes until
m1 and m2, acquisition rules and causal controls are frozen. If this division
proves inappropriate, invalidate and replace it before evaluation, not after
seeing results. The eventual sample size, paired effect thresholds, uncertainty
analysis and net-benefit deployment volume still require freezing before
acquisition. This development protocol alone does not authorize method search.
