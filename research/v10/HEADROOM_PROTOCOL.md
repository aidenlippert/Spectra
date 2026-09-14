# V10 development-only solvable-reference test

Use the identical 36 V7/V9 inputs: n=3,4,6; XXZ/mixed; gamma=0,1/5,2;
T=1/5,1/2; central Z; total operator-norm tolerance 1/1000. Do not change
Hamiltonian coefficients, damping, horizon or observables. Reserved n5/7
evaluation remains unopened. This is one diagnostic pass, not a positive
timing-confidence test. A promising result would require independent paired
repetitions before any acquisition gate.

Three supplied conventional arms: V8 fraction-free Taylor, V9 strengthened
BFS, and the V10 onsite-Z reference Duhamel expansion. No learned artifact is
present. All arms export a rational Hermitian Pauli operator at the requested
time. Time includes construction, successful and failed certificate attempts,
fresh checking, and endpoint evaluation. Serialization is outside calculation
time and reported separately through artifact sizes. Refusal is recorded, not
treated as impossibility.

Reference caps: order24, 512 live eigenoperator labels and cached columns,
50000 total layer entries, polynomial degree32 and rational bits8192. Reference
coordinate counts are reported separately from Pauli output support. The full
Pauli output must also fit512. Choose the first order with factorial bound at
most 999/1000000. Allocate 1/1000000 to outward endpoint readout so the common
total error remains <=1/1000. A scalar preflight may refuse when no order <=24
can meet this bound; its calculation counts. Actual numerical error never
substitutes for the certificate.

Checker verifies D0, exact correction differential identities, zero correction
initial conditions, conjugate Hermiticity and the fixed model split. It derives
K=2 sum|V_P| itself. It uses the proved contraction/factorial bound, and the
endpoint module encloses all exponential evaluation and grid errors. The
original V7 checker is unchanged; the separately implemented V10 checker is
not a learner and must be frozen before any future acquisition using it.

Arm order is shuffled with a fixed seed10 within each case. No tuning after
viewing this diagnostic. The comparison measures actual Python implementation
cost, not an optimality claim for exponential integrators as a class. A stronger
known bound or representation would belong to conventional baselines too.
