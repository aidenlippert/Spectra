# Constructive response investigation — result

**The requested field-level many-body breakthrough has not been established.**
This pass produced a sharp conditional energy bound, exact finite checks, and
a useful representation experiment. It did not prove an efficient construction
for large strongly correlated systems or resolve a new physical phase question.

Five agents investigated response mathematics, prior art, Hubbard coercivity,
counterexamples, and existing implementation seams. Their proposals were reviewed
and corrected before acceptance. The consolidated mathematics is in
[DERIVATION.md](DERIVATION.md); exploratory notes are background, not a substitute
for that derivation or the accepting checker.

## Mathematical result

For shifted blocks H-e=[[A,B*],[B,D]], a trial response X, and

\[
R=B-DX,\quad M=I+X^*X,\quad K=A-B^*X-X^*B+X^*DX,
\]

the certified premises

\[
D\succeq\delta I,\quad K\succeq0,\quad R^*R\preceq\rho^2M,
\quad 2\rho<\delta
\]

imply

\[
\boxed{E_0\ge e-\rho^2/(\delta-2\rho).}
\]

The constant is optimal given those premises: an exact two-dimensional family
defeats every smaller proposed allowance. If delta<=2rho, those facts alone
give no finite uniform lower bound. This is a precise limit on this collection
of scalar information, not on the physical Hamiltonian or all algorithms.

The residual is measured against the lifted state's norm. A large response
therefore need not force a large energy allowance. However, K positivity and
the residual-metric inequality remain full operator obligations. Making their
construction cheap is the unsolved many-body step.

The proof uses established Schur/congruence/residual techniques. Prior work
already supplies rigorous matrix resolvent enclosures, including
[Zimmerling, Druskin and Simoncini (2025)](https://doi.org/10.1007/s10915-025-02799-z).
The exact formula and optimality were derived here; their literature priority
has not been established. They are not being presented as a world-level discovery.

## Actual exact acceptance

A cold four-site square Hubbard calculation at U/t=8, N_up=N_down=2 accepts

\[
E_0/t\in[-1.320235241487274,-1.320234922543856],
\]

with exact rational width **3.1894341811e-7 t**. The independent checker runs
with the Python standard library and verifies every accepted inequality.

Discovery uses 16 block-Jacobi updates and two six-dimensional Ritz solves,
without a full-system eigenvalue/eigenvector solve. It **does enumerate all 36
balanced determinants**, builds 30-by-6 response columns, and checks retained
6-by-6 Gram matrices. This is a correctness control, not a cheap many-body
discovery claim. Its local eliminated-gap proof uses only a 16-state bond PSD
check and correct global occupation counting.

The relative residual allowance is about 12.75% smaller than the independently
verified absolute-residual allowance on this case. That comparison is a small
component result, not a competitive solver benchmark. Units are Hubbard t,
not molecular hartrees, and the four-site square is not a chemical application.

Exact receipt:
`results/constructive_response_20260916/metric_square/replay.json`.
The final fresh component replay and corruption checks are in
`results/constructive_response_20260916/final_exact_replay/`.

## Fully coupled representation experiments

The frozen local rotation was selected from the isolated t=1,U=8 bond, with an
exactly rational orthogonal gate. It was never fitted to the global ground
state. Numerical full-sector states and gaps were subsequently computed solely
for charged diagnostics. Continuations used lambda=0,1/4,1/2,1; lambda=1 is the
fully interacting open ladder. Fragment charges were not fixed.

At full coupling:

| Sites | Enumerated configurations | Bare retained weight | Rung-dressed weight | Stored bare matrix entries | Stored rung-dressed entries |
|---:|---:|---:|---:|---:|---:|
| 4 | 36 | 87.17% | 93.45% | 222 | 850 |
| 6 | 400 | 78.82% | 87.24% | 3,740 | 29,536 |
| 8 | 4,900 | 71.82% | 81.85% | 60,830 | 909,406 |

These are stored floating sparse entries, including any nonzero arithmetic
residue; they are not lower bounds on exact operator representation size.
Rung dressing improved overlap but did not establish an overall cost advantage.

Applying local rotations to every bond raised the six-site retained weight to
**98.93%**, while the expanded matrix stored 135,248 entries. The eight-site
expanded approach exceeded its declared two-million-entry cap after six gates.
That stopped run is retained as a limitation, not counted as a completed result.

The investigation then removed that avoidable matrix expansion: apply the
circuit, the original Hamiltonian, and the inverse circuit directly. This
completed all ten gates for the eight-site case, reaching **98.33% retained
weight**, with zero stored transformed-Hamiltonian entries. It stored 60,830
original-H entries, 97,000 gate entries, and 338,100 response entries. The process
took about 1.84 seconds and peaked at 164 MB on this machine. This timing is a
numerical diagnostic, not time to a complete certified ground-energy interval.

Crucially, the action still uses 4,900-dimensional enumerated vectors. It fixes
the expanded-matrix implementation failure; it does not remove determinant
enumeration, certify the eight-site residual, or settle many-body scaling.
The retained spin spaces also still grow: 6, 20, 70 states in these examples.

All twelve bare/rung continuation cases and the subsequent circuit experiments
are in `results/constructive_response_20260916/`. Half-filled small ladders are
correctness diagnostics, not the difficult doped two-dimensional physics
benchmark needed for a field-level claim.

## Exact structural limitation

For L disconnected U/t=8 dimers, a global bare no-doublon projector has ground
weight p=((1+2/sqrt(5))/2)^L. The full physical gap stays positive, while the
eliminated gap above the exact ground energy lies between

\[
(-e_d)p\quad\text{and}\quad(-Le_d)p/(1-p),\qquad e_d=4-2\sqrt5.
\]

The response norm is at least sqrt((1-p)/p). Exact rational replay evaluates
this family through 256 dimers without generating the global state space.
This rules out assuming uniform conditioning from a local charge cost alone.
It does **not** prove exponential computation at fixed energy tolerance or
obstruct dressed/local projectors. It is a disconnected analytic control,
not the fully interacting result the campaign ultimately needs.

## Verification and accounting

- 29 focused exact/numerical tests cover fermion signs, independent model
  construction, rational PSD checks, two successive Schur eliminations,
  normalized dressing, sharpness, invalid-premise refusal, and circuit action.
- The standard-library replay checks 24 exact all-occupation bond inequalities,
  dimer radical enclosures, response identities, composition, the four-site
  certificate, and exact sharpness counterexamples. Four corrupted certificate
  inputs are rejected.
- Earlier candidate errors in energy normalization, symmetrization, occupation
  counting, and operator-order statements were caught and corrected. Agent
  assertions and early tests were not treated as acceptance.
- Run receipts record elapsed subprocess time, peak child RSS and outcomes.
  `ACCOUNTING.json` includes failed/stopped work distinctions and exclusions.
  Reasoning, code authoring, source research, and early unmetered agent checks
  are not included in subprocess sums. No fresh end-to-end scaling claim follows.
- Existing molecular results were preserved. No paid compute, installation,
  commit, push, or GitHub operation was used.

## What remains unsolved

There is still no directly constructible, compact certificate for K and the
relative residual metric whose cost is shown to remain controlled through
strongly interacting reductions. Nor is there a new controlled low-temperature
2D phase result, broad many-body solution, or demonstrated advantage over the
field's strongest methods.

The exact lemma and its sharpness narrow the mathematical obligation: improving
the scalar constant cannot supply the missing scalability. A successful further
result needs additional correlation structure that makes the operator proofs
constructible and composable. This pass has not supplied that structure.
