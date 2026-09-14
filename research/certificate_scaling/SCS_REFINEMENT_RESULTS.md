# Latest verified milestone: H4–H10 accuracy ladder passes

Independent python -S replay confirms H10 width **.001157618253846804 Ha**,
below .0016. H4/H6/H8 widths are 1.53002e-7,5.49944e-5,7.50585e-4 Ha.
Machine scoreboard: `results/certificate_scaling/active_space_ladder_latest.json`.
The goal remains ACTIVE: no general structural accuracy theorem or scalable
compact-discovery condition for chemistry has been established.

H10 used the full chosen mixed-cubic dictionary: 4674900 Gram entries,
4810 exported factors and 2581696 nonzero factor coefficients. Its source
run took2222.95s; spectral residual proof construction/replay took255.11s;
independent local interval replay took182.88s. The original FCI upper-witness
discovery cost also remains. This is a finite frozen-model accuracy milestone,
not a competitive comparison with other research teams or a scaling proof.

Artifacts: downloaded/h10_evd_1800 under the Lambda SCS campaign, spectral
proof `scs_refinement/h10_evd_1800_spectral`, interval `intervals/h10_evd_1800.json`.
Certificate SHA397eb5bbed0dc09b24169d08307f448aea1d47749828d7d32f1470a912cdf7c0.
All20 new output downloads hash-match their remote originals. Both hosts are
now idle and intentionally warm; the H10 and H6 jobs below are historical.

The basis-conditioning experiment is documented in polynomial_basis_conditioning.md.
Its H6 Clarabel/SCS intervals .0560816/.210587 both miss accuracy; the H4
integration passes at1.23319e-6 with3.235s total time. Exact spin algebra and
nearby-H symmetry checks now establish a controlled starting point for the
next experiment; see spin_irrep_frontier.md. A spin-adapted discovery backend
has not yet been implemented.

--- Prior campaign details, superseded above where statuses differ ---

# SCS refinement and an exact commutator obstruction

The structural-scaling goal remains open. This campaign produced an exact
negative result for a specific Hamiltonian-derived dictionary, explicit SCS
iterate checkpoints, and a measured PSD projection improvement. None proves
general scalable discovery or solves the downstream chemistry objectives.

## An exact mathematical elimination of one candidate

The quadratic baseline was enriched with cubic parts of `[H,a_i]` and their
adjoints. A larger version splits each cubic commutator by creation-orbital
index. All generators come from H; no ground-state vector, source certificate
or physical upper bound enters dictionary discovery. The Gram maps retain
every coefficient through degree six and include a body-two number ideal.

The six completed numerical experiments are documented in
`commutator_dictionary.md`. On H4 the baseline and summed commutators gave
interval widths about .00448 Ha; creator channels improved that to .003564 Ha.
On H6 the corresponding widths were .013006, .013007 and .012208 Ha. All miss
.0016 Ha. The H6 creator dictionary has 13.6 times fewer Gram entries than
the full cubic dictionary but 10.9 times more coefficient-map nonzeros:
smaller Gram matrices alone do not establish cheaper discovery.

We then constructed an **exact dual witness** for the enlarged H4 creator
dictionary. The independent standard-library checker rebuilds all generators
directly from H, uses full charge blocks without symmetry pruning, checks all
443 body-two ideal equalities, verifies `y(1)=1` and `|y_w|<=1`, and performs
exact rational PSD elimination on blocks of dimensions 64,64,28,28,64.
Missing canonical moments are explicitly zero. A separate exact generator
comparison confirms every numerical generator is a scalar multiple of a
rebuilt generator, so the checked cone contains the numerical cone.

For `H = b I + SOS + ideal + r`, positivity gives
`b - ||r||_1 <= y(H)`. The accepted witness gives

```
y(H) = -51387812558113145851462206230521 / 14000000000000000000000000000000
     = -3.6705580398652247...
```

The checker independently replays the full cubic H4 primal certificate,
with the identical rational Hamiltonian, obtaining a ground-energy lower
bound of `-3.667000108966591573`. Therefore every certificate in the restricted
L1-penalized cone undershoots the true ground energy by at least
**.00355793089863313 Ha**, exceeding .0016 Ha. This is a proof of failure for
this candidate on this H4 model, independent of solver convergence. It also
excludes its quadratic and summed-commutator subfamilies. It does not exclude
other dictionaries, higher ideal degree, or different residual-bound methods.

The rational witness is 85,170 bytes with maximum moment bit size 64. Proposal
and exact-tested repair took 7.28 seconds; separate `python -S` replay took
about two seconds. Mixtures 0, 1e-8 and 1e-7 were refused by exact PSD checks;
mixing 1e-6 of the uniform fixed-N moment functional passed. No floating
eigenvalue is used by the accepting checker.

Artifacts:
`results/certificate_scaling/commutator_dictionary/h4_exact_dual/`
(`witness.json`, `receipt.json`, `independent_replay.json`, `span_inclusion.json`).
Implementation: `commutator_dual_witness.py`. The off-diagonal dual entry is
`y(p_i† p_j)`; using the symmetrized upper-triangle primal coefficient without
dividing by two would be incorrect. A focused test catches this convention.

## Explicit restart is available, but its numerical quality must be measured

Both installed CVXPY versions cache an SCS solution automatically only for
`optimal`, not `optimal_inaccurate`. The new `scs_checkpoint.py` saves canonical
SCS x/y/s after finite inaccurate solutions and verifies exact canonical data,
cone, version and solver-option identities before loading. Budget changes are
allowed; changed problems, other options or corrupt array files are refused.
This restarts iterates while rebuilding SCS scaling, acceleration history and
factorization; it is not continuation of the full internal solver workspace.

The H4 fresh-process restart integration passed exact replay at interval width
3.13026565e-6 Ha. Twenty focused tests cover checkpoints, existing adaptive
discovery, commutator algebra and dual acceptance/refusal paths.

The matched H10 runs both finished on the two authorized warm A10 hosts:

| Protocol | SCS iterations | Total solve wall s | Total discovery/export s | Coefficient-L1 lower bound |
|---|---:|---:|---:|---:|
| One nominal 1200 s solve | 325 | 1216.414 | 1638.083 | -12.406926470236069 |
| Two nominal 600 s stages | 175 + 175 | 1300.941 | 1734.723 | -13.731550073449739 |

Actual time overruns and restart overhead are charged. The restarted raw
primal/dual gap was smaller, but its exact coefficient residual was much worse
(1.40226 Ha versus .0367167 Ha). A smaller floating solver gap is not evidence
of a better certificate. Independent spectral strengthening and full
`python -S` two-sided replays completed for both. The continuous run improved
the best H10 interval width from .02334631 to **.0078930925603 Ha**. The
restarted run gave .2082457959935 Ha. Both still fail .0016 Ha. Spectral witness
construction/replay added 262.81/279.64 seconds remotely respectively; the
separate local end-to-end replay times are in `scs_refinement/intervals/`.

All source snapshots, raw canonical checkpoints and certificates are retained
under `results/lambda_runs/scs_refinement/`. Neither run used the reference
upper as discovery input. The H10 physical upper remains an expensive FCI
validation witness, not a scalable discovery result.

## Measured PSD projection opportunity

SCS attributed 95.3% of the first H10 stage's solve time to cone projections,
3.5% to the linear system. The installed direct SCS extension links bundled
OpenBLAS 0.3.15 (`Prescott SINGLE_THREADED`). The
[SCS 3.2.8 projection source](https://raw.githubusercontent.com/cvxgrp/scs/3.2.8/src/cones.c)
calls LAPACK `syev` for its PSD eigenstep.

Using that exact bundled library on eight actual H10 blocks, with one thread,
two warmups and three timed trials:

| Batch | dsyev median s | dsyevd median s | Ratio |
|---|---:|---:|---:|
| 8 x 725 | 2.274826 | 1.150596 | 1.98 |
| 8 x 225 | .090833 | .048784 | 1.86 |

Reconstructed matrices and PSD projections agreed within 2.43e-17. These
measurements exclude copies, workspace allocation, projection reconstruction
and full solver overhead. They justify an isolated solver/backend experiment;
they do not establish an end-to-end speedup or a structural scaling theorem.
Executable harness: `projection_kernel_gate.py`; source/library hashes and
all trial timings: `results/certificate_scaling/scs_refinement/projection_gate/`.

A controlled SciPy backend comparison (OpenBLAS 0.3.27.dev, SkylakeX, runtime
thread count verified as one) gave eight-block `evd` medians .448795 s at 725
and .025816 s at 225. Its corresponding `ev` times were 1.759284/.071982 s.
Cross-driver PSD agreement was within 3.21e-17. All batch trial timings and
library hashes are retained in `alternate_projection_gate/`. The earlier
agent benchmark is explicitly excluded: it changed thread settings after
imports and confused per-block medians with batch timings.

An isolated Linux LP64 preload adapter now calls the faster eigenbackend
without replacing installed libraries. Root corrected workspace allocation
and a separate OpenBLAS symbol collision exposed by gdb. Stock/preload small
SDP controls matched within 3e-14, with 504 intercepted calls and no failures.
An H4 molecular integration passed independent exact replay at width
3.89191924e-7 Ha. See `scs_evd_preload.md` for failure history, build command,
scope and receipts. A new cold H10 600-second run is active under the isolated
source snapshot `evd_source.tar.gz`; its result is not yet part of this table.

The model scope throughout is the same frozen rational electronic Hamiltonian
for straight H chains at 1.4 Angstrom in STO-3G/RHF canonical orbitals. Claims
do not cover basis error, reaction barriers, finite temperature or experiment.

## Latest completed H10 acceleration result

The 600-second accelerated PSD-backend run completed 825 SCS iterations. Its
independently replayed spectral-residual interval is **.00386819649650 Ha**,
improving .00789309256030 from the stock continuous1200 run, while still
missing .0016. The original coefficient-l1 interval is much weaker; the exact
spectral residual proof is part of the accepted artifact and its cost counts.
This is a convergence/backend result using full cubic discovery, not a compact
structural discovery theorem.

The certificate lower bound after spectral replay is -12.372167652819943;
the independently replayed reference upper is -12.368299456323442. The source
certificate has 4907 factors and 2656820 nonzeros. The full discovery retains
4674900 Gram scalar entries and 2347238 map nonzeros. Its recorded total time
is 1035.50 seconds, with 155.49 seconds map construction, 626.43 seconds solve
wrapper and 253.31 seconds exact export. Spectral proof construction and a
separate standard-library replay are additional costs. The reference upper's
FCI discovery is still an exponential validation cost.

Artifacts: `results/lambda_runs/scs_refinement/downloaded/h10_evd_600/`,
`results/certificate_scaling/scs_refinement/h10_evd_600_spectral/`, and
`results/certificate_scaling/scs_refinement/intervals/h10_evd_600.json`.
Certificate SHA256:
`e73a14d4cad8fe48f03a0737c697c9f5533fd44cbfd97dad7f8233c8a8a7b1aa`.
Twenty-three newly downloaded output files across this run and the H6 runs
matched their remote originals; see `new_download_hashes.json` in the Lambda
campaign directory. Logs were retained separately.

A fresh cold H10 solve with the same accelerated backend and an 1800-second
solver budget is still running on host B. It is not a restart of the saved
600-second iterate. Its result has not been counted. Host A's H6 experiments
are terminal, and both hosts remain warm by user authorization.

## A separate structural control now works

`dressed_quadratic_structure.md` derives a restricted matching-CZ recognizer
and a compact certificate compiler. Unknown matching recovery uses graph XOR
constraints and exact full-H regeneration. It produces both rational SOS and
rational dressed-Slater witnesses without FCI. Separate python -S replays pass
for 4,8,16,32,64 modes, with M factors and O(M^2) stored coefficients. At 64
modes the interval is 7.0554e-6 Ha and the full witness is 293111 bytes.

This is hidden-free physics with strong structural restrictions. The four
molecular fixtures explicitly refuse it. It is a functioning positive control
for structure-driven discovery and verification, not a general chemistry
breakthrough. Fixed-precision numerical discovery has no universal accuracy
or termination guarantee. The main goal remains active.
