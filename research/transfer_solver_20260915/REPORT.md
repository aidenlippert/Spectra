# Cold construction, transfer, and measured controls

The prior **0.767448135463 mHa H8 milestone remains achieved and preserved**.
The new input-derived constructor produced a fresh H8 certificate of
**0.605110981693 mHa**, including fresh integral generation, MPS discovery,
nonsinglet construction, coefficient preparation, optimization and exact replay.
The original H8 rational Hamiltonian and the new one are exactly unitarily
equivalent under two paired orbital sign flips. This is a complete cold run,
not a repeat conditional on an old checkpoint.

## Frozen cases and complete cost

The changed-input cases were frozen together before their outcomes were known.
H4 was the development control. A frozen procedure adapts its algebra to the
input; it does not reuse H8 dimensions or parity labels.

| Case | Certified total width (mHa) | <=1.6 mHa | Complete campaign seconds | Peak process MB | Last step |
|---|---:|---|---:|---:|---|
| h4_control | 0.001544748 | yes | 27.84 | 328.6 | replay |
| h8_cold | 0.605110982 | yes | 858.83 | 908.3 | replay |
| h6_asymmetric | 0.021516666 | yes | 411.12 | 388.6 | replay |
| water_asymmetric | 0.025861832 | yes | 540.93 | 474.5 | replay |
| h10_size | No complete interval | not established | 902.51 | 1572.5 | state |
| h10_bond64 | No complete interval | not established | 603.09 | 1082.5 | state |

If present, `h10_bond64` is a separately declared exploratory retry after the
original H10 state/upper budget failure. It uses a fresh random state capped at
bond 64 and six sweeps, plus an existing faster exact proposal-side spin twirl.
The complete lower family, optimizer and accepting checker remain the same.
This retry is not counted as a success of the original untouched frozen rule;
the original failed cost remains additional.

All completed certificates cover the entire fixed-N sector, using separate S=0
and M_S=1 proofs plus one exact spin-defect allowance. Their upper Rayleigh
quotients and lower remainders use integer/rational arithmetic, without numerical
packages on the accepting path. The asymmetric H6 case has only the trivial
total-parity class; water changes the molecule and electron/orbital counts.
These are a few deliberately selected cases, not a statistical reliability test.

The per-case times start with new integrals. Historical state, coefficient maps,
optimization checkpoints, and nonsinglet certificates are absent from their
constructor inputs. The fresh MPS supplies moments to discovery. Numerical FCI
references below were run only after discovery and never supplied as teachers.

The H8 search retains **343,424 singlet Gram entries**, with largest block 260,
versus 335,168 and 256 in the earlier winner. The declared difference is that
mixed linear words are retained uniformly at every size. This small enlargement
means the fresh run is not an identical-cone reproduction of the previous winner.
Fewer exported factors are not used as a proxy for discovery cost.

## Equivalent and restricted representation controls

| Control | Exact width (mHa) | Singlet Gram entries | Control prep/solve/replay seconds | Same stages in full highest-weight case |
|---|---:|---:|---:|---:|
| h4_control_magnetic | 0.000363332 | 19248 | 56.46 | 14.07 |
| h8_cold_magnetic | 94.301434775 | 1230272 | 794.92 | 548.01 |
| h8_restricted64 | unavailable | 55616 | 384.51 | 548.01 |

The magnetic controls use the same attainable spin-averaged cone, physical zero
initial Gram matrices, new MPS functional, optimizer and stopping rule. Shared
integrals, state, upper and nonsinglet costs are excluded from both sides of
these representation-only comparisons. The restricted64 control instead retains
at most the first 64 deterministic rational columns per mixed highest-weight
block. It tests a declared smaller family; it is not the best possible restriction.
No failed finite run is promoted to an exact family obstruction.

The mathematical cone equivalence, exact kernel checks and acceptance rule are
in [DERIVATION.md](/Users/aidenlippert/Documents/Spectra/research/transfer_solver_20260915/DERIVATION.md).
Symmetry reduction itself is established; see
[Gatermann and Parrilo](https://arxiv.org/abs/math/0211450).

## Comparisons at the same target and with numerical FCI

The enumerated baseline constructs its own integer upper and rational
Cholesky-remainder lower. It covers all fixed-N determinants, partitioned only
by alpha count, and independently replays the resulting certificate. It is
capped at 2,000 determinant labels. Its time includes the common integral
generation cost once for each method comparison.

| Case | Enumerated exact width (mHa) | Generation + construction + replay seconds | Enumerated labels |
|---|---:|---:|---:|
| h4_control | 0.800000036 | 1.980 | 70 |
| h6_asymmetric | 0.799999982 | 5.842 | 924 |
| water_asymmetric | 0.799999957 | 5.493 | 1001 |

Both methods target <=1.6 mHa; the table states their actual attained widths.
This is a small rigorous reference implementation, not a survey of the strongest
available certified solvers. It can directly reveal when eliminating enumeration
has not yet paid off. No overall competitive advantage follows from this campaign.

| Case | Numerical FCI solve seconds | Converged | Spectra upper minus FCI (mHa) | FCI minus Spectra lower (mHa) |
|---|---:|---|---:|---:|
| h4_control | 0.0032 | True | -0.000000 | 0.001545 |
| h8_cold | 0.0396 | True | 0.423752 | 0.181359 |
| h6_asymmetric | 0.0214 | True | 0.000000 | 0.021517 |
| water_asymmetric | 0.0232 | True | 0.000000 | 0.025862 |
| h10_size | 0.3110 | True | nan | nan |
| h10_bond64 | 0.2849 | True | nan | nan |

FCI uses fixed tight numerical convergence, not a rigorous lower guarantee.
These energies are for the floating integrals before rational rounding; the
recorded rounding envelope and numerical comparison slack are explicit in the
receipts. Its timings and numerical agreement must not be presented as exact
error bars. See the [PySCF FCI documentation](https://pyscf.org/user/ci.html).

## Independent checks and fragment composition

`independent_h4_oracle.json` records a separately implemented full bit-state
Hamiltonian, rational LDL positivity test and explicit integer MPS expansion.
It imports neither the original CAR nor tensor checker. Its enumeration is
charged as a validation control, not hidden in construction. The raised-lower
mutation must fail. This strengthens validation on a small system; it is not
a formal verification of every implementation path.

`fragments/summary.json` records exact composition of 2, 4 and 8 independent
H4 fragments with fixed local electron counts, including direct exact product-MPS
uppers. Removing those local constraints or adding cross-fragment hopping is
refused. This establishes the compositional rule, not size consistency of a
monolithic restricted optimizer. The latter remains a separate question.

## CH2 physical-model control

The existing exact fixed-geometry STO-3G CAS(6,6) gap near 25.49 kcal/mol remains
a certificate for that model, not a reliable experimental prediction. The new
small control optimizes the A1 singlet and B1 triplet CASSCF(6,6) geometries,
checks total spin at every evaluation, and adds SC-NEVPT2 single-point corrections.

| Basis | CASSCF E_S-E_T (kcal/mol) | SC-NEVPT2 E_S-E_T (kcal/mol) |
|---|---:|---:|
| cc-pvdz | 11.93775 | 13.10483 |
| cc-pvtz | 10.51181 | 12.30047 |

The experimental study reports 9.00 +/- 0.09 kcal/mol for the observed splitting
and an inferred electronic separation of 8.7 +/- 0.5 kcal/mol. The latter is the
appropriate rough electronic comparison here; no zero-point term was computed.
[Leopold, Murray, Miller and Lineberger](https://experts.umn.edu/en/publications/methylene-a-study-of-the-xsup3supbsub1sub-and-%C3%A3sup1supasub1sub-st/).
These new conventional calculations have no certified physical error bar.
Basis sensitivity is not an error enclosure, and agreement is retrospective.
The [CASSCF](https://pyscf.org/user/mcscf.html) and
[SC-NEVPT2](https://pyscf.org/user/mrpt.html) components are existing methods,
so the control does not establish a Spectra-specific predictive advantage.

## All attempts and preservation

The ledger contains **102 instrumented attempts**, summing to
**6875.312 process-wall seconds** and
**6472.750 CPU seconds**. Maximum single-process
RSS was **1798.111 MB**. One heavy process
ran at a time. This includes unsuccessful attempts and the controls; it is not
the cost of one fresh solution. Editing, browsing, reporting and preservation I/O
are outside this compute ledger, and earlier R&D costs stay in their earlier
ledgers rather than being silently zeroed.

- algebra_initial: failed, 2.734 s; log in `runs/algebra_initial.log`.
- h10_bond64_state: timeout, 600.046 s; log in `runs/h10_bond64_state.log`.
- h10_size_state: timeout, 900.086 s; log in `runs/h10_size_state.log`.
- h4_control_fci: failed, 0.900 s; log in `runs/h4_control_fci.log`.
- h4_control_fci_fixed: failed, 0.844 s; log in `runs/h4_control_fci_fixed.log`.
- h6_asymmetric_fci: failed, 0.524 s; log in `runs/h6_asymmetric_fci.log`.
- h6_asymmetric_fci_fixed: failed, 0.516 s; log in `runs/h6_asymmetric_fci_fixed.log`.
- h8_cold_fci: failed, 0.540 s; log in `runs/h8_cold_fci.log`.
- h8_cold_fci_fixed: failed, 0.504 s; log in `runs/h8_cold_fci_fixed.log`.
- h8_restricted64_stage1: failed, 140.495 s; log in `runs/h8_restricted64_stage1.log`.
- water_asymmetric_fci: failed, 0.548 s; log in `runs/water_asymmetric_fci.log`.
- water_asymmetric_fci_fixed: failed, 0.516 s; log in `runs/water_asymmetric_fci_fixed.log`.

All **12,497 previously protected files** and all frozen constructor
source hashes match. No older achieved result has been overwritten.

Cold success and limited transfer are now evidenced. Broad reliability,
manageable larger-system scaling, the attainable-family ceiling, monolithic
fragment size consistency, competitive advantage over strong methods, and
prospective experimental value require their own evidence. A numerical miss
must be diagnosed using its upper, lower, residual and budget receipts rather
than being called a proof of impossibility.

- [Complete compute ledger](/Users/aidenlippert/Documents/Spectra/results/transfer_solver_20260915/cost_ledger.json)
- [Frozen validation protocol](/Users/aidenlippert/Documents/Spectra/results/transfer_solver_20260915/validation_set_protocol.json)
- [Exact original-H8 bridge](/Users/aidenlippert/Documents/Spectra/results/transfer_solver_20260915/h8_original_model_bridge.json)
