# Molecular collective-completion results

This pass produced new compact molecular certificates and a complete transfer test. Exact receipts determine the outcomes below; numerical solver status never accepts an energy bound.

## Strongest complete intervals

- **h6: 1.034856 mHa**, 18,228 PSD entries in the retained solve. Certificate: `results/collective_completion_20260914/candidates/h6_actual_spin_r32_dual8/round_0/certificate.json`. Full interval: `results/collective_completion_20260914/candidates/h6_actual_spin_r32_dual8/round_0/interval.json`.
- **h8: 2.874956 mHa**, 54,932 PSD entries in the retained solve. Certificate: `results/collective_completion_20260914/candidates/h8_spin_r32_invariant_dual48/round_0/certificate.json`. Full interval: `results/collective_completion_20260914/candidates/h8_spin_r32_invariant_dual48/round_0/interval.json`.
- **fresh_h6: 1.034056 mHa**, 18,003 PSD entries in the retained solve. Certificate: `results/collective_completion_20260914/candidates/fresh_h6_2p23_spin_r32/round_1/certificate.json`. Full interval: `results/collective_completion_20260914/candidates/fresh_h6_2p23_spin_r32/round_1/interval.json`.

The H8 singlet construction additionally requires its separate nonsinglet proof: 18,128 PSD entries and its own discovery/replay cost. Do not omit this component when comparing total proof size.

## What changed

The paired cubic identities cancel sextic terms before moment preparation. The constructor now requests only 499 H6 or 1,525 H8 quartic moment coordinates. The H8 guiding-moment step took 0.890 s in the first direct run, compared with 24.145 s in the preceding full preparation; its whole direct preparation took 23.417 s. It still visited 1,206,592 cubic word pairs. This is a polynomial operator-dictionary dependency, not a full determinant enumeration and not a scalability theorem.

The search retains all cross terms in the chosen operator spans and uses the coefficient dual to find omitted positivity directions. It also tested independent linear responses and exchanging numerical subspaces. Subspace selection remains heuristic; exact reconstruction pays every actual error.

Exact balanced-spin lifting was implemented with the original Hamiltonian spin defect charged. On H8 a separate bound covers every nonsinglet of the spin-invariant Hamiltonian. Singlet constraints are therefore used inside a complete two-part proof. This is not a claim about the first excitation gap of the entire molecule.

The compact singlet proof can itself be spin averaged before reconstruction. An exact integer projector and modular row-rank certificate reduce 1,525 equations to 581. The first such run solved in 63.844 s; the preceding formulation took 257.895 s. The formulation and thread setting changed, so this is not an isolated one-variable timing experiment. See `MATHEMATICS.md` for the full soundness argument and primary-source context.

A cold original-H6 run then used the direct quartic preparation and the equivalent moment-dual SDP, without any preceding search checkpoint. Its two-round process took 19.632 s and certified 1.048913 mHa (1.048971 mHa after rechecking and outward-rounding the actual upper endpoint). Its 18,228 retained PSD entries are 11.92 times fewer than the old 217,268-entry construction. Direct preparation and trial-state discovery/verification are additional costs; this is not an end-to-end speed ratio.

Larger H8 attempts included global pure triples and sharing factors related by spin flip. The first two larger primal formulations hit their 480-second process limits without certificates. A one-second profile of the moment formulation found its active stack inside the optimizer's numerical LDL factorization. Smaller file or variable counts alone did not resolve that cost. Those timed-out attempts provide no new lower endpoint.

A response test that retained the old paired span reached its time limit with an invalidly indefinite numerical Gram proposal; exact repair gave a wide interval. The same cone in the SCS moment formulation completed 8,800 iterations, but also failed to improve the retained exact bound. Numerical status labels therefore remain separate from certificate acceptance.

The coherent-square variant avoids growing cubic PSD blocks. Four seeded H8 solves took 136.658 s as a process. The final singlet solve used 18,976 Gram entries and 509,898 projected-map nonzeros, plus the separate 18,128-entry nonsinglet proof. Its exactly verified interval was 2.880545 mHa. This is close to the best 2.874956 mHa interval, with a smaller retained optimization; the seed's earlier discovery remains a dependency. A separate cold trial uses only MPS quartic moments and records every round below. The two finite search cones differ, so a smaller representation is not evidence of equal expressive power.

## Independent transfer

The rule was frozen before generating a new H6 chain at 2.23 Angstrom, and no prior fixture at that spacing was found. Source copies are in `results/collective_completion_20260914/transfer/source_at_freeze`. The generator disabled FCI. A random-initialized bond-48 MPS, direct two-particle preparation, rank-32 lower search, and one dual refinement were constructed on this new input. The two exact intervals were 1.421065 and 1.034056 mHa. The two lower solves took 14.401 and 23.613 s internally. Fixture generation, MPS discovery, preparation, acceptance, and the failed output-directory launch are additional recorded costs.

All recorded fresh-case processes sum to 90.375 seconds, including fixture generation, trial-state discovery, direct preparation, both lower solves, failed launches, and additional replay checks. This includes the entire measured fresh-case computation; it excludes prior method development and environment setup. The integrated bundle independently rechecks the actual MPS upper and the new lower proof.

## Candidate ledger

| Candidate | Exact width (mHa) | PSD entries | Search + construction (s) |
|---|---:|---:|---:|
| fresh_h6_2p23_spin_r32/round_0 | 1.421065 | 14,244 | 15.056 |
| fresh_h6_2p23_spin_r32/round_1 | 1.034056 | 18,003 | 39.156 |
| h6_actual_spin_r32/round_0 | 1.718501 | 14,244 | 168.663 |
| h6_actual_spin_r32_cold/round_0 | 1.728479 | 14,244 | 6.736 |
| h6_actual_spin_r32_cold/round_1 | 1.048913 | 18,228 | 18.710 |
| h6_actual_spin_r32_cold_wedge/round_1 | 1.039047 | 18,228 | 27.094 |
| h6_actual_spin_r32_dual8/round_0 | 1.034856 | 18,228 | 30.200 |
| h6_actual_spin_r32_response/round_0 | 1.969413 | 15,924 | 18.545 |
| h8_spin_r1_nonsinglet/round_0 | Sector proof only | 18,128 | 8.324 |
| h8_spin_r32_cold_atoms/round_0 | No complete replay | 18,608 | 33.866 |
| h8_spin_r32_cold_atoms/round_1 | No complete replay | 18,820 | 65.375 |
| h8_spin_r32_cold_atoms/round_2 | No complete replay | 19,026 | 97.144 |
| h8_spin_r32_cold_atoms/round_3 | No complete replay | 19,212 | 131.210 |
| h8_spin_r32_cold_atoms/round_4 | No complete replay | 19,400 | 176.547 |
| h8_spin_r32_cold_atoms/round_5 | No complete replay | 19,582 | 223.162 |
| h8_spin_r32_cold_atoms/round_6 | No complete replay | 19,760 | 271.509 |
| h8_spin_r32_cold_atoms/round_7 | 11.406913 | 19,930 | 328.332 |
| h8_spin_r32_exchange/round_0 | 5.276027 | 27,706 | 256.966 |
| h8_spin_r32_global_atoms_ready/round_0 | 2.886618 | 18,710 | 32.988 |
| h8_spin_r32_global_atoms_ready/round_1 | No complete replay | 18,800 | 66.289 |
| h8_spin_r32_global_atoms_ready/round_2 | No complete replay | 18,888 | 100.973 |
| h8_spin_r32_global_atoms_ready/round_3 | 2.880545 | 18,976 | 135.303 |
| h8_spin_r32_invariant/round_0 | 7.521336 | 28,672 | 69.005 |
| h8_spin_r32_invariant_dual16/round_0 | 4.859466 | 37,440 | 84.776 |
| h8_spin_r32_invariant_dual48/round_0 | 2.874956 | 54,932 | 226.354 |
| h8_spin_r32_response/round_0 | 14.007361 | 30,976 | 151.752 |
| h8_spin_r32_retained_response/round_0 | 606.549275 | 60,212 | 383.958 |
| h8_spin_r32_retained_response_scs/round_0 | 3.895871 | 60,212 | 266.634 |
| h8_spin_r32_screened_response/round_0 | 11.612197 | 30,976 | 51.800 |
| h8_spin_r32_screened_response/round_1 | 7.288539 | 39,040 | 114.020 |
| h8_spin_r32_screened_response/round_2 | 4464.024962 | 46,666 | 236.055 |
| h8_spin_r32_singlet/round_0 | 7.556042 | 28,672 | 261.163 |
| h8_spin_r64/round_0 | 16.310828 | 53,248 | 246.334 |

Every seed dependency is recorded. Initial H6 development reused the preceding prepared moments; later direct preparation was checked against those moments to 2.23e-16, and the fresh transfer used the direct constructor. No successful full cubic Gram proof or determinant reference state entered the new lower discovery. The old full cubic certificates remain available as preserved baselines.

## Other angles

The supplied spectral-filter sufficient theorem was derived and its scalar error kernel checked in 1,376 exact rational cases. No compact molecular filter and deterministic norm/defect certificate was obtained. The referenced toy archives were not located. The molecular spectral-degree estimates in `alternative_angles.json` are diagnostics using a crude enclosure, not required degrees or impossibility results.

A bounded molecular ratio probe found an exact frustrated triangle in each fixture, querying 49 selected determinant labels and two Hamiltonian columns. This rules out making all transitions negative by a diagonal sign gauge in that occupation basis. It does not rule out general signed routing or another basis. These selected-label diagnostics are separate from the enumeration-free lower constructors. No spin-chain toy calculation was substituted for a molecule.

One concrete mathematical difference remains: the strongest preserved H8 cubic certificate has 2,848 quartic terms in its number multiplier. The present direct-quartic paired constructor uses quadratic number multipliers and cancels sextic terms by pairing. The stronger proof can use fixed-number identities of degree six. This identifies an omitted source of flexibility; it does not prove that those terms are essential, nor an exact obstruction for the paired family. No full Gram factors from that older certificate were used to seed this pass.

## Accounting and limits

All 10,357 inherited sealed files match their original hashes. There are 81 recorded process receipts, totaling 5211.687 process-wall seconds. Pending: []. Failed processes: 9; their logs and costs are preserved.

The final algebra suite passed 24 tests, including negative-energy cases, mutated inputs and wedge indices, spin-sector refusals, exact spin projection, response inclusion, and the coherent-square constructor. Intermediate cold-search candidates without an exact receipt are shown as unverified; their numerical objective is not an accepted endpoint.

The tables report retained Gram entries, not certificate bytes, total discovery work, or verifier work. Dense projected maps remain a substantial cost. The H8 upper MPS is a frozen inherited input whose discovery and exact verification costs were measured in the previous campaign; it is not a free end-to-end solver component. New integrated bundles recheck the actual MPS upper inequality and new lower proof together, charging a tiny outward upper rounding allowance.

Each process receipt stores its command and input paths. For an independent lower replay, use `spin_replay` for the balanced-spin certificates or `spin_screen` with both spin pieces; these command-line checkers refuse numerical-library imports. `bundle` additionally rechecks the actual MPS upper inequality. Write new replay outputs outside this sealed campaign.

The certified statements solve the supplied finite electronic models to the stated intervals. General accuracy with increasing size, an undisclosed Nooterra algorithm, and experimental predictive advantage have not been established. A failed bounded search is not an impossibility theorem.
