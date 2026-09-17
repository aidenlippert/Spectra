# Verified finite-time molecular control — scope and outcome

A new exact calculation steers the **specified rational H8 model from its supplied MPS initial state** to a negative orbital population contrast. Two finite-time protocols meet D<=-0.6, including nonzero bounds on waveform deviations, mixed-state initialization error, and local phase-flip dephasing.

This is an intervention-valid reduced dynamical calculation and a finite-model reachability result. It is **not** the universal many-body shortcut or physical-compiler breakthrough requested. The reduced basis was discovered from full balanced-sector simulations; both discovery and acceptance enumerate 4,900 basis states and build 903,620 Hamiltonian entries. That is an explicit unresolved computational dependency, not hidden preprocessing.

## Physical model and preparation assumptions

The Hamiltonian is the unchanged rational STO-3G H8 electronic model with eight electrons, in the balanced 4-alpha/4-beta sector. The initial state is the actual rounded MPS from the accepted H8 bundle, normalized by its exact norm. It is not assumed to be an exact ground state, and no physical method for preparing it is demonstrated.

The controls are coefficients of

D=n6+n7-n10-n11

W=a6^dagger a10+a10^dagger a6+a7^dagger a11+a11^dagger a7.

They are noncommuting, number-preserving operators. They are not an experimentally calibrated laser, field, geometry, or synthesis interface. All molecular interactions remain present. The nuclei and orbital basis are fixed.

The exact initial D expectation is approximately **1.390133812147592**. A control reverses that contrast; this is not a statement about macroscopic charge transport, a new material, superconductivity, or reaction yield.

## Accepted dynamic results

| Protocol | Time, atomic units | Reduced coordinates | Nominal final D interval | Robust final D interval | Target D<=-0.6 |
|---|---:|---:|---:|---:|---|
| Four-piece fast control | 2 | 24 | [-0.673065344, -0.664524560] | [-0.711065344, -0.626524560] | Proven |
| Constant slower control | 21 | 64 | [-0.644952032, -0.639180083] | [-0.681192032, -0.602940083] | Proven |
| No-control comparison | 21 | 32 | [1.360886861, 1.402259972] | same | Proven D>=1 |
| Eight-coordinate fast approximation | 2 | 8 | [-0.956154750, -0.380970893] | [-0.994154750, -0.342970893] | Not proved |

Displayed bounds are conservatively rounded; exact endpoints are in the replay receipts. Time 21 atomic units is approximately 0.508 femtoseconds; time 2 is approximately 0.0484 femtoseconds. The short durations and coefficient-level actuation are important limitations.

The four fast phases last 1/2 atomic time unit each. Their (u,v) controls, in Hartree, are:

(-88687/1000000,1/2), (1/2,1/2), (1/2,1/2), (311663/1000000,1/2).

The longer control is u=1/5, v=1/10 for 21 atomic time units.

## Continuous robustness, not sampled confidence

The fast protocol permits |delta u(t)|,|delta v(t)|<=0.001 Ha pointwise. The longer one permits 0.00009 Ha. Both allow initial trace distance up to 0.0005 from the normalized MPS in the same balanced sector. Both allow local phase-flip rates with sum-integral at most 0.001 over the complete evolution.

All time-dependent waveforms and nonnegative rate profiles within these bounds are covered by a norm theorem. The conditions do not state experimentally achieved calibration or preparation fidelity. This is a specified dephasing model, not arbitrary environmental physics.

The exact nominal state-error upper bounds are approximately 0.001067598 for the 24-coordinate path and 0.000721493 for the 64-coordinate path. Integration errors, basis rounding, initial representation error, and jumps between polynomial segments are included.

## What was new in this pass

- A reduced quantum trajectory accompanied by exact small Gram matrices of its **actual molecular equation defect**.
- Continuous-time acceptance using exact polynomial integrals, not sampled residuals.
- A guaranteed finite-time change in the model, rather than a static comparison of ground energies.
- A uniform error budget combining waveform, initialization, and dephasing uncertainty.
- A mandatory insufficient-proof outcome when an eight-coordinate representation cannot justify the target.

The residual theorem, quantum control, POD, and offline/online reduction have established predecessors. No priority claim is made for them. The new implementation and verified molecular output do not establish a new universally efficient quantum algorithm.

## Discovery failures and dependencies

An initial basis trained on constant-control trajectories over 20 time units gave a state-error estimate above one for an untrained switching protocol, even at 128 coordinates. The numerical reference state error was about 0.417. It was rejected as a useful control certificate. The source log is retained.

The successful bases were adapted to the selected interventions using numerical full-model trajectories. The fast control was proposed with local L-BFGS-B optimization in a reduced model and then re-adapted. The slower control came from an exploratory grid. Global optimality, minimal time, minimal field strength, and broad control-space coverage are not proved.

The long basis has 313,600 stored integer coefficients, larger than a single reference state. Cheap online matrices do not imply cheap discovery. Removing the full-state snapshot and embedding dependency remains open.

The complete control theorem does not use an inherited energy lower certificate. The inherited MPS supplies a specified initial condition and remains a state-discovery dependency. Its energy certificate does not establish its fidelity to the true ground state, and no such fidelity is assumed.

## Verification and comparisons

The stand-alone standard-library dynamics replay ran twice successfully, around 30 seconds for four cases after reading the original model and MPS. A complete packaged replay also runs 20 focused tests; its authoritative elapsed time is recorded in RESULT.json after that run.

The independent numerical expm_multiply calculation gives final D approximately -0.6687977225 and -0.6420483090 for the two protocols, within the exact intervals. It is a sanity check, not the acceptance argument. Additional control-box corners and an independently exponentiated two-level Lindblad example were checked numerically.

No complete cost advantage over optimized ordinary propagation is established. In particular, our exact result's shared preparation and validation costs must be included. Reusing small defect matrices reduces per-query dimension; it is not evidence of general scaling.

No H12 solve, new geometry, thermal phase calculation, laboratory experiment, or synthesis policy was completed. No repository write, external provisioning, paid compute, or package installation occurred.
