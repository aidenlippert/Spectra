# Measured campaign tables

All energy widths below come from exact accepted endpoints. A process that
finished successfully can still miss the 1.6 mHa accuracy target. Times are
single measured runs on the local M1 host, without an estimated variance.

| Fresh run | Width, mHa | Complete seconds | Target met | Selected Gram entries | Peak child RSS, MB |
|---|---:|---:|---|---:|---:|
| [h10_matched_direct](../../results/interacting_scaling_20260915/cold/h10_matched_direct/result.json) | 1.557099696 | 1837.661 | True | 702,332 | 1012.3 |
| [h10_matched_reference](../../results/interacting_scaling_20260915/cold/h10_matched_reference/result.json) | 10.754481671 | 1132.241 | False | 1,287,700 | 805.9 |
| [h10_matched_reference_extended](../../results/interacting_scaling_20260915/cold/h10_matched_reference_extended/result.json) | 0.479085103 | 1284.551 | True | 1,287,700 | 820.6 |
| [h12_heldout](../../results/interacting_scaling_20260915/cold/h12_heldout/result.json) | No accepted interval | 1502.800 | False | — | 1086.4 |
| [h4_631g_heldout](../../results/interacting_scaling_20260915/cold/h4_631g_heldout/result.json) | 0.156525949 | 328.437 | True | 181,268 | 694.5 |
| [h4_cold_integration](../../results/interacting_scaling_20260915/cold/h4_cold_integration/result.json) | 0.003055344 | 19.278 | True | 15,740 | 205.2 |
| [h4_reference_integration](../../results/interacting_scaling_20260915/cold/h4_reference_integration/result.json) | 0.004060572 | 11.791 | True | 6,368 | 200.5 |
| [h8_matched_direct](../../results/interacting_scaling_20260915/cold/h8_matched_direct/result.json) | 1.030971387 | 463.317 | True | 181,268 | 616.3 |
| [h8_matched_reference](../../results/interacting_scaling_20260915/cold/h8_matched_reference/result.json) | 0.227658227 | 324.278 | True | 343,424 | 594.4 |
| [water_heldout](../../results/interacting_scaling_20260915/cold/water_heldout/result.json) | 0.084534412 | 182.389 | True | 111,356 | 439.4 |

## Representation and retained witnesses

Selected Gram entries above describe the accepted representation. The sum
below also counts earlier main-proof levels prepared during that cold run;
it is not a peak allocation or arithmetic-work estimate. Magnetic-screen
work is charged separately in the component ledger. Retained witness size
includes the main and nonsinglet proof, MPS/upper and orbital transform,
but excludes shared Hamiltonian descriptions, checker code and replay logs.

| Run | Sum of prepared main Gram entries | Largest accepted block | Coefficient rows | Map nonzeros | Retained witness, MB | MPS nonzero integer entries |
|---|---:|---:|---:|---:|---:|---:|
| h10_matched_direct | 1,117,888 | 422 | 143,085 | 1,891,122 | 15.343 | 3904 |
| h10_matched_reference | 1,287,700 | 505 | 112,701 | 4,369,663 | 11.515 | 3904 |
| h10_matched_reference_extended | 1,287,700 | 505 | 112,701 | 4,369,663 | 11.528 | 3904 |
| h4_631g_heldout | 181,268 | 220 | 34,133 | 426,877 | 6.356 | 2334 |
| h4_cold_integration | 15,740 | 56 | 819 | 33,849 | 0.643 | 108 |
| h4_reference_integration | 6,368 | 34 | 435 | 14,338 | 0.352 | 108 |
| h8_matched_direct | 181,268 | 220 | 34,133 | 426,877 | 5.813 | 2400 |
| h8_matched_reference | 343,424 | 260 | 28,461 | 1,084,148 | 4.616 | 2400 |
| water_heldout | 111,356 | 170 | 18,056 | 259,080 | 3.897 | 1527 |

## Cold component costs

These child-stage sums are included within the fresh clock above; they are
not additional expenses. Rejected proof levels remain included.

| Run | Input/orbitals, s | State/upper, s | Magnetic screen, s | Proof preparation, s | Proof solves, s | Exact replay, s | Other, s |
|---|---:|---:|---:|---:|---:|---:|---:|
| h10_matched_direct | 5.513 | 265.709 | 63.353 | 189.827 | 907.615 | 403.457 | 1.321 |
| h10_matched_reference | 5.619 | 266.594 | 51.989 | 100.289 | 452.846 | 227.515 | 26.782 |
| h10_matched_reference_extended | 5.445 | 265.729 | 51.777 | 99.525 | 606.102 | 228.768 | 26.619 |
| h12_heldout | 12.328 | 835.220 | 654.488 | 0.000 | 0.000 | 0.000 | 0.000 |
| h4_631g_heldout | 4.073 | 131.260 | 19.035 | 28.316 | 94.371 | 50.720 | 0.083 |
| h4_cold_integration | 2.075 | 1.865 | 2.996 | 2.442 | 8.277 | 1.040 | 0.048 |
| h4_reference_integration | 1.271 | 1.447 | 2.816 | 1.725 | 3.244 | 0.621 | 0.088 |
| h8_matched_direct | 2.418 | 132.507 | 21.550 | 28.366 | 225.128 | 52.687 | 0.089 |
| h8_matched_reference | 2.462 | 131.304 | 17.911 | 27.446 | 84.814 | 53.503 | 6.267 |
| water_heldout | 1.946 | 27.675 | 13.665 | 16.853 | 98.032 | 23.583 | 0.090 |

## Coupling paths

Intermediate Hamiltonians are diagnostic models. The endpoint at lambda=1
includes the exact original-model rotation allowance. Each path includes
state continuation and failed proof attempts, but source integral/orbital
preparation is additional and remains in the campaign ledger.

### h4_coupling

Measured continuation clock: 65.604 seconds; [receipt](../../results/interacting_scaling_20260915/coupling/h4_coupling/results.json).

| Coupling | Width, mHa | Target met | Selected case |
|---|---:|---|---|
| 0 | 0.002783053 | True | [h4_coupling_l0](../../results/interacting_scaling_20260915/cases/h4_coupling_l0/design.json) |
| 1/4 | 0.007518891 | True | [h4_coupling_l025](../../results/interacting_scaling_20260915/cases/h4_coupling_l025/design.json) |
| 1/2 | 0.009877641 | True | [h4_coupling_l05](../../results/interacting_scaling_20260915/cases/h4_coupling_l05/design.json) |
| 1 | 0.012779484 | True | [h4_coupling_l1](../../results/interacting_scaling_20260915/cases/h4_coupling_l1/design.json) |

### h8_adaptive_coupling

Measured continuation clock: 1929.606 seconds; [receipt](../../results/interacting_scaling_20260915/coupling/h8_adaptive_coupling/results.json).

| Coupling | Width, mHa | Target met | Selected case |
|---|---:|---|---|
| 0 | 0.004003898 | True | [h8_adaptive_coupling_l0_level0](../../results/interacting_scaling_20260915/cases/h8_adaptive_coupling_l0_level0/design.json) |
| 1/4 | 0.018295612 | True | [h8_adaptive_coupling_l025_level0](../../results/interacting_scaling_20260915/cases/h8_adaptive_coupling_l025_level0/design.json) |
| 1/2 | 0.059878041 | True | [h8_adaptive_coupling_l05_level0](../../results/interacting_scaling_20260915/cases/h8_adaptive_coupling_l05_level0/design.json) |
| 1 | 1.031020918 | True | [h8_adaptive_coupling_l1_level0](../../results/interacting_scaling_20260915/cases/h8_adaptive_coupling_l1_level0/design.json) |

## H12 post hoc diagnostics

These continuations do not replace the failed frozen H12 run. They reuse
its state and selected maps; its complete 1,502.800-second clock remains
additional. Parent-clock sums below exclude intervening work on other
cases and are not a newly timed cold pipeline.

- [h12_recovery_result.json](../../results/interacting_scaling_20260915/h12_recovery_result.json): 452.545612138 mHa; 1304.418 additional seconds.
  Recorded diagnostic parent-clock sum through this point: 2107.001 seconds; source cold run additional.
- [h12_main_refinement_result.json](../../results/interacting_scaling_20260915/h12_main_refinement_result.json): 12.590172569 mHa; 906.084 additional seconds.
  Recorded diagnostic parent-clock sum through this point: 3013.085 seconds; source cold run additional.

## Disconnected controls

| Control | Width, mHa | Local Fock labels over every charge | Global determinants enumerated |
|---|---:|---:|---:|
| [h4](../../results/interacting_scaling_20260915/fragment_controls/h4/replay.json) | 0.000020000 | 32 | 0 |
| [h8](../../results/interacting_scaling_20260915/fragment_controls/h8/replay.json) | 0.000040001 | 64 | 0 |

## External algorithm controls

Numerical FCI estimates and exact enumerated certificates are separate
comparisons. FCI runs after proof discovery on the pre-rounding integral
model; it does not supply a rigorous two-sided certificate.

| Exact enumerated control | Fresh complete seconds | Accepted exact width, mHa | Literal original input match |
|---|---:|---:|---|
| [h4_631g_heldout_enumerated_control](../../results/interacting_scaling_20260915/models/h4_631g_heldout_enumerated_control/baseline_result.json) | 20.559 | 0.799999976 | True |
| [h4_cold_integration_enumerated_control](../../results/interacting_scaling_20260915/models/h4_cold_integration_enumerated_control/baseline_result.json) | 0.861 | 0.800000036 | True |
| [water_heldout_enumerated_control](../../results/interacting_scaling_20260915/models/water_heldout_enumerated_control/baseline_result.json) | 5.001 | 0.799999972 | True |

Numerical FCI timing receipts:

- [h10_matched_direct_numerical_fci_solve](../../results/interacting_scaling_20260915/runs/h10_matched_direct_numerical_fci_solve.json): 0.795 seconds, passed.
- [h12_heldout_numerical_fci_solve](../../results/interacting_scaling_20260915/runs/h12_heldout_numerical_fci_solve.json): 7.242 seconds, passed.
- [h4_631g_heldout_numerical_fci_solve](../../results/interacting_scaling_20260915/runs/h4_631g_heldout_numerical_fci_solve.json): 0.504 seconds, passed.
- [h8_matched_direct_numerical_fci_solve](../../results/interacting_scaling_20260915/runs/h8_matched_direct_numerical_fci_solve.json): 1.151 seconds, passed.
- [water_heldout_numerical_fci_solve](../../results/interacting_scaling_20260915/runs/water_heldout_numerical_fci_solve.json): 0.496 seconds, passed.

## Complete campaign ledger

There were 342 measured child attempts, totaling 15170.247 seconds, including failed attempts, tests, and dependency installation. This is a sum of measured compute stages, not one fresh solver run.

Editorial work and routine file inspection are excluded from that compute ledger.
Cold orchestration clocks are reported separately and must not be added again.
No new cloud instances or external spending were used.

See [the complete accounting record](../../results/interacting_scaling_20260915/accounting.json) for exact receipt hashes, failed process attempts, accuracy misses and representation sizes.
