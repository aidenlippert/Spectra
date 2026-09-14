# Wave 2 upper-state refinement

Source implementation: `/Users/aidenlippert/Documents/Spectra/research/all_angles_20260913/selected_refinement/refine.py`, reused without modification. It performs sparse Ritz solves and signed accumulated residual selection with a diagonal-denominator score. Wave 1 used dense projected diagonalization and a single heuristic batch; Wave 2 incrementally expands the basis at requested support budgets.

Discovery uses only each fixture Hamiltonian. With no reference file supplied, the adapter chooses the lowest-diagonal determinant by scanning the finite fixture sector, then records that setup cost (`full_fock_enumeration: false` in the sparse algorithm receipt refers to no full matrix/Fock-vector construction; the diagonal reference scan is still recorded). The selected basis and the nonzero integer witness support differ because rationalization rounds small Ritz coefficients to zero.

| system | requested | actual basis | nonzero witness | upper Ha | width vs existing lower Ha | pass |
|---|---:|---:|---:|---:|---:|---|
| H4 | 64 | 15 | 15 | -3.6573313621 | 0.0096687469 | no |
| H4 | 128 | 20 | 20 | -3.6669999560 | 0.0000001531 | yes |
| H6 | 256 | 60 | 60 | -6.3046674663 | 0.0284461543 | no |
| H6 | 512 | 191 | 191 | -6.3320509382 | 0.0010626824 | yes |
| H6 | 1024 | 200 | 200 | -6.3330586262 | 0.0000549944 | yes |
| H8 | 256 | 185 | 185 | -9.2025215477 | 0.0535334665 | no |
| H8 | 512 | 512 | 512 | -9.2423737192 | 0.0136812950 | no |
| H8 | 1024 | 1024 | 1024 | -9.2522115086 | 0.0038435056 | no |
| H8 | 2048 | 2048 | 2048 | -9.2551715971 | 0.0008834171 | yes |

The H8 2048 width is below 0.0016 Ha and therefore passes. Its exact replay witness is [step_03_upper.json](/Users/aidenlippert/Documents/Spectra/results/wave2_20260913/upper_states/h8/step_03_upper.json), with receipt [receipt.json](/Users/aidenlippert/Documents/Spectra/results/wave2_20260913/upper_states/h8/receipt.json). The H6 1024 result also passes. H8 1024 was measured and fails, so this run does not claim a 1024-support pass.

The matched amplitude-score run is in `results/wave2_20260913/upper_states/h8_amplitude_hf2/`; the separate PT2 adapter is [refine_pt2.py](/Users/aidenlippert/Documents/Spectra/research/wave2_20260913/upper_states/refine_pt2.py), with results in `results/wave2_20260913/upper_states/h8_pt2_hf/`. Both start from the same explicit HF determinant (state 255), use budgets 32, 64, 128, 256, 512, 1024, one BLAS thread, identical rounding, and exact original-H replay. At 1024, amplitude gives `-9.25234502273209 Ha`; PT2 gives `-9.252483015042671 Ha`, an exact-Fraction difference of approximately `-0.00013799231 Ha`; PT2 is better. Source hashes: amplitude `ef307748076127ec080c5c9806da5b7c4374eee86eb923b0e164fc26f5186ac4`; PT2 `72dfc39e733b541ae5d6fe982ba6bc092f86f253aeb3261b0d4ec014f7c6ecfc`.

The 2048 H8 run took 1.77 s total, including sparse matrix assembly, Ritz solves, candidate generation, and exact replay; it accumulated 5,073 candidates and 321,030 projected nonzeros. All source/setup and replay costs are in the machine receipt. No FCI coefficients or prior Wave 1 coefficients were used as seeds. The method remains heuristic selected-CI discovery with an exact finite-H Rayleigh upper endpoint and no scaling theorem.
