# H8 target closed: 0.7674481354629528 mHa

This bundle contains a **new, independently accepted lower certificate for the actual uploaded H8 rational Hamiltonian**, the unchanged rational MPS upper, the separate inherited nonsinglet certificate, the exact checker source, and discovery code/evidence.

## Exact result

- Full fixed-N lower: `-277669368076014353 / 30000000000000000` Ha.
- Upper: the same exact rational MPS expectation used in the supplied comparison, approximately `-9.2548781543983488138677` Ha. Its complete fraction is in `RESULT.json` and the receipts.
- Certified interval width: **0.7674481354629528 mHa**, below **1.6 mHa**.
- Supplied compact baseline: **2.8491970326490175 mHa**, independently replayed before discovery.
- Model: electronic rational STO-3G H8 chain, 1.4 Angstrom spacing, 16 spin orbitals, 8 electrons. This is not a statement about basis-set, geometry, nuclear-motion, or experimental error.

The final singlet proof contains **203 factor rows**. An earlier accepted export had 580 rows; 377 tiny rows were removed with a total discarded-square norm bound below 1.6e-12 Ha, then the smaller proof was independently re-expanded and accepted. Both versions are retained.

The new lower was accepted by the **unmodified supplied CAR/SU(2) checker**, using the ordinary coefficient L1 remainder. It does not inherit the stronger full-cubic lower or its residual witness. The separate nonsinglet certificate and MPS upper remain inherited components and were independently replayed. No fixed-N determinant Hamiltonian or state vector was enumerated in this construction or acceptance.

## Replay the entire interval

Unzip this bundle, enter its directory, and run with Python:

```sh
python3 -B -S replay.py --out ../H8-new-independent-replay
```

Only the standard library is required. The output directory must not already exist. The command hashes the bundle inputs/source, recomputes the actual MPS expectation, expands and checks the new positive squares and residual, verifies the nonsinglet piece, charges the spin defect, and calculates the interval. It never accepts a stored numerical optimizer objective. The included receipts are comparison data, not substitutes for those computations.

The final standalone replay of this exact packaged proof passed in **90.82 seconds** in this runtime.

The earlier separate accepting runs in this runtime took 32.16 s for the MPS and 55.57 s for the final complete lower. They are observations, not runtime guarantees. Running both proofs can take a few minutes on another machine.

## What changed

The prior compact frames excluded useful operator couplings. The successful construction **does not heuristically truncate the mixed cubic dictionaries**. It reorganizes them by spin into highest-weight doublet and quartet multiplicity spaces. Exact, sparse raising maps give rational doublet bases; the coefficient map then stays sparse.

For each particle charge and conserved spatial parity, four magnetic-component dictionaries of dimensions `112, 368, 368, 112` are represented using one `256`-dimensional spin-1/2 multiplicity Gram block and one `112`-dimensional spin-3/2 block. Cross terms within these multiplicity spaces remain available. The already admitted spin-average and sector ideal rules are unchanged.

The resulting numerical solve has:

- 335,168 singlet Gram matrix entries; 353,296 including the inherited nonsinglet construction.
- 1,370,610 nonzeros in its projected coefficient maps.
- An 8,533-dimensional normal matrix with **221,779 stored nonzeros**.
- Sparse LU factors with 954,952 nonzeros; the successful path does not use the supplied approximately 582 MB dense normal-factor cache.

This is a tradeoff, not an across-the-board compression win: the old compact construction used 115,032 complete Gram entries. The new construction uses more than that failed compact model, while substantially fewer than the roughly 1.2 million-entry large reference. The sparse representation, not the smallest possible rank, enabled the successful search. No optimality or impossibility claim about the old compact spans follows.

## Discovery dependencies and costs

The successful solve used the uploaded compact candidate as a warm start, transported into the spin-adapted basis. It also used explicitly contracted three-body moments of the uploaded MPS to initialize a physical proposal-side functional. It did not read the stronger full-cubic factor rows as a teacher.

The winning path has a recorded stage-1 checkpoint at 117.8004 s, followed by 58.4562 s of refinement and export. Building the sparse representation/normal matrix from the supplied full coefficient maps took 1.4540 s. These figures exclude original map/MPS/nonsinglet discovery and numerous exploratory attempts in this session. Available failed/stopped-run logs are preserved in `discovery_evidence/run_records`. There is **no claim of an audited total fresh-problem runtime or a matched end-to-end speedup**.

Two other branches tested larger MPS-guided spans and various numerical solvers. They did not supply the accepted bound here. A test initially asserted 1,920 mixed words instead of 3,840 across both particle charges; all word comparisons had succeeded. The corrected count test passed, and both logs are retained.

## Optional numerical reproduction in an existing Spectra repository

The `discovery/` scripts are prototype numerical code. They need NumPy, SciPy, and (for MPS moments) Numba, plus the original handoff's prepared H8 coefficient maps. They are separate from the standard-library accepting path. Set a fresh output directory:

```sh
export SPECTRA_ROOT=/path/to/Spectra
export SPECTRA_OUTPUT=/path/to/new-h8-discovery
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMBA_CACHE_DIR="$SPECTRA_OUTPUT/numba_cache"
mkdir -p "$SPECTRA_OUTPUT"
python3 -B discovery/mps_six_moments.py
python3 -B discovery/moment_probe.py
python3 -B discovery/build_spin_sparse.py
python3 -B discovery/solve_spin_sparse.py stage1 --seconds 130 --mu 2
python3 -B discovery/solve_spin_sparse.py stage2 --seconds 180 --mu .03 \
  --restart "$SPECTRA_OUTPUT/stage1/checkpoint.npz"
```

A fresh numerical trajectory need not reproduce the identical candidate or meet the target in the same wall-time budget. Independently replay its exported certificate with the original `spin_screen` checker. The accepted certificate included here is independently reproducible without repeating numerical discovery.

The dictionary decomposition is specialized to the supplied H8 frame/group organization. Generalizing or testing it on other models requires additional work. Spin decomposition and semidefinite optimization are not claimed as new mathematics.

## Scope

This closes the **specified H8 numerical accuracy task** and supplies a new successful lower construction. It is not a proof of polynomial scaling, a general solution for arbitrary molecular Hamiltonians, a transfer result, or experimental validation. No original project file was modified, and nothing was uploaded to GitHub or published externally.
