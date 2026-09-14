# Asymmetric quartic discovery and solver precision

The asymmetric column generator now searches all 581 quartic dictionaries using only the surviving pair charges and global left/right exchange. New directions are exported as an identity/exchange pair of direct positive squares. They are not averaged over the broken flavor-permutation symmetry.

The first broad run reduced the exact interval to **1.4499595698358421e-5**, from the native-polished baseline of **1.929574335056606e-5**. Negative quartic directions therefore can improve this asymmetric proof, but the tested search has not located the quartic optimum. A separate [spin-sector Schur proof](marginal_schur_transfer_results.md) now gives a much tighter physical bound on this structured family.

## Search and exact acceptance

`experiments/marginal_asymmetric_columns.py` starts from existing positive-square directions. Its restricted LP optimizes their nonnegative weights and a fixed-number multiplier, with an explicit weighted residual penalty. HiGHS solves a minimization problem, so its equality marginal is negated to obtain the physical moment functional. A scalar test verifies this sign by adding an improving column.

For each quartic dictionary, the generator finds a negative eigenvector of its moment matrix, rounds it to a common rational denominator, and recomputes its coefficient column. The exporter stores the square and its global-exchange image as direct squares with half weights. Exact tests compare these columns with rational exchange averages and check that asymmetric components survive.

The construction has 866 projected coefficient rows, 316 multiplier directions, and 581 pricing dictionaries with maximum dimension 186. It still builds complete candidate maps. This is a finite-model discovery experiment, not a scaling result.

| Run | Added cuts, including resumed cuts | Best exact interval width | Observed total seconds |
|---|---:|---:|---:|
| Initial narrow SciPy run, 16 rounds x 8 cuts | 128 | 2.9806291662429152e-5 | 108.26 |
| Native precision, tight dual centering | 384 | 1.9420766410602812e-5 | 322.58 |
| Relaxed dual centering, eight further rounds | 576 | 1.940764878057709e-5 | 496.70 |
| Broad pricing, four further rounds | 2347 | **1.4499595698358421e-5** | 299.54 |

These runs use different starting cuts and settings, so their times are not a controlled speed comparison. Their directories are under `results/marginal_asymmetric_columns/` with names `epsilon001_16`, `epsilon001_centered_fixed`, `epsilon001_centered_relaxed`, and `epsilon001_breadth4`.

The broad run adds the lowest negative direction from every violated dictionary in each round. It finishes with 13 orbit and 395 direct squares, lower 3.280990278287705, and residual L1 about 8.64e-10. Independent standard-library replay checks the lower certificate and recomputes the upper witness against the exact epsilon=1/1000 Hamiltonian.

## What caused the earlier plateau

The narrow run concentrated 384 cuts in only 41 dictionaries. These are different directions; repeated dictionary indices do not establish duplicate rays. There were still more than 500 negative dictionaries on many iterations. Broadening the search improved the energy after four rounds, while narrow dual centering produced little gain. This identifies limited coverage as one contributor; it does not prove a unique cause or convergence rate.

The tight centering QP was accepted on seven of seventeen LP states; rejected proposals fell back to the raw LP dual. Widening its energy-face tolerance to 1e-5 produced more stable numerical eigenvalues, but little energy improvement. The QP is a search heuristic, and its floating feasibility is never accepted as an exact moment witness.

The weight cap is observable in the receipts. Raising the cap on the same 1,196 raw directions from 100 to 1,000 and 10,000 did not close the gap and increased exact export error. A nearly physical-null polynomial reaches the cap. Its small but nonzero rational action is accounted for by exact replay, rather than discarded as mathematically zero.

The production LP now uses native HiGHS with internal small-matrix threshold 1e-12. A controlled same-backend comparison isolated the previous coefficient-dropping loss. See [precision measurements](marginal_lp_precision_probe_results.md). Neither the threshold fix nor a solver `Optimal` status certifies an exact LP or quartic-cone optimum.

## Remaining cone question

The search still needs an energy-relevant set of directions that stays manageable as the symmetry and system size change. The failed full quartic numerical reference does not establish a hierarchy gap. Persistent negative eigenvalues do not establish one either: those duals have not been repaired into exact feasible quartic witnesses. The Schur certificate resolves a physical precision target for this specific family while leaving that general discovery question open.

A further broad continuation (`epsilon001_breadth8`) reached three completed update rounds, with numerical objective 3.2809904829. The delegated process was interrupted after a recorded 747 seconds, before its planned 900-second deadline. Its final exact export did not complete. Existing checkpoints and the corrected process record are preserved; no new accepted interval or solver-timeout conclusion is drawn from this partial run.
