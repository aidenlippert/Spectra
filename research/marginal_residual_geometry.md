# Complete pair pricing, residual rescaling, and an orbital test

The strongest exact complement bound in this checkpoint is **−6.530499343449799 Ha**. It comes from adding rational pair factors to the residual of an earlier four-coordinate certificate. The required reference threshold remains **−6.26489910104 Ha**, so the reference-gap construction is still incomplete.

The subsequent [compact residual and native-basis extension](marginal_compact_residual.md) reduces that proof's file size by55% while slightly strengthening the exact bound to−6.530499343207549 Ha. It also implements basis continuation with a matched first-LP control.

## The residual has its own useful geometry

Keep the original diagonal metric W fixed and write

\[
R=WH_QW-\sum_a\lambda_a v_av_a^T,\qquad m_i=W_{ii}^2.
\]

For any strictly positive rational vector t, each nonzero residual edge can be represented by a positive pair factor:

\[
z^{ij}=t_j e_i+\operatorname{sgn}(R_{ij})t_i e_j,
\qquad \mu_{ij}=\frac{|R_{ij}|}{t_it_j}.
\]

Subtracting these factors cancels every off-diagonal entry exactly. The remaining diagonal is

\[
d_i=R_{ii}-\sum_{j\ne i}|R_{ij}|\frac{t_j}{t_i}.
\]

Therefore H_Q≥γI follows whenever d_i≥γm_i for all i. The new t is independent of the fixed W. Its numerical proposal comes from the smallest generalized comparison-matrix eigenvector, C(R)t=γ diag(m)t, separately on each connected residual component. There is no nonlinear change of the original metric.

The implementation rounds positive t to bounded integers, normalizes each pair vector, and compensates its scalar coefficient by the square of that normalization. It merges duplicate directions, floors the combined coefficients to a common rational scale, reconstructs the entire residual exactly, and preserves the prior certificate if the proposed bound is weaker. This stays within rational factor width four: the new factors have support two.

| Source dictionary | Exact bound before polish | Exact bound after polish | Exported atoms, blocks 200 / 168 | Certificate bytes |
|---|---:|---:|---:|---:|
| Earlier rational active-pair run | −6.570674583080039 | **−6.530499343449799** | 10433 / 7665 | 4592165 |
| Complete-pair LP run | −6.567664856776602 | −6.5423889672490025 | 10837 / 8183 | 4812477 |
| Support-balanced dictionary restart | −6.563216121229496 | −6.546726489875678 | 10912 / 3281 | 3671835 |

The earlier branch polishes better even though its original DD bound is lower. Optimizing the original DD objective does not optimize the quality of this subsequent rescaling. Both branches are retained. The best polished block bounds are −6.530499343449799 and−6.479583689079265 Ha.

Artifacts: `results/marginal_h6/rational_atom_earlier_residual_polish` and `results/marginal_h6/rational_atom_residual_polish`. The implementing module is `experiments/marginal_residual_polish.py`; the existing versioned rational-atom replay checks the resulting certificates without NumPy or a solver.

The additional `results/marginal_h6/rational_atom_restart_residual_polish` branch is also weaker than the best earlier branch. Its168-state proposed polish loses against the prior exact bound, so the constructor retains that block's initial certificate. A numerical proposal or a newer search branch does not supersede a stronger accepted proof.

The pair factors admit a more compact description: the original atoms plus t determine them all. This checkpoint explicitly lists them. The subsequent [compact verifier](marginal_compact_residual.md) implements the displayed weighted residual inequality directly and avoids thousands of pair records.

## Discovery now distinguishes scanning from selection

Rational pricing can scan every two-coordinate support: 19,900 pairs for the200-state block and14,028 for the168-state block. Diagnostics separately report the number of tested supports, the minimum local eigenvalue, and the number of negative eigenspaces for support sizes two, three, and four. These are numerical diagnostics. Complete pair-support coverage is not an exact positivity certificate, and triples/quadruples remain sampled.

The complete-pair run imported the preceding rational certificate. It reached exact common bound −6.567664856776602 Ha, with block bounds −6.567664856776602 and−6.4866107156708575. Both blocks stopped at their240-second time limits. At one200-state round, 2,802 pair matrices were negative, but all256 selected directions came from larger supports. Scanning the pairs did not force the selection rule to use them.

An optional selection rule now reserves roughly one third of each batch for each support size, then fills unused slots with the strongest remaining violations. A deterministic fixture checks that this selects negative pairs even when quadruples rank ahead of all of them. Local eigenspaces that round to the same smaller-support rational direction are deduplicated before selection.

Search state is stored separately from the proof. The dictionary includes all loaded directions, zero-weight directions, and final priced directions that were not yet loaded. Imports bind the exact Hamiltonian, retained states, ordered complement blocks, metric weights, and vector family. Every pair touched by any dictionary direction is present before its constraint is reconstructed. The initial weighted proof remains separate, so dictionary ordering cannot corrupt its coefficients. This retains the cut pool; it does not persist the native solver basis or solution.

The complete-pair output is `results/marginal_h6/rational_atom_complete_pairs`; the support-balanced dictionary restart is `results/marginal_h6/rational_atom_dictionary_restart`.

The actual restart loads8,099 and7,771 directions while retaining weighted proofs with4,455 and3,419 atoms separately. Its first selected batches contain85 pairs,85 triples, and86 quadruples. The exported block bounds improve to−6.563216121229496 and−6.4467539341315065 Ha, with4,462 and3,281 positive atoms. Both block searches hit180-second limits after limited new rounds. The saved dictionaries grow to8,611 and9,307 directions. This confirms real dictionary reuse and support-balanced selection, while showing that the restart still pays for a fresh LP solve.

## A physical orbital change exposes a different tradeoff

The new orbital probe applies the exact rational rotation with cosine3/5 and sine4/5 to spatial orbitals2 and3, identically for both spins. Exact CAR, inverse recovery, and spin-commutation checks pass. The full400-dimensional spin-zero numerical spectra agree within3.4e−14 Ha. The original and rotated Hamiltonians each receive a fresh P32 selection using the existing constructor.

| Quantity | Original canonical orbitals | Rotated orbitals |
|---|---:|---:|
| Complement block dimensions | 200, 168 | 368 |
| Physical Q spectral floor, numerical | −6.26479910103936 | −6.238780309243709 |
| Comparison-matrix floor, numerical | −7.985247211264101 | −8.445200650256595 |

The physical complement floor improves by about0.0260 Ha, while the comparison bound worsens by about0.4600 Ha. A physically equivalent orbital basis can improve the actual gap and simultaneously make this certificate family harder to use. This single rotation neither optimizes the basis nor supplies a new energy certificate. The recorded result uses freshly selected P spaces for both Hamiltonians and distinguishes physical from comparison eigenvalues.

Files: `experiments/marginal_orbital_rotation.py` and `results/marginal_h6/orbital_rotation_probe.json`.

## Remaining scope

The exact positive decompositions still miss the required complement threshold. Complete support-four pricing and full-FW4 feasibility remain unresolved. All368 complement configurations remain explicit, and full replay references400 spin-sector configurations. The new results do not add a ground-energy interval to the accepted energy/witness ledger, which remains114. They identify a stronger constructive direction, a concrete selection problem, and a distinct basis tradeoff; they do not establish general representability or scalable quantum chemistry.

The full marginal regression suite passed **264 tests in392.049 seconds**. Focused checks cover missed pair support discovery, duplicate rounded directions, balanced selection, inactive dictionary fill pairs, source binding and malformed imports, fixed-metric residual rescaling, disconnected components, duplicate pair merging, exact orbital inverse recovery, CAR, and spin commutation.
