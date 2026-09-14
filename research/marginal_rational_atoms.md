# Rational amplitudes and active pair constraints

The strongest independently replayed complement certificate from this extension is **−6.570674583080039 Ha**, still below the required reference threshold **−6.26489910104 Ha**. The constructor now discovers unequal rational amplitudes on two to four determinant coordinates. Both H6 blocks finish twelve pricing rounds within their individual 180-second budgets. These are valid positive decompositions, not an optimality result for factor width four.

Subsequent [complete pair pricing and residual rescaling](marginal_residual_geometry.md) improve the exact bound to **−6.530499343449799 Ha** and implement complete dictionary restarts. That report also records the first physical orbital-rotation diagnostic.

## Exact representation

For the fixed positive diagonal metric W, the exported proof is

\[
WH_QW=\sum_a\lambda_a v_av_a^T+R,
\qquad \lambda_a>0,\quad |\operatorname{supp}v_a|\le4.
\]

Every amplitude and weight is rational. The verifier reconstructs the entire residual and checks

\[
R_{ii}-\sum_{j\ne i}|R_{ij}|\ge\gamma W_{ii}^2.
\]

Thus R−γW² is diagonally dominant with nonnegative diagonal, hence positive semidefinite, and H_Q≥γI. The physical verifier checks complete complement coverage. Numerical eigenvectors and LP multipliers only propose factors; they are not trusted as certificates.

The new format `spin_rational_atom_complement_v1` explicitly permits unequal amplitudes. The legacy signed format continues to reject them. The local pricing eigensolver uses sampled principal matrices, rounds candidate amplitudes to rationals, and tests the rounded direction before accepting it for numerical discovery. Final acceptance depends on exact residual arithmetic.

With arbitrary amplitudes, the full factor-width-four cone is invariant under invertible positive diagonal congruence: W preserves every factor's support, and W inverse gives the converse. Consequently, changing W affects conditioning and the finite discovered dictionary; it cannot enlarge the full cone. This differs from the earlier restriction to flat signed vectors.

## Active pairs

The LP initially represents every nonzero Hamiltonian pair and every pair used by an imported atom. A new atom activates any additional pairs it touches, including pairs with zero Hamiltonian entry. An omitted pair is zero in both H and every selected atom, so its residual is identically zero. This preserves the finite-dictionary problem while reducing its initial number of pair variables. The full exact replay still checks every residual entry.

New pair-cap rows can be inserted between atom rows. The implementation tracks explicit atom row IDs, ensuring that pair-cap multipliers cannot become atom weights. A regression fixture forces a new zero-H pair to activate and checks that the same exact optimum is recovered in all-pair and active-pair modes.

## Actual H6 runs

| Run and block | Initial → final pair variables | Exported atoms | Exact bound, Ha | Discovery seconds |
|---|---:|---:|---:|---:|
| Signed, 200 | 6012 → 6027 of 19900 | 4093 | −6.578099830700481 | 180.714 |
| Signed, 168 | 4264 → 4939 of 14028 | 3062 | −6.491386413608022 | 65.511 |
| Rational, 200 | 6012 → 6573 of 19900 | 4259 | −6.570674583080039 | 126.661 |
| Rational, 168 | 4415 → 4941 of 14028 | 3163 | −6.489746592835086 | 85.320 |

The rational run imports the preceding signed certificate. Its 200-state block exports 261 unequal-amplitude vectors: 248 on four coordinates, 12 on three, and one on two. Its 168-state block exports 252: 239 on four coordinates and 13 on three. Remaining vectors use the legacy signed encoding. The rational certificate occupies 2,048,815 bytes.

The rational block bounds are exactly

\[
-\frac{379215390111231648997528301}{57713311672402500000000000},
\qquad
-\frac{122560936870840741455541157}{18885319344541500000000000}.
\]

Both blocks stop at the twelve-round budget, with negative sampled pricing values still present. These runs do not establish a stationary dictionary, complete pricing, or a factor-width-four obstruction. The timings are actual run observations, not a controlled speedup comparison: the starting dictionaries differ from earlier runs.

Artifacts are in `results/marginal_h6/signed_atom_active_pairs` and `results/marginal_h6/rational_atom_active_pairs`. Each contains a certificate, proposal history, receipt, and independent `python -S` replay. Every independent replay field matches its corresponding saved receipt. Both constructor receipts explicitly set `requested_target_certified` to false. The replay's `target_certified` field certifies the exported weaker gamma, not the requested reference threshold.

## What remains

The immediate unresolved issue is separating incomplete discovery from inadequate factor width. The current search samples supports, and a restart imports positive factors without preserving the full inactive cut dictionary or native solver basis. There is no proof yet that all four-coordinate directions can reach the needed threshold, or that they cannot.

A separate read-only provenance check establishes that H6 is still expressed in RHF canonical molecular orbitals (`experiments/marginal_h6_fixture.py`). No orbital rotation or localization experiment has been implemented. A physical orbital transformation is a distinct coordinate change from diagonal determinant weighting and is a plausible next test, not an accepted improvement.

An independent numerical probe tested a simple full-FW4 dual family. For any real vector u, B=4 diag(u_i²)−uuᵀ has every principal submatrix on at most four coordinates positive semidefinite by Cauchy–Schwarz. With A=H_Q−γI and positive D=diag(A), this family has negative trace against A precisely when the largest eigenvalue of D⁻¹ᐟ²AD⁻¹ᐟ² exceeds four. The actual block diagnostics were 2.4800516991085217 and 2.12195959956127. Thus this probe found no obstruction; it supplies neither an exact nonexistence certificate nor a constructive decomposition. Its explicitly numerical record is `results/marginal_h6/fw4_rankone_dual_probe.json`.

The cheapest omitted pricing coverage is all two-coordinate supports: 19,900 pairs in the200-state block and14,028 in the168-state block. Current rational pricing uses neighbor and random support selection. Exhaustive pair pricing could remove that omission, while still leaving three/four-support completeness open. No present evidence identifies rounding precision as the limiting issue: the last sampled negative values remain substantial.

All 368 Q configurations remain explicit; replay references 400 spin-sector configurations. This extension supplies no new ground-energy interval and no general scaling theorem. The accepted energy/witness ledger remains 114.
