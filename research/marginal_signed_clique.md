# Cancellation-enabled clique squares: exact obstruction and next directions

The saved numerical proposal has now become an independently replayable rational dual. It rules out the required H6 reference threshold even when balanced clique squares may overspend Hamiltonian edges and opposite-sign residual edge squares compensate.

The exact ceiling on the168-state support is

\[
\gamma\le -\frac{389793629361447613489}{57777922000000000000}
=-6.746411360405928\;\mathrm{Ha},
\]

below the reference target−6.26489910104 Ha. The prior sign-matched residual ceiling was−6.798478663107022 Ha. Permitting cancellation enlarges the primal family and therefore raises its best threshold; it still does not reach the target in these metrics. This result does not exclude arbitrary signed three/four-coordinate atoms or arbitrary local amplitudes.

## The exact dual

Let A'=WHW, m_i=W_ii², b_i=A'_ii−sum_{j!=i}|A'_ij|, and c_ij=|A'_ij|. For nonnegative integer weights y_i,z_ij, define

\[
B_{ii}=y_i,\qquad
B_{ij}=\tfrac12\operatorname{sign}(H_{ij})(z_{ij}-y_i-y_j)
\]

on nonzero Hamiltonian edges, with B_ij=0 on absent edges. The signed residual edge constraints are

\[
0\le z_{ij}\le2(y_i+y_j).
\]

They enforce nonnegative pairing with both (e_i+e_j)(e_i+e_j)^T and (e_i−e_j)(e_i−e_j)^T, regardless of the Hamiltonian edge orientation. Nonnegative diagonal weights handle diagonal residual squares. The balanced-clique constraints remain

\[
\sum_{e\subset C}z_e\ge(|C|-2)\sum_{i\in C}y_i.
\]

Consequently every decomposition using those balanced support-three/four squares plus a real diagonally dominant residual obeys

\[
\gamma\le
\frac{\operatorname{Tr}(A'B)}{\operatorname{Tr}(W^2B)}
=\frac{\sum_i b_iy_i+\sum_e c_ez_e}{\sum_i m_iy_i}.
\]

This is an upper bound on the threshold achievable by the declared proof family, not an upper spectral bound for the physical Hamiltonian.

## Rational repair and support

The numerical node and edge weights are rounded at scale10^9. Edge weights are clipped to the exact interval[0,2(y_i+y_j)]. If a clique inequality is violated, add an integer t to every y_i and2t to every z_e. This adds tI to B. Every r-coordinate clique slack increases by rt and every opposite-edge slack increases by2t, so a sufficiently large exact t repairs all constraints without invalidating earlier ones. On this saved H6 proposal, t=0: rounding and exact checking already suffice.

Only y at determinants3192 and3252 is nonzero, and those two values are equal. There are76 nonzero edge weights. Replay checks32,013 balanced cliques and4209 two-sign edge inequalities, evaluates168 source actions, and references200 determinants. The certificate is280,785 bytes.

The dual is extended by zero outside the declared Q support. A global balanced clique intersects that support in at most four coordinates; its restriction is a checked balanced clique, a checked edge square, or a nonnegative diagonal. Thus complete Q component coverage is unnecessary for this negative result. The metric may be extended positively outside the support. A positive complement proof, by contrast, still needs complete complement coverage.

The previous unit edge-cover obstruction is rejected by the new verifier: its six covered edges away from its active row violate the opposite-sign pair gate. This confirms that the expanded family is being checked rather than merely relabeling the old result.

```sh
python -S -m experiments.marginal_signed_clique \
  --verify results/marginal_h6/clique_signed_residual_obstruction/certificate.json
```

The independent standard-library replay matches every saved receipt field. The accepted energy/witness ledger remains114. The unresolved construction must change the local atom signs or amplitudes while preserving an inexpensive exact residual certificate; using a dense PSD residual would restore the old dense-factor dependency.


## A square outside the balanced-clique dictionary

Exact replay finds a four-coordinate separator on determinant states[123,183,3192,3252], with transformed amplitudes[+1,−1,+1,−1]. Its unnormalized pairing with the accepted dual is−299554717546, or−2 times either active node weight. After normalizing the physical positive atom and physical dual to unit trace, the pairing is

\[
-\frac{125000000000000000}{250834572067659521}
=-0.4983364094096359.
\]

The pairs(183,3192) and(123,3252) have zero Hamiltonian entries. Thus a constructor restricted to complete Hamiltonian cliques cannot introduce this direction. Adding it requires tracking off-diagonal fill and cancellation even on zero Hamiltonian edges.

```sh
python -S -m experiments.marginal_signed_clique \
  --verify-separator results/marginal_h6/signed_clique_separator/certificate.json
```

This replay checks the complete parent dual as well as the positive separator. Rejecting this one dual does not establish a sufficient enlarged certificate.

## New constructive family

`marginal_signed_atoms.py` now verifies explicit signed rank-one atoms with three/four-coordinate support, permitting several distinct sign patterns on the same support. It subtracts these atoms from WHW exactly, then accepts only when the entire residual minus gamma W² is diagonally dominant with nonnegative diagonal. There is no dense PSD residual factor.

The native proposer uses a pair variable for **every** coordinate pair, including zero Hamiltonian entries. Atom row multipliers propose nonnegative square coefficients. Row-based, eigenvector-based, and fixed-seed random supports supply bounded heuristic pricing, with all anchored sign patterns tested on each selected support. The final rational residual determines the bound independently of the numerical LP status. A negative-threshold rank-one control verifies the row-multiplier sign convention and the metric target; an independent test permits exact cancellation between different sign patterns on the same support.

This is still an explicit configuration construction and heuristic atom search. No global pricing theorem, compact sufficient H6 gap, or all-size bound is implied by these implementation checks.


A separate boundary-geometry avenue was checked analytically. If A_gamma=W(H−gamma I)W is nonzero positive semidefinite, has kernel u, and every signed support-at-most-four vector s has nonzero u^T s, then epsilon=min_s (u^T s)^2/|supp(s)| is positive. B=uu^T−epsilon I is nonnegative on every flat signed sparse atom, yet Tr(A_gamma B)=−epsilon Tr(A_gamma)<0. This would obstruct the entire flat sparse family at that boundary. It is a conditional theorem, not an H6 certificate: the numerical H6 metric-168 target matrix has smallest eigenvalue about2.57e−4, and sampled four-component signed sums of its low eigenvector are as small as5e−17. Symmetry-related exact cancellations remain possible. No positive epsilon or exact kernel is established here; this route supplies no additional accepted obstruction.


## First H6 construction

The first bounded all-pair run exports a common exact complement threshold−6.584937444700064 Ha. It uses5107 atoms on the200-state block and3458 on the168-state block. The respective exact block bounds are−6.584937444700064 and−6.501896415893476 Ha. The168-state proof contains46 positive atoms whose supports are not balanced Hamiltonian cliques; its improved bound therefore concretely leaves the previously excluded family.

The200-state initial solve reaches its180-second limit before an optimal LP permits the next pricing round. The168-state run completes several adaptive rounds, then also reaches its180-second budget. Neither stopping event proves a cone limit. Recorded per-block proposal/export times are181.159 and180.756 seconds, including exact export overhead. The certificate is2,271,282 bytes and replays368 source actions, referencing400 determinants. It eliminates the dense positivity factor from this weaker gap proof, but neither reduces the configuration frontier nor supplies a smaller byte representation than all earlier proofs.

```sh
python -S -m experiments.marginal_signed_atoms \
  --verify results/marginal_h6/signed_atom_adaptive_gap/certificate.json
```

A subsequent active-atom restart imports only the positive directions and their exact weights, bound to the same Hamiltonian, P, ordered blocks, and metric. It recomputes every residual and preserves the original exact bound if numerical proposals deteriorate. The restart is a new bounded search, not continuation of a saved HiGHS basis: inactive directions and the native basis are not currently persisted. This distinction matters for the remaining optimization cost.


The restart finishes with exact common bound**−6.578099830702214 Ha**. Its block bounds are−6.578099830702214 and−6.501896415859300 Ha, using4315 and2960 atoms. The certificate is1,972,929 bytes. Both block searches reach their180-second budgets after limited new pricing. On the200-state block, a later time-limited multiplier proposal would give a much worse residual bound near−23.25 Ha; it is not selected, and the earlier better exact result is preserved. The requested reference target is explicitly reported as unmet.

```sh
python -S -m experiments.marginal_signed_atoms \
  --verify results/marginal_h6/signed_atom_active_restart/certificate.json
```

These constructions establish that independently signed local squares can leave the old balanced-clique family and improve an actual positive complement proof. They do not establish that the full flat support-four cone reaches the target, nor that the numerical master LP has reached its optimum. The full reference complement still contains368 explicit configurations and the physical ground-energy ledger is unchanged.

## Next unresolved work

The measured bottleneck is the master LP after new atom constraints are added, especially on the200-state block. An active-direction restart saves atom rows but still allocates19,900 pair variables. A justified next implementation can activate a zero-Hamiltonian pair only when a selected atom first fills it: all Hamiltonian nonzero pairs must be present initially, and every pair in every selected atom must be active. Other residual entries are then identically zero. New pricing can extend the current dual by B_ij=0 on absent pairs, activate all pairs of each chosen atom, and solve the enlarged exact model. Omitting filled pairs would be invalid; the final full rational DD replay must remain unchanged.

Other concrete options are persisting the complete discovered dictionary and native LP basis, or testing an equivalent unscaled-y coordinate system to reduce the large1/m_i coefficients. These are unimplemented optimization proposals, not established speedups. Independent local amplitudes and a broader support family remain available if the full flat family proves insufficient. No exact obstruction to that full family has been obtained here.

The final regression suite passed **250 tests in391.462 seconds**. All four new independent standard-library replays match their saved fields. These are restricted-family obstructions, a separator, and weaker complement bounds; the accepted energy/witness ledger remains114.
