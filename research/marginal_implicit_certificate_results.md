# Both energy endpoints in implicit coordinates

The defect-Dicke representation now supports complete finite-model energy intervals: lower positivity proof, variational upper witness, numerical direction proposal, and targeted second enrichment. None of these paths needs a Fock configuration list. Targeted constructions with short Hamiltonian-action upper witnesses give widths **9.976325470374013e-8** for the cycle and **9.562520469916092e-8** for the mixed interaction, using **20 and 42 retained dimensions** in the lower positivity checks.

This extends [the implicit lower-bound construction](marginal_defect_dicke_results.md). All current complete-interval experiments remain M=10,N=5 at perturbation strength 1/100 around the matched reference. They do not establish general chemistry or large-system accuracy.

## Exact upper witness and stable proposal

For retained implicit columns U, compute G=U^TU and A=U^THU using exact combinatorial overlaps. Any nonzero integer vector x gives a physical state Ux and the exact variational upper bound

\[
E_0\le\frac{x^TAx}{x^TGx}.
\]

The verifier requires a strictly positive exact denominator. It reconstructs the basis from the Hamiltonian and the stored targeting recipe; it does not trust a supplied Gram matrix or a numerical energy.

Direct floating calculations with G can lose definiteness. The proposer instead factors G=L D L^T exactly and constructs the rational transform T=L^-T. It verifies T^T G T=D, then solves the floating problem in the normalized orthogonal coordinates

\[
C=D^{-1/2}T^TATD^{-1/2}.
\]

The eigenvector only proposes a state. Rational orthogonal coefficients are mapped through T exactly, and the resulting Rayleigh quotient is recomputed exactly. A test uses a metric containing entries of order 10^40 whose unit correction is lost in floating representation; the exact transform still recovers the intended two-dimensional problem.

### Compressing the witness coordinates

Exact coordinate conversion initially produced 4,061-bit integers for the mixed model. The revised proposer rounds the original-basis coordinates to increasing significant-digit precision and accepts that rounding only when

\[
\frac{(x-\widehat x)^T G(x-\widehat x)}{x^TGx}\le10^{-24}
\]

holds in rational arithmetic. This controls physical-state error despite the ill-conditioned coordinate system. On the mixed first space, 18 significant digits give targeting integers of at most **71 bits**, while retaining the displayed Rayleigh value 3.2805844370809547. The error certificate is stored separately from the final energy certificate. The latter remains valid for any nonzero proposed integer vector, regardless of proposal quality.

## Targeted second enrichment

For a proposed integer x, form the exact leakage

\[
w=HUx-U G^{-1}Ax.
\]

The checker verifies U^Tw=0. It retains the exact H0-Krylov closure of w, verifies that closure is orthogonal to U, and then rebuilds the enlarged Gram, Hamiltonian, and leakage matrices. The lower endpoint is accepted only after exact Schur positivity; the upper endpoint is another exact Rayleigh quotient in the resulting space.

The cycle grows **6→14→20**. The mixed targeting recipe closes at **6→36→42**. Its grouped construction completed in 41.68 seconds and independently replayed in 41.13 seconds, before the upper-only refinement.

## Verified intervals

| Fixture / construction | Retained dimension | Certified width |
|---|---:|---:|
| Cycle, first enlargement | 14 | 2.8519887820773048e-5 |
| Mixed, first enlargement | 36 | 2.8906380954868705e-5 |
| Cycle, targeted second enlargement | 20 | 1.3046477061915772e-7 |
| Mixed, targeted second enlargement | 42 | 1.1119102727038293e-7 |
| Cycle, targeted plus two upper-witness steps | 20 | **9.976325470374013e-8** |
| Mixed, targeted plus one upper-witness step | 42 | **9.562520469916092e-8** |

These first-space intervals are wider than the earlier intervals using an explicit upper witness, because the new upper states are confined to the retained space. The targeted cycle improves its upper state as well as its lower proof. Its interval is approximately

\[
3.2796999517\le E_0\le3.2797000821647706.
\]

Exact endpoints are in the certificate. Energies use the model's interaction units.

The refined upper witness is a rational vector in span{psi,H psi,...,H^k psi}, with k=2 for the cycle and k=1 for the mixed model. Only the original retained coefficients and two or three extra scalar coefficients are exported. Its Gram and Hamiltonian matrices are rebuilt with the implicit algebra; the lower proof is rechecked unchanged. Independent tests compare these moments with explicit Fock actions through H^3. The strongest intervals have rounded endpoints 3.2796999517 to 3.2797000514632546 and 3.2805774206 to 3.280577516225205.

Two proposal failures were resolved without weakening proof gates. An unnormalized witness metric overflowed during floating conversion; an exact common scaling of all projected matrices fixes that. A very large rounding-error fraction then exceeded Python's decimal-string limit; the comparison remains exact, while diagnostics now record the rigorous 1e-24 error bound instead of that enormous fraction. Both failed directories contain diagnostic receipts and no accepted certificate.

The cycle certificate was subsequently regenerated with the checked coordinate compression enabled throughout. Its initial upper coordinates use at most 102 bits, and the final certificate shrank from 24,937 to 9,355 bytes. Its independently replayed width is 9.976325470380547e-8, differing from the earlier version only at the shown trailing digits. Both variants are preserved.

## Eliminating zero overlap work

Atoms with different empty/double pair patterns belong to mutually orthogonal pair-occupation sectors. Inner products can group by that pattern and skip cross-sector contractions exactly. This bookkeeping does not assert that the perturbed Hamiltonian preserves those sectors; the Hamiltonian action still connects them.

All rational pivots and endpoints match before and after grouping. In separate observed replays, the cycle targeted case went from 66.44 to 4.60 seconds and the mixed first-space case from 60.90 to 10.08 seconds. These are individual runs, not an isolated hardware benchmark. Earlier receipts and the ungrouped search diagnostics are preserved. The ungrouped mixed targeted search was stopped at an observed 300-second wall budget with no accepted certificate. Its 42-dimensional closure and targeting proposal remain diagnostics; the stop is not an infeasibility result.

## Verification and next blockers

The complete marginal suite passed **136 tests in 43.131 seconds**. The seven new tests check the exact orthogonal transform, physically bounded coordinate compression, rejection of bad metrics and zero or malformed upper vectors, Hamiltonian tampering, independent equality with a Fock-space Rayleigh calculation, and reconstruction of the targeted space. Replay tests also replace the old explicit action and upper-witness helpers with functions that throw if called.

The next blockers are efficient higher-order enrichment, energy certificates beyond ten modes, discovery of useful reference spaces, and accuracy control outside the certified perturbative region. A separate [order-versus-error theorem](marginal_enrichment_error_bound.md) now bounds full enrichment by (B-b)[eta/(c-u)]^(2r+2) when eta<c-u. It does not guarantee the accuracy of selective targeting or efficient cost when r grows.

All nine accepted certificates passed independent standard-library replay and were matched to their declared fixture Hamiltonians. The partial ungrouped search and two failed proposal/serialization attempts are excluded from that count.

Implementation: `experiments/marginal_implicit_certificate.py`. Certificates and independent receipts: `results/marginal_implicit_certificate/`.

```sh
python3 -S -m experiments.marginal_implicit_certificate --verify results/marginal_implicit_certificate/cycle_1_100_targeted/certificate.json
```
