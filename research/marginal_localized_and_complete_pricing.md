# Localized orbitals, complete support pricing, and a rank obstruction

The subsequent [valence-reference result](marginal_valence_reference.md) resolves the finite gap and transfer requirements by changing P to the entire20-state singly occupied manifold. The fixed-P cone and norm obstructions below remain valid.

The localized H6 search now supplies an exact complement lower bound of **−6.347558453269706 Ha**, compared with its exact empty diagonally dominant seed of **−6.412577307516451 Ha**. Its exploratory threshold is **−6.323 Ha**, leaving **0.024558453269706 Ha** uncertified. This is a new Hamiltonian representation and freshly selected retained space; its complement differs from the original-basis complement. It is not a new ground-energy interval.

The original-basis continued certificate reaches **−6.519147644269328 Ha**; two complete-pricing rounds further improve it to **−6.516166792109046 Ha**, leaving0.251267691069046 Ha to the original requested threshold. A new exact dual bounds every factor-width-three certificate on that same complement above by **−6.5334123152169274 Ha**. Thus the achieved four-coordinate proof is strictly outside the entire three-coordinate family by an exact separation of approximately0.0172455231078814 Ha. The certificate Hamiltonians and retained spans are checked to match in `factor_width_strict_separation.json`. The stronger original target **−6.26489910104 Ha** remains unmet.

An exhaustive numerical scan has now checked all **66,018,250** two-, three-, and four-coordinate supports of an actual 200-dimensional search dual. It found 37,120,439 supports with minimum eigenvalue below −1e−8, including 36,689,554 quadruples. This directly identifies substantial omitted constraints in that search state. It does not prove either feasibility or infeasibility of the complete factor-width-four cone.

## What changed under localization

The target is the Löwdin symmetric orthogonalization of the six atomic orbitals. Regenerated RHF integrals are aligned to the saved rational Hamiltonian using the exact original two-electron index convention and exhaustive spatial phase signs. The best signs are `[1,-1,1,-1,1,-1]`. The coefficient L1 mismatch is 1.8855206146052513e−10 Ha; the largest individual mismatch is 4.924441930620471e−13 Ha. These are numerical provenance checks, not integral-error certificates.

For the aligned canonical orbitals, the convention is

\[
U=C_{\rm MO}^{T}S^{1/2},\qquad a_{\rm old}=Ua_{\rm new},\qquad h_{\rm new}=U^Th_{\rm old}U.
\]

The corresponding localized atomic-orbital coefficient matrix is S^(−1/2). Fifteen rational Givens rotations approximate this target, with maximum entry error 4.903811074985498e−7. Their assembled entries have at most 327 numerator bits and 328 denominator bits. The rational matrix is exactly orthogonal. Independent standard-library replay checks CAR, every transformed Hamiltonian coefficient, exact inverse recovery, and spin symmetry. No coefficient truncation changes the Hamiltonian during this transformation.

The original and rotated 400-state numerical spectra agree within 2.7533531010703882e−14 Ha. Both reference selections use the same minimum-diagonal seed policy and fresh residual-driven determinant selection, capped at 32. Here P is explicitly the span of those determinant basis vectors; deleting their labels gives exactly Q=I−P, not a complement of a single Ritz vector.

| Numerical diagnostic | Canonical basis | Localized basis |
|---|---:|---:|
| Selected seed determinant | 63 | 1386 |
| Q block dimensions | 200 + 168 | 368 |
| Ground energy | −6.3330586262329955 | −6.333058626233006 |
| Physical Q floor | −6.26479910103936 | −6.311577514960223 |
| Physical Q gap above ground | 0.06825952519363554 | 0.021481111272783515 |
| FW2 comparison ceiling | −7.985247211264101 | −6.412543371350798 |
| FW3 comparison ceiling | −6.533412315216928 | −5.957259586030358 |
| FW4 comparison ceiling | −6.127741516801051 | −5.8500976359225385 |

Localization improves the comparison geometry but reduces the physical gap available to this P32. A localized default-seed control starting at determinant63 gives a numerical Q floor −6.306887745920947 and a FW2 comparison ceiling −6.4074459983582805. Seed choice matters, and the canonical P labels cannot simply be reused as though they represented the same physical subspace.

Artifacts: `results/marginal_h6/localized_basis`, `localized_orbital_target_reproduced.json`; implementing modules: `experiments/marginal_lowdin_target.py`, `marginal_localized_basis.py`.

## Exact comparison-family duals

For a real symmetric physical Q block H and a positive rational vector v, define

\[
X_{ii}=v_i^2,\qquad
X_{ij}=-\operatorname{sgn}(H_{ij})\frac{v_i v_j}{k-1}\quad(i\ne j).
\]

Every principal block of X of order at most k is PSD: diagonal congruence removes v, leaving diagonal1 and absolute off-diagonal row sum at most1. Padding the witness by zero outside its physical support preserves this property. Therefore, if H−γI has a factor-width-k decomposition,

\[
\gamma\le\frac{\langle H,X\rangle}{\operatorname{tr}X}
=\frac{\sum_i H_{ii}v_i^2-\frac{2}{k-1}\sum_{i<j}|H_{ij}|v_i v_j}
{\sum_i v_i^2}.
\]

Numerical comparison eigenvectors propose integer v; the quotient is then reconstructed exactly from Hamiltonian actions on validated Q determinants. Production replay uses the analytical PSD proof, not enumeration of principal minors.

| Physical certificate | Exact ceiling, rounded to float | Target | Excludes target? |
|---|---:|---:|---|
| Original basis, k=3 | −6.5334123152169274 | −6.26489910104 | Yes |
| Localized basis, k=2 | −6.412543371350795 | −6.323 | Yes |
| Localized basis, k=3 | −5.957259586030356 | −6.323 | No |

All three independent replays match the saved receipts. These ceilings are upper bounds on achievable proof-family thresholds, not physical spectral lower bounds. Failure to exclude does not establish feasibility. The original three-coordinate ceiling and the stronger achieved four-coordinate lower certificate give an explicit finite H6 separation between these cones at the achieved threshold.

The earlier files under `comparison_family_duals` are invalidated as physical certificates: their replay trusted a supplied matrix without reconstructing it from H/P, and their stored exclusion flag had the wrong comparison sign. Accepted replacements are exclusively under `comparison_family_duals_checked`.

## Localized positive certificates

The seed constructor uses complete physical Q coverage, positive integer weights proposed from the comparison ground vector, and balanced three/four-coordinate directions drawn from each row's six strongest couplings. It does not enumerate every initial clique. The initial dictionary has6,421 directions. Its empty exact DD proof has threshold −6.412577307516451; the requested threshold remains separate from this achieved baseline.

| Localized export | Exact Q lower, Ha | Explicit atoms | Certificate bytes |
|---|---:|---:|---:|
| Empty weighted DD seed | −6.412577307516451 | 0 | 1,451,608 |
| Eight adaptive rounds | −6.354266190861379 | 3,729 | 2,307,087 |
| Eight-round scaled residual | −6.354184984289176 | 3,729 | 2,315,310 |
| Thirty-two further rounds | −6.347558453269706 | 7,265 | 3,237,083 |
| Resumed scaled residual | −6.347558453269706 | 7,265 | 3,241,169 |

The first search takes29.47 seconds in its proposal routine; the continuation takes165.69 seconds. These are observations on different LPs and dictionaries, not controlled speedup claims. The later scaled form adds bytes without improving the exact endpoint, so the ordinary resumed certificate is preferred. Significant negative sampled pricing values remain.

The connected 368-dimensional block required an explicit bounded extension of the atom verifier's block limit from256 to384. The old default remains256 for other callers. Physical-state validity, complete coverage, uniqueness, and absence of interblock coupling remain mandatory. Tests cover those gates at the larger size.

## A basis-independent obstruction to reducing the hopping norm

Let

\[
\Delta=t\sum_{\sigma=\alpha,\beta}(a^\dagger_{p\sigma}a_{q\sigma}+a^\dagger_{q\sigma}a_{p\sigma})
\]

act in the whole fixed-(Nα,Nβ) sector of m spatial orbitals, with t≠0 and 1≤Nσ≤m−1. Diagonalizing this orbital pair gives single-spin eigenvalues ±t. Occupying the appropriate extremal orbital, leaving the other empty, and choosing the remaining Nσ−1 electrons among m−2 spectator orbitals gives an extremal eigenspace of dimension

\[
D=\binom{m-2}{N_\alpha-1}\binom{m-2}{N_\beta-1}.
\]

If an orthogonal retained projector P has rank r<D, dimension counting gives a nonzero intersection between that eigenspace and ker P. On this intersection QΔQ retains an eigenvalue of magnitude2|t|. Compression cannot increase operator norm, so

\[
\boxed{\|Q\Delta Q\|=2|t|\quad\text{for every rank-}r\ P\text{ with }r<D.}
\]

For H6, D=binom(4,2)^2=36. Consequently **every rank32 retained projector leaves the hopping norm exactly0.04 Ha at t=1/50**, independent of orbital basis or the method used to select P. Raising the retained rank to at least36 is necessary to make a strict norm reduction possible; it is not a sufficiency theorem.

Combined only numerically with the localized physical Q floor, the best scalar-shift threshold would be at most about −6.351577514960223 Ha, below the existing connected ground interval near −6.3332381. This combined comparison is a numerical diagnostic because that Q floor has not been exported as an exact Q Rayleigh witness here. The norm/rank theorem itself is exact. Thus localization's smaller proof error must not be mistaken for repairing the original connected scalar-shift pipeline. A direction-sensitive perturbation inequality, a different/larger P, or a different reference mechanism is still needed.

## Complete support scan and remaining limits

`marginal_complete_pricing.scan_supports` streams every principal support in batches, retains bounded negative eigendirections by order, and records actual evaluated counts and near-boundary cases. The final partial batch is explicitly included; a regression catches its former omission. The selection step keeps only each batch's best candidates, since any discarded candidate already has enough better candidates in its own batch to be excluded globally. This preserves complete coverage while avoiding Python work on millions of retained-candidate objects.

The actual original-basis LP is integrated through `--support-pricing all_supports`. Each negative numerical direction is rounded and checked again before it becomes an LP constraint; final positive decomposition acceptance remains exact. Even a scan with no negative eigenvalue would be a numerical observation, not an exact dual-PSD certificate. A bounded top-candidate list also need not contain every novel rounded direction.

For the first200-state scan, counts are19,900 pairs,1,313,400 triples, and64,684,950 quadruples. The worst observed eigenvalue is−0.8432743000700084;161,528 supports lie within1e−8 of zero. Both native blocks completed two rounds. The168-state scan checks32,809,154 supports per round; its first round finds21,540,499 negative supports. The final exact block bounds are−6.516166792109046 and−6.353908051232628 Ha. The common native certificate is2,439,025 bytes. Full enumeration resolves missing finite support coverage, but its combinatorial cost is not a scalable representation of general matter. The measured run used the pre-optimization heap loop; no measured speedup is claimed for the later per-batch candidate optimization. Complete scans can overrun the between-round time budget; the first200-state two-round block used437.63 seconds against a360-second soft budget.

The accepted ground-energy/witness ledger remains114. Every current molecular positive proof still references400 spin-sector determinants, and the localized Hamiltonian has1,818 rational terms instead of918. The open problem remains constructing a small, useful reference and positivity certificate with controlled error, state count, and rational bit growth as molecular size increases. General marginal representability, finite-temperature response, kinetics, and synthesis do not follow from these finite ground-sector certificates.

The full295-test suite passes in362.213 seconds. Independent standard-library replays cover the rotation, physical duals, native positive certificates, and localized positive certificates; receipt matching is recorded in `results/marginal_final_validation.json`.
