# Independent campaign audit

Audit date: 2026-09-13. This is a read-only review of artifacts currently present; it does not treat the 72-route brainstorm as 72 executed experiments.

## Passing claims

- The Schur toy certificate is valid where its independently certified complement bound has `mu > theta`; the refusal gate correctly prevents using residual size alone. The H4 attempt is an honest refusal, not a molecular result.
- Relative-continuation's two-dimensional rational PSD premise and resulting lower bound are valid for the abstract family. Its widths use the exact dense eigenvalue as a post-freeze oracle and are not certified molecular intervals.
- CISD witnesses are valid variational upper endpoints after exact rational replay. They are cheap witnesses, not evidence of scalable discovery: the projected CISD diagonalization has basis sizes 53/262/849 and is already quadratic in the selected subspace.
- Spin ablation correctly establishes a representational basis-count comparison and exact replay of the original H; it does not establish SU(2) invariance of the unrestricted operator family.

## Failed or misleading claims to correct

- `sparse_dual/sparse_chordal_probe.py` comments that treewidth is bounded by the largest block. This is false in general: local/maximal block order does not bound treewidth of the global interaction graph, and the reported 52 edges on H4 is only a measured graph statistic.
- Sparse dual REPORT/receipt replays the old certificate with lower `-4.459009...` (the older H4 fixture). The active ladder H4 is approximately `-3.667`; these are different Hamiltonian fixtures and must not be compared or presented as an active-ladder result. Any sparse-dual scaling claim is therefore unverified until rerun on active H4/H6/H8.
- `spin_weight2(w)==0` encodes a zero spin-weight / Sz-conserving word filter. It is not by itself proof that every operator is SU(2)-invariant; full invariance requires commutation with Sx, Sy, and Sz (or an equivalent scalar-coupling construction). Rename the route description accordingly.
- Relative-continuation widths are widths against an eigval oracle in abstract units. They are not verified trial-upper intervals and do not demonstrate molecular continuation, reaction-energy cancellation, or extensivity control.
- Schur's `oracle naive residual lower` labels in toy artifacts must remain explicitly uncertified whenever `mu<=theta`; a residual-norm endpoint is not a lower bound without a valid global operator inequality/complement certificate.

## Coverage accounting

Executed substantive families visible here: sparse H4 replay only; abstract relative bound; Schur toy plus H4 refusal; H4/H6 spin-basis/one matched solve; CISD H4/H6/H8 witnesses; spinless-ring information diagnostic; channel-rank metrics. Frontier and other route directories are design/proposal documents. No Lambda/GPU run is evidenced in these receipts. Existing route numbers are coverage labels, not executed experiments.

## Required follow-up gates

1. Rebuild sparse-dual candidate and replay against the active-ladder fixture, recording fixture hash/path, full discovery scan, graph treewidth (or a certified upper bound), and verifier result.
2. Replace SU(2) wording with Sz/number-conserving wording unless commutator checks are added.
3. Report every interval endpoint with provenance: certified lower, exact rational trial upper, or oracle-only diagnostic; never combine oracle lower/upper into a claimed certificate.

## Follow-up audit (current tree)

- Information refinement now passes the asymmetric coherence regression and uses the valid `min(sqrt(n_i n_j),sqrt((1-n_i)(1-n_j)))` bound. The ring result remains a diagnostic, not a molecular certificate.
- Exact H4 oracle passes its corrupted-diagonal and legacy-fixture hash refusals, but is explicitly an exponential oracle baseline.
- Sparse-dual is now correctly tied to active H4 and falsifies the 64-atom budget: lower `-14.02` versus paired upper `-3.5257`. The 512-atom result is a direct H-ranked quadratic LP with 512 retained atoms, not the full quadratic/SOS cone; its `full_Gram_constructed=false` receipt must remain prominent. It is also not a treewidth experiment.
- Schur code's `diagonal_mu` is metadata only; the endpoint uses full complement Gershgorin. A diagonal-only minimum is not a lower bound on an off-diagonal complement and must not be described as one. Focused tests pass when run from the Schur directory (`2 passed`); the repository-root module invocation fails due to import-path setup.
- Relative interacting-cycle code uses a dense rank-one vector with fixed per-coordinate amplitude. Its vector norm and rank-one spectral scale grow with dimension, and its fixed `eps` is a total-energy slack. Therefore the apparent width behavior cannot be promoted to physical size scaling without normalization/extensivity analysis.
- Selected-CI H6 receipt records a valid exact replay and width near `0.00106` against the current lower endpoint, but its dense projected matrix (`154449` entries) remains a bounded experiment; it does not establish scalable discovery.

## Selected-refinement implementation audit

The new signed-residual sparse-CI implementation passes both focused tests. Integer compilation agrees with the independent CAR action on sampled states, and every rounded witness is independently replayed against the original rational Hamiltonian. Each stored upper endpoint is therefore sound, subject to the explicit replay/refusal path.

The algorithm remains heuristic: candidate ranking uses signed residual accumulation and a denominator floor; each iteration builds a projected matrix and retains cached state actions. `full_fock_enumeration=false` is accurate, but gives no polynomial or scalable-discovery result. `connected_support_saturated` means only saturation of explored support. Projected residual norms are numerical diagnostics, not rigorous error bars. Requested budgets and post-rounding witness supports are correctly reported separately.

H4/H6 receipts provide fixture hashes, exact replay, matrix counts, candidate counts, and cache sizes. H8/H10 and GPU/Lambda outcomes remain unverified until raw receipts and witness hashes are present. GPU reporting should distinguish CPU-subset discovery from GPU-all-state enumeration.

## Final evidence correction

The preceding provisional wording is superseded. Fresh H4/H6/H8/H10 interval artifacts are locally present and all referenced files were rehashed successfully against their receipts: each upper reference, certificate, and (where present) spectral proof matched exactly. The remote-looking `/home/ubuntu/.../replay/` prefixes resolve to the corresponding local workspace paths.

H10's `273.108 s` is the complete two-sided proof replay time and excludes the `66.56 s` selected-CI discovery time and broader campaign cost. GPU evidence is dense SciPy `eigh` on the CPU-selected subset versus the full dense GPU eigensystem; it is not sparse selection and not GPU all-state enumeration.
