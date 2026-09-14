# Wave 2 bounded audit

## Passing

- Wave 2 upper progress files support the corrected claims: H4 passes at requested 128, H6 at 512, and H8 at 2048. H8 at 1024 fails; the report correctly says so. Requested basis and nonzero rounded witness supports are separately recorded.
- Upper witnesses are exact replay artifacts against the original fixtures; no interval claim is made from projected Ritz residuals.
- Sparse no-prune H4 and H6 receipts are correctly treated as lower-only results. H4 no-prune is slightly better than the earlier 512-atom baseline; H6 is catastrophically weak. Numerically, the conservative H4 result `-7.067...` is a *stronger* lower bound than no-prune `-8.408...` (larger lower endpoint), though still far from useful accuracy. The width-4 run ends at `-4.907...`, which is stronger still but remains weak and is not a matched fixed-budget comparison.

## Provenance and interpretation findings

- Sparse dual receipts do not contain fixture path or fixture SHA-256 fields. The associated report names H4/H6 hashes, but machine-readable provenance should be added before treating remote reruns as independently identity-bound.
- The Wave 2 upper report's `1.77 s` H8 total is a run-local observation including exact replay, not a full campaign cost. It correctly avoids a scaling claim, but diagonal reference scanning and all setup costs should remain visible in any cross-method timing table.
- The no-prune comparison is a matched atom-cap comparison only where the same fixture, atom family, and seed construction are verified. The width-4 history begins at 512 atoms and ends at 1792 after 20 rounds (1280 additions), so the count is internally possible; it is not a 2048-retained result. Report `atoms_retained`, `rounds_solved`, and seed counts explicitly.
- “Full quadratic” terminology is unsafe for the 512 atom result: receipts identify a retained positive quadratic atom family with `full_Gram_constructed=false`, so it is an inner-cone LP, not the full PSD/SOS cone.
- The H4 conservative lower `-7.067...` is weaker than H4 no-prune `-8.408...` because lower bounds are ordered numerically; less negative is stronger. Any prose calling this “degraded” should state that pruning degraded the certificate relative to no-prune while both remain far from useful accuracy.

No blocker was found in the upper endpoint arithmetic or the corrected H8 pass/fail interpretation. Remaining blocker for strong sparse-dual provenance is missing fixture identity in the receipts.

## Final code-level checks

- `complement/selected_schur.py` uses exact CAR reconstruction and exact rational LDL/PSD checks. Its norm-radicand ceiling is conservative (the exact rational post-check increments any under-ceiling), and the selected-subspace outputs are valid lower bounds under the stated block bound. The diagonal-only value is explicitly diagnostic.
- `theory/parent_test.json` is correctly interpreted as a global HF projector diagnostic. Fractions confirm the eta-zero width is `262389/250000 = 1.049556`, not `0.262389`; eta-one gives `141281/1000000`. These are weak but valid exact relative bounds.
- Homotopy checks use the corrected H4 hash and exact LDL zero-pivot regression. Their Gershgorin rows use a determinant trial upper, so they are valid finite-H intervals but not continuation/scaling results.
- Tensor-network H4 reports an exact replayed witness with optimizer failure (`maximum iterations exceeded`); it must not be called converged DMRG. H6 is correctly recorded as incomplete.
