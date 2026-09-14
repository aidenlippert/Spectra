# Wave 2 information route: molecular H4 functional envelope

The assigned original active-space H4 fixture was loaded directly from `results/certificate_scaling/active_space_ladder/h4/fixture.json` (SHA-256 `8e64505deb4cc06a87f75e669b0f3073d15ebca72daa6d086e06eeea0b0f5120`). Exact Fraction matrix construction gives diagonal minimum `-3525719766889/10^12 = -3.525719766889`. Applying full-row Gershgorin gives the rigorous lower bound `-597155777003/125000000000 = -4.777246216024`, via `H_ii - sum_{j!=i}|H_ij| <= lambda_min(H)`. This is a fresh baseline calculation from the original H4 input; no legacy H4 certificate is used.

I tested the proposed 1-RDM-plus-pair-channel envelope on the molecular Hamiltonian conceptually at the level of valid operator inequalities. For each hopping matrix element, positivity of the 1-RDM and its hole matrix gives
`|gamma_ij| <= min(sqrt(gamma_ii gamma_jj), sqrt((1-gamma_ii)(1-gamma_jj)))`.
Adding diagonal pair variables does not sharpen this unless D/Q/G (or equivalent) overlap constraints couple them to the hopping channels. Therefore this envelope cannot honestly claim an H4 lower endpoint beyond the already replayed certificate. The prior ring experiment is not used as the molecular result.

This Gershgorin result is a valid molecular lower bound but is not a new functional-envelope improvement. The proposed 1-RDM plus pair relaxation remains **numerical-relaxation-only / deferred**, since no rational H4 dual was derived.

## Novel conditional proposal

Add a small H4 2-RDM block containing the hopping and pair channels, enforce D/Q/G principal-minor positivity, and export a rational dual after interval rounding. A decisive test is whether the rounded dual improves `-4.509041922449` while preserving the existing 70-state Hamiltonian identity. Failure at the same block size would falsify the expected benefit of pair channels without higher-order consistency.

## Scope and limitations

No ground-state oracle was used to construct a lower endpoint. The source H4 certificate's upper witness and lower replay were read only for comparison. Full 1-RDM data do not define the interaction functional by themselves; diagonal occupation data are strictly weaker. Any future recovery or transport route must retain fermion number/parity and provide an independently checked dual.
