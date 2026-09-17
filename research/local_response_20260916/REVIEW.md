# Review of the local exact certificate pass

## Findings

The local CAR construction is internally consistent on inspection.  The rung
Hamiltonian uses the two spin legs with the expected Jordan-Wigner signs; the
two-rung patch includes both directions on both longitudinal legs; half-rung
weights telescope across the three overlapping patches; and the boundary
messages cancel as `+J,-J` on adjacent patches.  The chemical shift
`-4(N-8)` is zero in the declared `(N_up,N_down)=(4,4)` target sector.

The PSD primal uses Clarabel's standard symmetric-vector scaling: diagonal
entries are unscaled and off-diagonal entries carry `sqrt(2)`.  The exact
replay, rather than the floating solver, decides acceptance, so floating
eigenvalues only choose candidate shifts.

The dual replay correctly requires PSD local density blocks, unit total trace,
and equality of adjacent one-rung partial traces.  Because the boundary
messages telescope against those equal marginals, its objective uses the bare
local blocks.  This is the expectation of the original four-rung Hamiltonian
in a compatible local-marginal assignment; it is not a physical global-state
claim.

## Important interpretation

`verify_dual` accepts one feasible local marginal assignment and therefore gives
a valid *upper ceiling* on any lower-shift sum: every PSD certificate must obey
`sum(shifts) <= objective(rho)` for that assignment.  It does not compute the
minimum over all compatible marginals, so `family_lower_ceiling_over_t` should
not be called the optimal ceiling unless an outer minimization is also solved.
The direction is nevertheless correct for proving that a particular accepted
primal cannot exceed the chosen dual ceiling.

The ceiling is relative to the fixed mathematical four-rung Hubbard model and
the fixed local PSD/message family.  It is not a physical ground-energy upper
bound, and it cannot by itself establish a minimum certified interval against
an unrelated or differently normalized upper endpoint.

## Conditions worth retaining in the report

State explicitly that all local charge sectors are included, while the global
target sector is fixed only through the chemical-potential identity.  The local
Fock enumeration is finite and local; `global_determinants_enumerated=0` is
accurate but does not mean the calculation is free of local basis enumeration.

The exact-family obstruction is conditional on the declared three overlapping
patches, two-rung supports, spin-number-preserving real symmetric boundary
messages, and the fixed chemical potential.  It says nothing about larger
patches, indefinite/cross-patch terms, adaptive projectors, or a different
response family.

I found no sign, transpose, or cone-direction error in the inspected code.
Numerical solver status and any claimed quantitative H8 conclusion still
require the focused test/replay receipts; source inspection alone is not an
independent computational verification.
