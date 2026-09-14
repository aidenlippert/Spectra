> Superseded exploratory report. Do not use its conclusions as accepted evidence. See [the final campaign report](/Users/aidenlippert/Documents/Spectra/research/correlated_pair_20260913/REPORT.md) and [rejected prototypes](/Users/aidenlippert/Documents/Spectra/research/correlated_pair_20260913/REJECTED_PROTOTYPES.md).

# Direct self-consistent guide factors

This pass implements a Hamiltonian-input-only particle-hole guide and a paired
odd linear/mixed-cubic construction. For each mode the paired polynomial is
`B_i = a_i + a†_{i+1} a_{i+2} a_i`; its explicitly paired adjoint gives the
anticommutator closure. The exact CAR map confirms cancellation of all
degree-six terms, and the production verifier confirms residual degree at most
four. No moments, FCI vector, full Gram solution, or many-body basis enters
acceptance.

The builder passed exact replay on H4 and H6. The bounds are deliberately weak:
H4 lower `-31.025777711746 Ha`, H6 lower `-74.341206001716 Ha`, with residual
L1 equal to the magnitude of the lower because the guide contains no number
multiplier. This is a useful implementation and refusal baseline, not a useful
chemical certificate. The paired cubic construction was then run independently:
H4 lower `-31.025777711778 Ha`, H6 lower `-74.341206001764 Ha`; it passes exact
degree-six cancellation but does not improve the L1 residual. This isolates the
remaining issue as coefficient matching/self-consistent weighting, rather than
an absent cubic cancellation identity.

All numerical calls were run through `budget.py` with `sc_` names. The builder
is exposed through `guide_certificate(..., correlated_moments=...)`, allowing a
future correlated-state proposal to rank factors without changing the exact
acceptance path.

The requested creator-distribution map is implemented as
`tau_from_perturbation(V, weights)`. It assigns each canonical monomial only to
its creator labels, divides by the exact positive sum of guide weights, applies
the equal-degree half convention, and checks the cross identity through
`tau_cross`. A bounded `self_consistent_update` performs up to eight damped
updates with the full induced anticommutator retained. It refuses nonpositive
weights and unsupported constant/annihilator-only terms; it never clips them.
The exact cross-map and one-step update tests pass. Multi-step updates expose a
real algebraic issue: the induced anticommutator creates lower-degree terms
outside the creator-distribution domain, so the current restricted map refuses
iteration rather than silently discarding them. This is the concrete next
closure requirement for a complete self-consistent implementation.
