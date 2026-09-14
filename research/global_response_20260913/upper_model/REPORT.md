# Compact Slater upper controls

The accepted implementation is now `research/global_response_20260913/slater.py`.
It contracts the exact rational occupied-orbital projector using Wick's theorem,
with the norm computed from the occupied columns. No many-body state vector is
created during this acceptance path. The earlier support-expanded ansatz in
`pair_rotation.py` remains a discovery diagnostic, not the final compact checker.

Four hydrogen fixtures have independently replayed rational upper certificates.
Their exact energies and timings are in `compact_replay.json` and the campaign's
`integrated.json`. On H6, the optimized Slater improves HF by 1.080400 mHa;
on the 1.6 and 1.73 Angstrom fixtures it improves HF by about 70.315 and
148.390 mHa. H8 stays at HF to numerical precision. The numerical H proposer
records precision-loss termination and has redundant inactive parameters;
none of these candidates is certified to be an optimum.

CH2 uses common zero-rotation orbitals with nested alpha/beta occupations.
These are exactly pure S=0 and S=1 controls, not optimized molecular results.
They do not reproduce the prior tight model gap certificate. Solver bounds,
model error and experimental accuracy remain distinct.

The final tests compare exact nonzero-rotation H4 energy and norm against an
independent explicit determinant expansion and compiled CAR energy evaluator.
They also check a doubly occupied rotated pair, nonnormal CAR contractions,
and refusal of an invalid pure-spin claim. The expanded test oracle is not
called by the accepted compact contraction.

The main campaign REPORT.md and integrated.json supersede earlier numerical
summaries. The separate lifted_upper branch is the actual response-lift control.
