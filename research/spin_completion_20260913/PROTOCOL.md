# Spin-changing joint operators on the fixed ten-pattern H6 model

The question is whether richer joint operators from the same ten spatial
patterns remove the accepted full-spin-frame counterexample and improve a
compact certified energy bound. The frozen Hamiltonian, collective tail,
quadratic baseline, body-one number ideal and coefficient-L1 residual rule
remain the same. No wavefunction or prior SOS factor is used for discovery.

First complete the density spin indices: Q_k^(s,t) = sum_pq L_kpq
a^dagger_(p,s) a_(q,t), for all four (s,t). Form a_i and Q_k^(s,t) a_i,
partitioned by exact spin charge and molecular parity. Reuse the existing
exact full-spin dual as the diagnostic input, with its discovery ancestry
reported. Replay that dual and find negative directions in the completed
frame. Rationalize and exactly replay each selected separator; a negative
floating eigenvalue alone establishes nothing. Diagnostic budget: 180 seconds
including imports, construction, old dual replay, pricing and export.

Conditional on exact separation, run a fresh adaptive energy search under a
240-second outer watchdog, at most four rounds and 64 combinations. The first
round uses the exact counterexample's directions; subsequent rounds price
the current energy dual. Test one and two negative independent directions
per symmetry block, couple the selected directions, and choose continuation
by exact bound gain per measured trial plus pricing cost. All candidate work
counts. Main solve cap: 40 seconds, with 25 seconds reserved for export/replay.
Require at least 10^-6 Ha accepted gain. Save the best result independently
of the continuation path. Discovery of the original counterexample is shared
ancestry, not free end-to-end discovery.

Controls: complete spin-changing frame under 240 seconds, solver cap 150
seconds; the selected small directions with diagonal Gram weights only under
90 seconds, solver cap 50 seconds. Preserve failures and do not retry by
silently changing the budgets. Refuse automatic numerical basis truncation.
Any exact dependence removal must preserve and verify the full operator span.

All new lower certificates reconstruct exact CAR polynomials, tie each
operator to its exact adjoint, and pass the existing rational SOS and tail
verifiers. No surviving sixth-degree terms may be omitted. Recheck accepted
intervals in a fresh standard-library-only process against the rational
reference upper. Report proof size, coefficient counts, complete costs and
the 1.6 mHa target. Failure to reach it is not a proof of cone insufficiency.

If spin completion does not separate or does not support a sufficient bound,
use the evidence to choose the next operator family. Further families and
exact obstructions require a separately recorded bounded protocol. Keep all
earlier sources and artifacts unchanged.
