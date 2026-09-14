# Small spin-resolved subspaces on the fixed ten-pattern H6 model

Keep the frozen H6 Hamiltonian, ten spatial density patterns, and exact
collective tail from the preceding molecular pass. A fresh quadratic
fermionic baseline is the common proof backbone. The preceding spin-summed
joint constraints have an exact accuracy obstruction and are not carried as
costly extra blocks. The Hamiltonian representation is unchanged.

Construct all a_i and Q_k,spin a_i generators from the ten patterns. Partition
them by the exact molecular spin-charge and parity symmetries. These four
frames each have 63 generators. No wavefunction, prior proof factor, or old
dual proposal is used for discovery. The old obstruction is a comparison
certificate only.

The adaptive arm gets 180 seconds end to end, including imports, initial
construction, the baseline solve, all pricing, candidate solves (including
rejected candidates), export, exact replay, and writes. It has at most four
enrichment rounds and 32 retained combinations. At each round, propose one
or two negative, numerically independent directions per symmetry class from
the current dual. Test both bundles, adding up to four or eight directions.
Directions are rationalized before any candidate solve. Within each class,
the selected combinations share a full anticommutator Gram block.

Choose a bundle by positive exact lower-bound gain divided by its complete
trial cost plus the round's pricing cost. Require at least 10^-6 Ha gain.
All tested valid lower certificates remain in the evidence ledger; distinguish
the best certificate found from the gain-per-cost continuation path. Stop
when the round/direction budget or time reserve is exhausted, or neither
bundle gives an accepted gain. The main solver cap is 30 seconds per trial,
reduced as necessary to reserve 20 seconds for export/replay. Numerical
moments guide discovery; they are not accepted dual witnesses.

Two separate controls each get 90 seconds end to end. One uses all four
63-generator spin frames. The other uses the best discovered small subspace
but removes off-diagonal Gram entries between its combinations. The latter
reuses the selected descriptors, so their full discovery cost remains charged
as ancestry. Neither control is a fresh cheap-discovery claim. A watchdog
records timeouts and preserves any earlier completed certificates; it never
counts an unfinished solve as a passing result.

Every certificate uses rational factors, exact adjoint pairing, the unchanged
CAR SOS verifier, and the unchanged collective tail verifier. No surviving
degree-six term can be omitted. Report compact and expanded proof size, map
work, all timing stages, and fresh standard-library replay. Compare energies
with the separately replayed existing rational upper witness; its original
exponential wavefunction discovery is outside the lower-bound method. The
success target is a 1.6 mHa certified interval. A failed solve or missed target
does not establish insufficiency of the full spin-resolved cone.
