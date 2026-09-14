# Direct dual obstruction for the complete quadratic square cone

This is a family-specific obstruction experiment for the eight-mode, four-
particle H4-square Hamiltonian.  The cone contains every singleton square and
both signed two-word square polynomials from `dictionaries(8, "quadratic")`:
136 words, 18,496 square atoms, and 5,093 canonical CAR coordinates after
including the degree-six number-ideal products.

The sparse HiGHS LP minimized the Hamiltonian pairing over dual vectors `y`
with `y(1)=1`, exact annihilation of every `(Nhat-N)X` basis element, nonnegative
pairing with every square atom, and `||y||_infty <= 1`.  It solved in about 24
seconds.  The numerical objective was approximately `-4.4590073065`.

The rationalized dual was checked against every constraint, rather than
reported as a numerical certificate:

- constant constraint: exactly `1`;
- maximum absolute ideal pairing: exactly `0`;
- minimum square-atom pairing: exactly `0`;
- infinity norm: exactly `1`.

Thus this is a valid exact dual obstruction for this *specific complete
quadratic diagonal/coefficient-l1 norm cone*, assuming the supplied exact
Hamiltonian encoding.  It is not an obstruction to cubic, quartic, local,
nonquadratic, or general SOS certificates.  No FCI or many-body sector matrix
was used.

The full rational dual vector and all exact checks are in
`results/certificate_scaling/direct_dual_obstruction/receipt.json`.
