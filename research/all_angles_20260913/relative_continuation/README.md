# Relative perturbation and certificate continuation

This route tests the bound

`H = Eref I + P + R`, `P >= 0`, and `R + eta P + eps I >= 0`

which implies `H >= (Eref-eps)I` for `0 <= eta <= 1`.  The relative term can
be much sharper than charging `||R||` against the whole operator.  The test
below uses exact rational matrices and an independent dense eigensolve.

The construction is deliberately finite dimensional and uses abstract energy
units (no Hartree interpretation). It is evidence for this certificate route,
not a chemistry scalability theorem.
