# Sparse wedge multiplier repair audit

The root replacement uses a sparse LP with absolute-value off-diagonal
variables and one Gershgorin lower endpoint per exterior-power body. Its
objective includes the exact constant residual and binomial(N,k) lifting.
Both old and proposed certificates undergo exact replay; a weaker proposal
is discarded. Factors, Hamiltonian, and scalar b remain fixed.

The strengthened M4,N=2 control perturbs ONLY the number-ideal multiplier by
Hermitian hopping of magnitude 1/10. Its initial wedge lower is negative;
the LP restores the exact lower bound zero with H, b, and every factor
unchanged. A separate unperturbed control remains exactly zero. This replaces
the earlier weaker test that perturbed the Hamiltonian instead.

The LP objective is a numerical proposal. Only exact residual replay
certifies its result. The old coefficient-L1 verifier may become much weaker
when the multiplier introduces large cancelling terms; the stronger reported
bound requires the named wedge replayer and must not be attributed to the
old verifier alone.
