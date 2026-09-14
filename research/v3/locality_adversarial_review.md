# v3 locality certificate adversarial review

Audited `experiments/v3_locality.py`, `tests/test_v3_locality.py`, and executed `python -m unittest tests.test_v3_locality -v`: all 4 tests passed. Independently ran `run()`, obtaining representation sizes 30, 98, 109, 109, 109 for (n=3,6,12,24,48), with the uniform bound (0.0002389124397103113) in each case.

## Mathematical checks

The admitted source class is explicit: one- and two-site Pauli Hamiltonian terms with exact rational coefficients, bounded incident strength \(\kappa\), and local depolarization damping each nonidentity Pauli factor at rate \(\gamma\). For a Pauli column of weight (w), the commutator column sum is bounded by (2\kappa w), so \(\gamma>2\kappa\) gives \(\lambda=\gamma-2\kappa>0\). The implementation computes this threshold from the actual Hamiltonian and rejects failures.

The neighborhood is the exact commutator closure through depth (R). The short-time truncation is a conservative full-path bound: each omitted depth-(k) coefficient \(\ell_1\) mass is bounded by \((2\kappa)^k k!\), while the Dyson simplex contributes (t^k/k!\), leaving \((2\kappa t)^k\). Summing the geometric tail gives the implemented \((2\kappa t)^{R+1}/(1-2\kappa t)\). Since omitted paths are a subset of full paths, no extra factor of two is needed there.

For late times, both full and reduced semigroups have \(\ell_1\) norm at most (e^{-\lambda t}); their difference is bounded by (2e^{-\lambda t}). Taking the maximum of the early Dyson tail and late decay at the cutoff provides a uniform all-time bound. The code verifies the rational exponent enclosure and independently checks the structural certificate.

For any density matrix, each Pauli expectation has magnitude at most one, so an \(\ell_1\) coefficient error directly bounds the expectation error without Hilbert–Schmidt dimension factors. The reduced prediction additionally reports a separately bounded rational-uniformization error. The dense eigensolver tests are numerical cross-checks, not certificate construction.

## Residual scope limits

The theorem is conditional on known local depolarization, bounded degree/incident strength, exact admitted Hamiltonian terms, and the chosen one-site seed. It does not apply to arbitrary correlated noise, long-range/unbounded incident interactions, or a Hamiltonian supplied only approximately. The certificate is a compact local-observable guarantee, not a global state trace-distance guarantee.

The code selects a commutator neighborhood algorithmically, but `run()` explicitly reports `compounding_learning_demonstrated: false`. No stateful learning or baseline advantage is claimed. Large γ can produce a trivially frozen system, so the result is not by itself a useful-control or scientific-discovery result; any future capability claim must include nonzero signal and apparatus/control costs.

No blocking path-count, arbitrary-state norm, or certificate-verification defect was found in the current implementation and test receipt.
