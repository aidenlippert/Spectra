# Symbolic verifier review

The symbolic certificates replay successfully with `PYTHONPATH=. python3 -S` (the bare `python3 -S` invocation needs `PYTHONPATH=.` because `-S` suppresses the repository path setup):

- `full_mixed_matched.json`: lower bound (0.550936271779), residual coefficient ℓ1 norm (6.4370835\times10^{-5}).
- `matched_symmetry_compressed.json`: lower bound (0.5509129170803333), residual coefficient ℓ1 norm (8.77255337\times10^{-5}).

Both certificates pass without NumPy, SciPy, Fock-space enumeration, or eigenvalue computation. The implementation checks exact CAR re-expansion, Hermiticity, homogeneous charge, integer factor rows, and rational residual coefficients. Since every ladder monomial has operator norm at most one, the reduced coefficient ℓ1 norm is a valid operator-norm upper bound. The lower bound (b-\eta) is therefore sound for every (M,N) admitted by the verifier.

The number-ideal term is multiplied as ((\hat N-N)X). The verifier requires (X) to be Hermitian and charge zero, so it commutes with \hat N and the ideal term is Hermitian and vanishes on the declared sector. The orbit compressor also re-expands its recipe exactly before verification. No material soundness bug was found.

The observed residual degree six is expected: cubic squares generate degree-six symbols, while the target Hamiltonian has degree at most four. The coefficient ℓ1 bound is conservative and basis dependent, but it remains valid after canonical reduction. These receipts establish algebraic correctness of the finite certificates, not compactness or scalable certificate search.

