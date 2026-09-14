# Symbolic CAR certificate extractor: verification requirements

For an exact normal-ordered identity

\[
H-bI=\sum_\alpha B_\alpha^\dagger B_\alpha+(\hat N-N)X+R,
\]

the sign is correct: on the fixed-(N) sector, ((\hat N-N)X) vanishes, so (H\succeq(b-\eta)I) whenever (R\) has operator norm at most η. The argument works for arbitrary mode count (M) and particle number (N) without sector enumeration, provided the identity is checked in the CAR quotient itself.

## Required checks

1. **Hermiticity.** Require (H=H^\dagger), (X=X^\dagger), and (R=R^\dagger). If (X) is number-conserving, it commutes with \hat N, so ((\hat N-N)X) is Hermitian. Without number conservation, use the explicitly symmetrized ideal term
   (\tfrac12\{\hat N-N,X\}); otherwise it can be non-Hermitian even when (X) is Hermitian.

2. **Canonical reduction.** Expand all products with the CAR relations and reduce to one fixed ordered normal form, including a fixed ordering convention for creation and annihilation labels. Combine duplicate words exactly. Coefficient equality before reduction is meaningless because CAR produces signs and delta contractions.

3. **Residual norm.** If (R=\sum_w r_w W_w) in canonical ladder monomials and every ladder operator has norm at most one, then
   [
   \|R\|\le\sum_w|r_w|=:\eta.
   ]
   This is a valid, basis-dependent but dimension-independent bound by the triangle inequality. It can be extremely loose; report both η and the residual coefficient ℓ1 norm. Do not use a coefficient norm before CAR reduction, and do not silently discard repeated-index contractions.

4. **Gram positivity.** Store (Q\succeq0) as an exact rational (LDL^\mathsf T) or rational factorization (Q=L^\mathsf TL). Then (w^\dagger Qw=\sum B_\alpha^\dagger B_\alpha\) is a literal positive operator. A floating-point eigendecomposition rounded to rationals is only a proposal until the reconstructed (Q) is re-expanded and its exact PSD factorization is checked.

5. **Word normalization.** Fix whether (w) contains ordered words, antisymmetrized pair/triple words, or normalized exterior-algebra bases. The same convention must be used in (Q), in (Lw), and in coefficient matching; otherwise hidden factors of (2) or (k!) invalidate the certificate while preserving an apparently plausible matrix shape.

6. **Number ideal.** It is sufficient to exhibit one exact multiplier (X) with the stated identity; no completeness theorem for all polynomials vanishing on the (N)-sector is needed. To use the ideal term, the target sector must be explicitly declared as (\hat N=N). For a direct sum of sectors, replace it by the appropriate polynomial/projector constraints or verify each sector separately. The identity itself must hold in the CAR algebra (or after a stated quotient), not merely after testing sampled Fock states.

7. **Edge sectors.** Enforce (0\le N\le M). The proof remains valid for (N=0) and (N=M), but dictionaries containing annihilators or creators may become rank-deficient; this is harmless if the exact identity and PSD Gram check still pass.

## Minimal soundness theorem

If the extractor passes the seven checks and the canonical residual satisfies (\|R\|\le\eta) by the reduced coefficient ℓ1 bound, then for every normalized state ψ with (\hat Nψ=Nψ,

\[
\langle\psi,H\psi\rangle
 = b+\sum_\alpha\|B_\alpha\psi\|^2+\langle\psi,R\psi\rangle
 \ge b-\eta.
\]

This is a universal lower-bound theorem, independent of Hilbert-space dimension. It certifies correctness, not that the symbolic certificate is short, that η is tight, or that finding the identity scales efficiently.

