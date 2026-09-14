# Audit design for coefficient-only SOS discovery

The coefficient-only experiment is meaningful if discovery and verification both operate on canonical CAR polynomials, with no import or call path through `marginal_hopping`, occupation-basis matrices, sector projectors, or exact diagonalization. A six-mode replay demonstrates the extractor; an (M=8,N=4) solve tests whether coefficient matching itself remains tractable.

## Required algebraic contract

Represent every polynomial as a sparse map from canonical normal-ordered ladder words to exact rationals. Construct

\[
H-bI=\sum_a B_a^\dagger B_a+(\hat N-N)X+R.
\]

The coefficient solver may use floating-point SDP or least squares to propose (b,Q,X), but export must reconstruct rational factors and re-expand every product in the CAR quotient. Verification must check (H,X,R) Hermitian, every square Gram factor exactly, (X) charge zero, and the residual coefficient map exactly. The safe lower bound is (b-\eta), with η equal to the reduced residual coefficient ℓ1 norm.

Do not treat a coefficient-array norm before CAR reduction as η: duplicate words and contraction terms can change it. Do not infer global PSD from sampled states. A rational (LDL^\mathsf T) or explicit rational factor list is the accepted PSD certificate.

The number ideal is valid because (X) is Hermitian and charge zero, hence commutes with \hat N. The implementation should either enforce those properties structurally or reject the certificate. A free (X) can absorb much of the matching problem, so report its term count, degree, coefficient bit length, and a control with (X=0).

## Controls that distinguish algebra from copied structure

Use two four-mode blocks, (M=8,N=4), with intra-block density repulsion and bounded inter-block hopping. Solve at (t=0), a small rational (t), and a generic perturbation that breaks block exchange symmetry. The Hamiltonian should be generated independently in polynomial form. Do not use a known ground energy as an optimization constraint; if an upper witness is supplied, label it separately from the coefficient-only lower result.

At minimum compare:

- quadratic words only;
- quadratic plus local cubic words;
- a bounded mixed cubic dictionary;
- the same dictionaries with (X=0);
- randomized dictionaries with matched word count and degree.

For (M=8), record canonical equation count, unknown Gram entries, (X) terms, solve time, rationalization time, residual ℓ1, and safe lower bound. A successful lower certificate is useful only when its size is reported and its residual is independently replayed.

## Scaling and completeness traps

Increasing (M) by disjointly copying a solved block is not a scaling experiment. The (M=8,N=4) instance must contain cross-block terms or asymmetric coefficients and must be solved from its own coefficient equations. Full cubic dictionaries can approach sector completeness even without explicit sector matrices; therefore cap dictionary support and report the fraction of all generated words included. A large (X) or dense Gram factor can hide the same exponential description in the multiplier.

The strongest next test is a sequence of genuinely coupled (2L)-mode systems at fixed local dictionary budget, with coefficient-only verification. Plot η and certificate description length against (L). An apparent constant gap with growing coefficient bit length or word count is not compression.

## Soundness receipt

For every accepted artifact, include a hash of the source Hamiltonian polynomial and the exported certificate, exact counts of canonical terms, maximum ladder degree, Gram factor rows/nonzeros, (X) support, and η. Include a deliberately corrupted coefficient and a non-Hermitian/non-charge-zero (X) negative test. This makes the result a reproducible algebraic certificate rather than an SDP transcript.

