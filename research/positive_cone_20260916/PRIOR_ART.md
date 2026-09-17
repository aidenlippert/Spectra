# Positive-cone formulation for Hubbard lower bounds: prior art and scope

Date: 2026-09-16

## Exact structural precedents

Lieb's original Hubbard theorems establish the relevant half-filled/bipartite spin structure and uniqueness statements: [Two theorems on the Hubbard model, Phys. Rev. Lett. 62, 1201 (1989)](https://doi.org/10.1103/PhysRevLett.62.1201). The proof uses spin-reflection positivity after a particle-hole transformation on one spin species. This is the key reason a repulsive Hubbard problem can sometimes be represented by an operator acting on one spin sector.

The proposed map

\[
\mathcal L(C)=KC+CK^T-U\sum_i n_i C n_i+UN_\uparrow C
\]

is therefore best understood as a concrete matrix realization of an already established reflection-positive Hubbard reduction, not as a new positivity principle. The exact signs, particle-hole convention, boundary conditions, and half-filling assumptions must be checked against the chosen Hubbard Hamiltonian. The repulsive half-filled reduction used here does not supply the corresponding certificate for doped or frustrated repulsive systems, arbitrary molecular Hamiltonians, or general long-range two-electron integrals. Lieb's separate attractive-model theorem has a different filling scope.

## Why the proposed certificate is mathematically valid under its hypotheses

Assume:

1. \(\mathcal L\) is self-adjoint for the Hilbert--Schmidt inner product.
2. The lowest eigenmatrix satisfies \(C_0\succeq0\) and \(\mathcal L(C_0)=E_0C_0\).
3. A trial \(C\succ0\) satisfies \(\mathcal L(C)-\ell C\succeq0\).

Then

\[
0\leq \operatorname{Tr}\!\left[C_0(\mathcal L(C)-\ell C)\right]
=(E_0-\ell)\operatorname{Tr}(C_0C),
\]

and \(\operatorname{Tr}(C_0C)>0\) whenever \(C_0\neq0\) and \(C\succ0\). Hence \(E_0\geq\ell\). This is a cone-duality/Collatz-type certificate. It does not require a positive semigroup for this one-line implication if the PSD ground eigenmatrix and self-adjointness are already established. A positive semigroup is useful for proving or computing the Perron ground object, but is an additional hypothesis rather than a substitute for verifying the displayed matrix inequality.

Positive-operator Perron--Frobenius theory is established for cone-preserving maps, including trace-ideal settings: [Perron--Frobenius theory for positive maps on trace ideals](https://arxiv.org/abs/math-ph/0007020). For a Hubbard-specific discussion of spin-reflection positivity and its consequences, see Boretsky, Cohn, and Freericks, [arXiv:1712.02694](https://arxiv.org/abs/1712.02694). These results support the conceptual interpretation, but do not make a finite trial matrix \(C\) cheap or guarantee strict positivity/irreducibility for a Hubbard sector.

## Complexity caveat: one-spin reduction is not yet compression

For (N=8) at balanced spin, the one-spin configuration space has dimension

\[
d=\binom{8}{4}=70,
\]

so a dense \(C\) has exactly \(d^2=4900\) amplitudes. For \(N=2m\), the same construction has

\[
d=\binom{2m}{m},\qquad d^2=\binom{2m}{m}^2\sim \frac{16^m}{\pi m},
\]

before any additional symmetry. Thus reflection reduction can avoid assembling the full spinful Hamiltonian while leaving an exponentially large operator matrix containing all \(d^2\) amplitudes. It is a real algebraic reduction, but not evidence of scalable compression until \(C\), the action of \(\mathcal L\), and the PSD residual are all represented and verified with subexponential cost. The current finite target is reported as a certified width of \(4.77409316436\times10^{-7}t\); existence of a compact tensor residual-PSD certificate remains open.

## What would be new

The defensible advance is not the particle-hole transform, spin-reflection positivity, or the Collatz inequality separately. It would be a directly constructible family of PSD trial operators \(C\) with a compact tensor/MPO representation for which

\[
\mathcal L(C)-\ell C\succeq0
\]

is certified without forming the full \(d\times d\) matrix, and whose representation/error bounds remain controlled under increasing lattice size, doping or molecular perturbations. A local purification can guarantee \(C\succeq0\), but it does not automatically certify the residual PSD inequality. A sampled quadratic-form check is insufficient unless accompanied by a certified operator-norm or covering argument.

The strongest immediate test is therefore a matched benchmark: construct the reflection-positive \(C\) for a half-filled bipartite Hubbard chain, report \(d\), tensor/MPO bond dimensions, residual-PSD verification cost, and the gap between \(\ell\) and the exact ground energy. Then perturb away from the theorem's domain (doping, frustration, or non-bipartite hopping) and require the method to return a mathematically explicit failure or weakened bound. Success only inside Lieb's exact regime would be a useful certified solver component, not a general many-body breakthrough.

## Source limits

I inspected the linked primary records and abstracts/pages. I did not reproduce Lieb's proof, verify the exact convention of the proposed \(\mathcal L\) against a supplied Hamiltonian, or perform an exhaustive search of reflection-positivity and positive-map literature. The novelty assessment is therefore conditional on the stated hypotheses.
