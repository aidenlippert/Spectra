# Hunting the physical 2-RDM boundary

## Thesis

For a fixed particle number `N` and one-particle basis of size `M`, let

\[
 K_N=\{\Gamma(\rho):\rho\succeq0,\;\operatorname{Tr}\rho=1,\;
 \rho\text{ on }\wedge^N\mathbb C^M\}.
\]

The useful object is not a global description of `K_N`. For a two-body Hamiltonian
`K`, compute the exposed face

\[
 F_K=\arg\min_{\Gamma\in K_N}\langle K,\Gamma\rangle
\]

and its normal cone

\[
 N_{K_N}(F_K)=\{A:\langle A,\Gamma-\Gamma_*\rangle\ge0
 \text{ for every }\Gamma\in K_N,\;\Gamma_*\in F_K\}.
\]

The compiler’s dual output should be a short, verifiable element of this normal
cone. This is the directional version of N-representability: discover the few
supporting inequalities that matter for one physical Hamiltonian.

## Established geometry and algebra

The energy of every number-conserving two-body Hamiltonian is linear in the 2-RDM
(with one-body terms folded into the contracted 2-RDM). Thus variational 2-RDM
optimization is an SDP over outer approximations to `K_N`; its dual searches for
positive operator combinations. Mazziotti’s dual-cone formulation explicitly
computes lower energy bounds and identifies the 2-RDM as a Lagrange multiplier in
the generalized Hellmann–Feynman formulation:

https://arxiv.org/abs/2103.17155

The physical cone is dual to the cone of two-body operators nonnegative on every
`N`-fermion state. A valid support certificate is therefore an operator identity

\[
 \widehat A-bI = \sum_\alpha B_\alpha^\dagger B_\alpha
                 +\sum_j\lambda_j C_j,
\]

where the `C_j` vanish on the fixed-`N` sector (number constraints, CAR identities,
and symmetry identities). Taking expectations gives
`<A,Gamma> >= b` for every physical marginal. Low-degree `B_alpha` produces a
constraint visible to a 2-RDM after normal ordering. The familiar `D,Q,G,T_1,T_2`
conditions are instances of this positive-polynomial construction, rather than
the final object itself.

The complete N-representability problem is computationally hard, so a claim that
the whole cone has a short description is the wrong conjecture. The directional
claim is narrower: the normal cone needed by chemically relevant `K` may contain
short certificates. The dual-cone variational paper above and Mazziotti’s complete
hierarchy construction provide the relevant starting points:

https://arxiv.org/abs/2304.08570

## The attack: primal state -> nullspace -> dual certificate

Given a high-quality physical primal state `|Psi>` and its 2-RDM `Gamma_Psi`:

1. Solve a low-order SDP to obtain `Gamma_*` and a lower bound. Record every nearly
   active positivity block (`D,Q,G,T_1,T_2`, symmetry blocks) and its small-eigenvalue
   eigenspaces.
2. Use `|Psi>` to identify near-annihilators. Search for low-degree fermionic words
   `B` with `||B|Psi>||` small, while requiring that `B^dagger B` normal-orders into
   the allowed two-body operator space plus fixed-sector identities.
3. Combine those words into `S=sum B^dagger B` and solve for coefficients so that
   `S` is aligned with the Hamiltonian direction modulo the current SDP constraints.
   This is a small least-squares/eigenproblem in operator space, followed by exact
   CAR expansion and interval/rational residual verification.
4. Add only certificates whose inequality is violated by `Gamma_*` (or whose dual
   objective raises the lower bound), then re-solve. Stop only when a separately
   evaluated physical state obeys
   `E_lower <= E_0 <= E_upper` with a certified gap.

The key signal is complementary slackness. At an exact optimum, a dual positive
operator `S` satisfies `<Psi|S|Psi>=0`, hence every positive summand obeys
`B_alpha|Psi>=0`. The primal tensor network is therefore not merely an upper-bound
generator: it supplies approximate kernels in which the next short certificate
should live. Conversely, a dual certificate identifies the face on which the
primal ansatz must concentrate.

### The gap identity and the one-sided-kernel trap

Suppose the verified identity on the `N`-particle sector is

`H - b I = sum_alpha B_alpha^dagger B_alpha + R`, with `||R|| <= eta`.

For a normalized variational state `|Psi>`,

`E_Psi - b + eta >= sum_alpha ||B_alpha Psi||^2`.

More sharply, `E_Psi-b = sum ||B_alpha Psi||^2 + <R>`, so the certified width
from this certificate is `E_Psi-b+eta`. This is the exact reason a primal state
can seed a near-annihilator dictionary: any small gap forces the weighted sum of
one-sided residual norms to be small.

There is a sharp trap. A construction based on an anticommutator, or on trying to
make a mode both empty and full, can demand that *both* `B Psi` and `B^dagger Psi`
be small. A bare fermion cannot satisfy this: `{a,a^dagger}=I` implies
`||a Psi||^2+||a^dagger Psi||^2=1` for every normalized state. The dictionary must
use homogeneous charge-changing dressed operators, such as `a + a^dagger a a`,
or paired excitation operators, and accept one-sided kernels selected by the
reference sector. Number conservation is a useful block symmetry, not a
requirement on every certificate word.

The practical replacement is to learn a positive Gram matrix over dressed words
`B_mu = sum_w c_{mu,w} w`, with the state-dependent objective

`min_C Tr(C G_Psi) + lambda * support(C)`, `C >= 0`,

subject to explicit Hamiltonian-matching and scale constraints; without them the
optimizer chooses `C=0`, making sparsity vacuous.

where `(G_Psi)_{uv}=<Psi|w_u^dagger w_v|Psi>`. Cubic-word Gram entries generally
contain 3-RDM information and are not available from a 2-RDM alone. Retain only
combinations whose CAR reduction is two-body (or whose degree-6 lift is admitted),
then verify the identity. Orbital rotations improve only this restricted sparse
search basis; the full fixed-degree SOS cone is basis invariant.

Hastings gives the decisive calibration: degree-4 SoS/2-RDM fails to reproduce
second-order perturbation theory for quartic fermion systems, while degree-6
SoS/3-RDM does, and a useful degree-6 fragment can be cheaper than the full
hierarchy (https://arxiv.org/abs/2205.12325). His later self-consistent perturbative
dual ansatz improves speed and accuracy, while explicitly warning that the tested
ansatz is unsuitable for quantum-chemistry Hamiltonians without modification
(https://arxiv.org/abs/2412.03564), and his June 2026 work extends the perturbative
SOS analysis (https://arxiv.org/html/2606.31765v1). The attack should therefore include a
*selective degree-6 escape hatch*: promote only the near-kernel words whose
two-body projection cannot certify the observed gap.

## What to learn

The learned object should be a sparse, symmetry-adapted map

\[
 (K,\Gamma_*)\mapsto \{B_\alpha,\lambda_\alpha\}
\]

with three hard filters: (i) CAR-normal-ordering lands in the two-body observable
space; (ii) positivity is exact or has a certified operator-norm residual; and
(iii) the resulting cut improves the directional dual bound. Train on violated
SDP points, not on energies alone. The loss should reward bound improvement per
certificate support and penalize operator degree, rank, and spatial/orbital support.

Orbital adaptation is essential. Localize or rotate orbitals using the primal
1-RDM/2-RDM, then search fragments, active spaces, and interfragment pair modes.
Certificates should be assembled from symmetry irreducible blocks (particle number,
spin, momentum, point group), because positivity and nullspaces decompose there.
The right representation is consequently a *directional normal-cone atlas*: each
Hamiltonian sector gets a small set of local algebraic charts, rather than one
universal list of inequalities.

## The theorem target

For a family of electronic Hamiltonians with target accuracy `epsilon`, prove or
measure that there exists a certificate set of total description length `L(K,epsilon)`
such that

\[
  \min_{\Gamma\in\mathcal R_L}\langle K,\Gamma\rangle
  \ge E_0(K)-\epsilon,
\]

where `R_L` is the SDP region cut by those certificates, and `L` grows gently with
system size for the chosen family. The strongest form asks for certificates whose
`B_alpha` supports remain bounded by correlation length or fragment complexity.
This is a conjecture, not an established theorem.

## Failure condition

The bet fails if chemically relevant exposed faces require certificates with degree,
rank, or support growing exponentially, or if near-annihilators of a good primal
state do not yield valid two-body projections. Degenerate ground spaces give a
larger exposed face but do not obstruct energy in the single Hamiltonian direction;
a normal-cone basis is needed for state selection or nearby observables. Finally,
small SDP duality gaps do not certify physicality unless the cuts are proven valid,
and a 2-RDM alone does not close real-time dynamics or finite-temperature evolution.

## Primary anchors

- Liu, Christandl, Verstraete, “N-representability is QMA-complete”: https://arxiv.org/abs/quant-ph/0609125
- Mazziotti, “Dual-Cone Variational Calculation of the 2-Electron Reduced Density Matrix”: https://arxiv.org/abs/2103.17155
- Mazziotti, “Quantum Many-body Theory from a Solution of the N-representability Problem”: https://arxiv.org/abs/2304.08570
- Delgado-Granados and Mazziotti, “Direct Variational Calculation ... via Semidefinite Machine Learning”: https://arxiv.org/abs/2603.05524

The 2026 semidefinite-machine-learning preprint is especially close in spirit: it
learns boundary information from molecular data while retaining an SDP structure.
It supports the direction, but does not establish scalable exact certificates.
