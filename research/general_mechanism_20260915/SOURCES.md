# What “solve the general quantum many-body problem” can mean

This note keeps the broad mission intact while separating three targets: (i) ground energy, (ii) a compact ground-state or reduced-state representation, and (iii) dynamics. A universally correct method may have exponential cost. The stronger claim of universally efficient classical computation for arbitrary Hamiltonians faces conditional complexity barriers. These sources do not prove that every physically relevant family is hard.

## 1. Arbitrary local Hamiltonians

**Kempe, Kitaev, and Regev, “The Complexity of the Local Hamiltonian Problem,”** [arXiv:quant-ph/0406180](https://arxiv.org/abs/quant-ph/0406180).

The promise problem asks whether the smallest eigenvalue of a polynomial-size sum of constant-locality terms is below `a` or above `b`, with an inverse-polynomial promise gap. It is QMA-complete already for 2-local Hamiltonians. “Local” here is algebraic, not necessarily geometric or chemical.

This does not rule out efficient algorithms for restricted physical families, special geometries, average cases, quantum algorithms, or certificates exponential only in a small correlation parameter. It also says nothing directly about dynamics. A general Spectra result must declare a Hamiltonian class and a parameter `r` for which construction, verification, and energy error are controlled, then show `r` stays manageable.

## 2. Fermionic consistency

**Liu, Christandl, and Verstraete, “N-representability is QMA-complete,”** [arXiv:quant-ph/0609125](https://arxiv.org/abs/quant-ph/0609125).

The paper proves QMA-completeness of deciding whether a proposed reduced density matrix is compatible with an `N`-fermion state. Thus local moments can satisfy CAR and positivity constraints while still failing global fermionic representability.

This does not show every molecule is hard or that fixed-order RDM methods cannot work on structured classes. It cautions against treating arbitrary locally consistent moments as a physical state. A sound SOS lower certificate needs no physical representation of relaxed moments; a promised tight, complete, or state-reconstructing scheme needs additional representability or controlled-error guarantees covering the omitted global consistency constraints.

## 3. Complete noncommutative positivity hierarchies

**Klep, Magron, Massé, and Volčič, “Upper bound hierarchies for noncommutative polynomial optimization,”** [arXiv:2402.02126](https://arxiv.org/abs/2402.02126).

This work studies minimizing the eigenvalue of a noncommutative polynomial subject to noncommutative polynomial inequalities. It gives SDP moment/SOS lower bounds converging under stated assumptions and complementary upper bounds from generalized eigenvalue problems. It is close in formal language to Hamiltonian SOS certificates.

The lower SDP is an outer relaxation: its feasible pseudo-moments need not be globally representable by a physical state. Soundness of a lower energy bound comes from relaxing the minimization domain in the correct direction. Global representability matters for tightness, convergence, and interpreting a pseudo-moment as a state; it is not required for every sound lower certificate. Completeness still means convergence as order increases, not a system-size-independent order or cheap level. A useful Spectra hierarchy must be generated and checked without materializing the full global moment/Gram object, while proving that its order, block size, or boundary rank suffices for a declared error class.

## 4. Noncommutative SOS and collective constraint modules

**Helton, Klep, and McCullough, “The convex Positivstellensatz in a free algebra,”** [arXiv:1102.4859](https://arxiv.org/abs/1102.4859).

For positivity on a matrix-convex set described by a monic linear pencil, the paper proves a weighted noncommutative sum-of-squares representation with degree bounds under its hypotheses. This is a rigorous precedent for SOS terms plus collective constraint multipliers.

The theorem is not automatically a fixed-number fermionic result: CAR, spin ideals, sector restrictions, and molecular residuals must be embedded and checked. Degree control is not coefficient-count or optimizer-cost control. The missing Spectra construction is a sound quotient/SOS translation with bounded generator support, residual, and discovery cost.

## 5. Renormalized collective lower bounds

**Kull, Schuch, Dive, and Navascués, “Lower Bounding Ground-State Energies of Local Hamiltonians Through the Renormalization Group,”** [arXiv:2212.03014](https://arxiv.org/abs/2212.03014).

Given a renormalization scheme, this paper constructs a tractable convex relaxation of feasible local density matrices. Coarse-graining maps eliminate many consistency constraints, while the remaining constraints yield rigorous lower bounds by linear optimization. The paper demonstrates bounds for 1D translation-invariant spin models and states extensions to other many-body settings as a direction.

The quality depends crucially on the chosen, target-tailored renormalization scheme. This is directly relevant to Spectra's local blocks, response channels, and collective boundary terms, but it is not a universal fermionic construction. The missing step is a fermionic/CAR-compatible coarse-graining map whose outer relaxation remains sound, whose omitted constraints have a certified effect on energy, and whose boundary complexity is provably smaller than the global operator-pair construction.

## What is already known, and what is missing

There are already complete, in-principle routes to universal correctness: full diagonalization, complete moment/SOS hierarchies, and sufficiently high-order exact certificates. Their cost can be exponential. There are also sound outer relaxations and terminal dual witnesses: they can certify a lower energy without representing a physical global state. What is not yet supplied is a general, efficiently constructible terminal certificate whose size tracks only the important correlations for arbitrary interacting fermionic Hamiltonians.

“Solve every finite many-body Hamiltonian with a polynomial classical cost” is not a defensible claim absent a major complexity-theoretic breakthrough in light of QMA-hardness and fermionic N-representability. That conditional barrier does not narrow the mission to molecules or a physical subclass; it identifies the distinction between universal correctness and universal efficiency. The research target remains to find whether Spectra can establish a broader constructive regime, and to state precisely where it cannot.

The required theorem and evidence are:

1. a sound quotient/SOS or response construction;
2. completeness or a controlled approximation theorem for the chosen regime;
3. a direct constructor whose work avoids the full global Gram matrix;
4. an explicit terminal condition or obstruction when compression stops helping; and
5. independent replay on coupled, changed, and held-out systems with discovery and verification costs charged.

Energy certification, state representation, and dynamics remain distinct. An energy interval does not provide a wavefunction; a compact state does not automatically certify energy; and neither supplies real-time dynamics without separate control of propagation and accumulated error.
