# Structural conditions for compact fermionic certificates

The proposed mechanisms are useful only when their hypotheses are made explicit. “Sparse,” “local,” and “gapped” do not automatically imply a small SOS certificate or a safe operator-pruning rule.

| Mechanism | What the primary theorem actually says | Certificate consequence |
|---|---|---|
| Chordal/low-treewidth structure | For a chordal sparsity graph, a partial Hermitian matrix has a PSD completion iff its fully specified clique principal submatrices are PSD (Grone–Johnson–Sá–Wolkowicz, 1984, [DOI](https://doi.org/10.1016/0024-3795(84)90207-6)). | This exactly decomposes a *PSD-completion problem with a known sparsity pattern* into clique constraints. It does not say a dense fermionic Gram matrix has a chordal pattern, nor that a sparse SOS truncation equals the original SOS cone. CAR normal ordering and ideal quotienting can create fill-in. A valid route is to prove a chordal pattern for the chosen operator dictionary and separately bound all omitted words. |
| Exponential clustering | Hastings–Koma prove decay from a nonzero spectral gap for short-range spin/fermion systems, with locality and observable-separation assumptions ([paper](https://arxiv.org/abs/math-ph/0507008), [DOI](https://doi.org/10.1007/s00220-006-0030-4)). | A Lieb–Robinson light cone alone is a finite-speed commutator bound; it is not exponential clustering. Pruning long-range dual operators needs a gap (and interaction assumptions), plus a quantitative conversion from correlation decay to energy/certificate error. Gapless systems, critical points, long-range Coulomb interactions, and degeneracies invalidate the simple implication. |
| SU(2) symmetry | If the Hamiltonian and constraints commute with the group action, representation theory block-diagonalizes operators into invariant sectors; this is an exact change of coordinates. | Exact symmetry can reduce PSD block sizes and variable counts. It does not reduce the number of irreducible sectors enough to guarantee polynomial/linear scaling, and broken symmetry, spin–orbit coupling, or approximate numerical symmetry removes exactness. Report the sector and multiplicity growth explicitly. |

## Additional mechanisms with concrete sufficient conditions

**Approximate Markov/recovery.** Fawzi and Renner show that small quantum conditional mutual information controls distance to a recovered state ([CMP DOI](https://doi.org/10.1007/s00220-015-2466-x), [arXiv](https://arxiv.org/abs/1410.0664)); Sutter, Fawzi, and Renner give a universal recovery map ([arXiv](https://arxiv.org/abs/1504.07251)). For a tripartition (A!:!B!:!C), a measured or certified bound (I(A:C|B)\le\delta) can justify replacing a long-range marginal by a reconstructed one, with a trace-distance error controlled by (delta). This is a state/marginal statement, not automatically a dual SOS truncation theorem; one still needs an observable-norm bound translating state error into energy error.

**Frustration-free plus a local gap.** The detectability lemma applies to local projectors whose common ground space is frustration-free and gives convergence of alternating local projectors as a function of a spectral gap; the local-gap-to-global-gap route is sharpened by Gosset and Mozgunov ([paper](https://arxiv.org/abs/1512.00088)). This can support compact local certificates when the Hamiltonian is genuinely frustration-free or perturbatively close to one. Generic ab-initio molecular Hamiltonians are not frustration-free, so this is a family hypothesis, not a universal chemistry principle.

**Tensor-rank/area-law structure.** In one dimension, a uniformly gapped local Hamiltonian has an area-law ground state and efficient matrix-product-state approximations under the assumptions of Hastings’ area-law theorem ([paper](https://arxiv.org/abs/0705.2024)). This supplies a constructive low-rank state representation, but a low-rank primal state does not automatically produce a small *dual* SOS certificate. A bridge would require a theorem bounding dual witness complexity from the MPS bond dimension and the Hamiltonian’s local terms.

## What the H4/H6/H8/H10 ladder can establish

The ladder is a good falsification program, not evidence of an asymptotic theorem by itself. For each size, freeze the structural rule before solving, measure dictionary size, clique/treewidth, PSD variables, discovery time, certified interval width, and omitted-operator remainder. Compare against the full rank-2 problem wherever feasible. A valid scaling claim needs a uniform rule and a proved or empirically stress-tested remainder bound across held-out geometries and interaction strengths.

The strongest near-term theorem target is therefore conditional:

> Given a specified operator graph with bounded treewidth, an exact symmetry action, and either a certified gap/correlation-decay bound or a certified recovery error, the restricted fermionic certificate has clique-local verification cost and an explicit bound on the contribution of omitted operators.

Each condition is meaningful separately. None of the cited results establishes this combined fermionic theorem automatically.
