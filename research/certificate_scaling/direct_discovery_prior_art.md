# Direct-discovery prior art for sparse SOS certificates

This note classifies the closest established baselines for replacing a large PSD Gram cone with sparse, cheaply checkable pieces. It is not a novelty claim for the project; fermionic coefficient matching modulo CAR/number/spin ideals needs a dedicated prior-art search.

| Prior work | Established construction | Relevance to a fermionic certificate search |
|---|---|---|
| Ahmadi & Majumdar, 2019, DSOS/SDSOS ([SIAM DOI](https://doi.org/10.1137/18M118935X)) | DSOS restricts the Gram matrix to be diagonally dominant, giving an LP. SDSOS restricts it to be scaled diagonally dominant, giving an SOCP. These are inner approximations of SOS and trade solution quality for speed. | The “sparse ± two-word rank-one” cone is close in spirit to DSOS: a diagonally dominant PSD Gram matrix decomposes into nonnegative diagonal terms and (2\times2) pair blocks, each representable by weighted squares of sums/differences. Variable positive scaling corresponds to SDSOS. This is a known cone family, not by itself a new structural theorem. |
| Ahmadi & Hall, 2017, SOS basis pursuit ([AMS DOI](https://doi.org/10.1090/conm/685/13712), [arXiv](https://arxiv.org/abs/1510.01597)) | Iteratively solves LP/SOCP inner approximations, then pursues a better basis in which relevant SOS polynomials become DSOS/SDSOS. It applies to general SDP/SOS relaxations and reports discrete-optimization examples. | This is direct prior art for adaptive basis discovery. Our differentiator cannot simply be “iteratively add useful SOS atoms”; it must exploit fermionic algebra and prove a residual/omitted-operator bound or a family-specific complexity guarantee. |
| Waki, Kim, Kojima & Muramatsu, 2006, structured sparse SOS ([SIAM DOI](https://doi.org/10.1137/050623802)) | Uses correlative sparsity graphs to construct smaller SOS/SDP relaxations for sparse polynomial optimization. | Supports locality/block decomposition as a known route to scaling. A fermionic version would need to account for CAR normal ordering, ideal quotient collisions, particle/spin constraints, and the fact that Coulomb interactions are not strictly local in an orbital basis. |
| Wang, Li & Xia, 2018, SparseSOS ([preprint/PDF listing](https://www.researchgate.net/publication/327979820_Exploiting_Sparsity_in_SOS_Programming_and_Sparse_Polynomial_Optimization)) | Builds a cross-sparsity graph from polynomial support, chordally extends it, and solves block SOS problems over maximal cliques. | A useful baseline for support-driven block selection. It is polynomial/commutative SOS; it does not establish a CAR-ideal analogue or an energy error for omitted fermionic operator words. |
| Permenter & Parrilo, 2014, basis selection via facial reduction ([paper](https://www.mit.edu/~fperment/pdf/fr_sos.pdf)) | Facial-reduction and polyhedral procedures remove monomials that cannot occur in any feasible SOS representation and identify a smaller face of the cone. | This is prior art for provably eliminating basis elements before optimization. Any claim that ideal-aware basis pruning is unprecedented should be qualified: the fermionic quotient and certificate semantics may be new, but basis reduction itself is established. |

## Practical classification

There are three distinct ideas that should not be conflated:

1. **Fixed inner cone:** DSOS/SDSOS replaces PSD by LP/SOCP-checkable sufficient conditions. It is cheap, but can lose the optimum.
2. **Fixed structural sparsity:** chordal/correlative sparsity or facial reduction shrinks a known SOS basis. It can be exact for the selected relaxation, but does not automatically bound what was omitted.
3. **Adaptive discovery:** basis pursuit or column-generation-like loops search for a better basis. This can improve bounds, but the cited SOS basis-pursuit work does not give the fermionic, chemistry-specific omitted-operator certificate we need.

The project’s technically defensible opening is therefore narrower: specialize these mechanisms to a noncommutative fermionic coefficient map modulo CAR and physical ideals, then prove a **computable bound on the improvement available outside the selected dictionary**. The specialization may be novel; the general LP/SOCP cones, sparse support graphs, facial reduction, and adaptive basis pursuit are not.

## Minimal benchmark that would establish a real differentiator

For small fermionic Hamiltonians where the full rank-2 SDP is feasible, hide operator blocks, run the adaptive LP/SOCP or sparse-PSD search, and compare (i) certified lower-bound improvement, (ii) number of generated words/blocks, and (iii) the predicted omitted-block remainder against the exact full-SDP gap. A valid result must replay every certificate independently and fail loudly when the remainder bound is violated. Without that test, “direct discovery” is an implementation description rather than a scaling theorem.
