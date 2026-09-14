# Certificate scaling frontier

## Scope and evidence

The relevant comparison is not “who first used SOS.” The test is whether a method gives a *small, discoverable, independently checkable* certificate with a controlled error as orbital/system size grows. Public evidence supports several bottlenecks, but it does not establish that Google’s methods or Rubin–Low–DePrince have reached a universal hard wall. Claims about competitors should therefore be phrased as measured limits of a stated implementation or benchmark family.

| Work | What is actually demonstrated | Scaling implication / open gap |
|---|---|---|
| Rubin, Low & DePrince, 2026, *Near-frustration-free electronic structure Hamiltonian representations and lower bound certificates* ([arXiv](https://arxiv.org/abs/2602.05069), [JCTC DOI](https://doi.org/10.1021/acs.jctc.6c00318)) | Connects weighted SOS to the dual v2RDM program; gives explicit Hubbard and electronic-structure constructions, including rank-2 algebras; reports molecular and Fe–S numerical benchmarks. The JCTC text says the SDP used libSDP and PySCF integrals and compares full versus subset algebras. | The paper demonstrates useful representations and lower-bound behavior, not a theorem that optimized certificates stay small with system size. The remaining questions are basis/operator selection, omitted-operator error, solver cost, and discovery without enumerating a large dictionary. |
| Veeraraghavan & Mazziotti, 2015, *Semidefinite programming formulation of linear-scaling electronic structure theories* ([PRA DOI](https://doi.org/10.1103/PhysRevA.92.022512)) | Exploits Hamiltonian sparsity for an effective one-electron SDP and demonstrates linear-scaling behavior on H chains up to H1500. | This is strong prior art for locality-aware SDP scaling, but it is a different effective one-electron problem; it does not certify a general interacting molecular SOS certificate. |
| Nakatani et al./Mazziotti large-scale v2RDM line ([DePrince lab overview](https://www.chem.fsu.edu/~deprince/research/rdm.html)) | v2RDM-CASSCF implementations report active spaces around (50e,50o) and many external orbitals. | “v2RDM cannot scale” is too broad. The unresolved comparison is rigorous certificate width and worst-case representability versus practical approximate energies and solver throughput. |
| Li et al., 2024, *Scalable semidefinite programming approach to variational embedding* ([JCP DOI](https://doi.org/10.1016/j.jcp.2024.113041)) | Uses localized clusters, parallel local updates, global consistency updates, and translation invariance to reduce PSD projection cost; reports iteration counts independent of cluster/system size for the studied setting. | Locality and embedding provide a credible competitor to global SOS. The gap is a chemistry-grade, non-periodic, ab-initio error theorem linking cluster size/operator truncation to a certified energy interval. |
| Chaykin et al., 2020, *Rigorous Lower Bounds for the Ground State Energy of Molecules* ([JCTC DOI](https://doi.org/10.1021/acs.jctc.0c00497)) | Demonstrates rigorous numerical lower bounds with interval arithmetic and rounding control. | Independent verification of numerical SDP outputs is established prior art; our possible distinction must be the certificate format, automation, or scaling/error guarantee, not “rigorous lower bounds exist.” |

## Where a real advantage could live

The most defensible target is an **adaptive local certificate with an a priori omitted-operator bound**. Start with a locality/symmetry-selected dictionary, solve the restricted dual, and compute a residual witness for every excluded operator block. A theorem would relate the residual norm, interaction tail, spectral gap or correlation decay assumptions, and the resulting energy error. The crucial deliverable is a computable stopping rule: enlarging the dictionary must either improve the lower bound or produce a certified bound on the remaining improvement.

That target is stronger than saying the optimum happens to have few nonzero coefficients. It addresses the failure mode where a compact decomposition exists but the optimizer cannot find it, or where coefficient sparsity is basis-dependent and gives no control over omitted directions.

## Falsifiable differentiators

1. **Discovery versus post hoc compression.** On a predefined family of Hubbard and ab-initio Hamiltonians, compare adaptive dictionary growth against the full rank-2 SDP. Record wall time, peak memory, number of PSD variables, final certified gap, and the number of generated operators. A successful claim requires subquadratic or near-linear growth in generated local blocks at fixed target error, plus independent replay of every certificate.

2. **Omitted-operator certificate.** Hide a fixed set of operator blocks during optimization. The method must output an upper bound on the possible lower-bound improvement from those blocks. Test the bound on small instances where the full SDP is available. It fails if the hidden-block improvement exceeds the predicted remainder, even when the restricted certificate looks tight.

3. **Transfer.** Train/select the structural rule on one geometry or interaction-strength range, then apply it unchanged to held-out geometries, basis sizes, and a small Fe–S cluster. Require monotone certified intervals and report failure cases. A rule that only works after per-instance hand selection is a useful heuristic, but not the claimed scalable structural condition.

## Bottom line for the race

The opportunity is real: combine locality, symmetry, low-rank structure, and exact verification into a method whose *search cost and truncation error are both measured*. We should say “beat the baseline on a specified benchmark and certificate budget,” rather than “Google has stalled,” until a public benchmark establishes that comparison. A successful result would be a theorem-plus-implementation showing compactness, efficient discovery, and an omitted-operator error bound together; any one of those alone is insufficient.
