# Primary-source boundary map for the general matter mission

Scope: bounded audit of the claims in the mission attachment (read 2026-09-10). The links below are to the original paper or an author/institution repository for the paper. Dates are publication/version dates reported by those primary records.

## MatterGen: conditional generation and a synthesis proof of concept

**Primary source and date.** Claudio Zeni, Robert Pinsler, Daniel Zügner et al., “A generative model for inorganic materials design,” *Nature* 639, 624–632 (published 2025-01-16; version of record 2025-03-12), DOI [10.1038/s41586-025-08628-5](https://doi.org/10.1038/s41586-025-08628-5). The Nature article is the authoritative journal record; the author list begins Zeni, Pinsler, Zügner, Fowler, Horton, Fu, Wang, Shysheya, Crabbé, Ueda, Sordillo, Sun, Smith, Nguyen, Schulz, Lewis, Huang, Lu, Zhou, Yang, Hao, Li, Yang, Li, Tomioka and Xie.

**What it proves/demonstrates.** MatterGen is a diffusion model for periodic inorganic crystal structures (atom types, coordinates, lattice). A pretrained base model is fine-tuned with adapters/classifier-free guidance for conditions including chemical composition, space group, and scalar mechanical, electronic, and magnetic properties. The authors report improved rates of stable/unique/new (SUN) generated structures and DFT relaxation behavior against selected prior generative baselines. They experimentally synthesize **one sampled structure** and measure a property within 20% of its target. They also report rediscovering experimentally verified ICSD structures held out from training.

**Boundary.** These are candidate-generation, computational screening/relaxation, and one material-specific synthesis/measurement demonstrations. The paper does not prove a general inverse map from arbitrary desired behavior to composition, structure, preparation history, operating policy, and reliable physical performance. It does not provide a universal synthesis planner, guarantee that generated candidates are reachable by a specified process, or give a distribution-free bound on the probability that a generated design will meet its target after fabrication. “Synthesizable” in the paper is evidence from the one proof-of-concept synthesis plus rediscovery statistics, not a general synthesis theorem.

## Huang–Chen–Preskill: efficient prediction of local properties of unknown quantum processes

**Primary source and date.** Hsin-Yuan Huang, Sitan Chen, and John Preskill, “Learning to Predict Arbitrary Quantum Processes,” *PRX Quantum* 4, 040337 (published 2023-12-06), [CaltechAUTHORS record and paper](https://authors.library.caltech.edu/records/kd4c1-par39), DOI [10.1103/PRXQuantum.4.040337](https://doi.org/10.1103/PRXQuantum.4.040337).

**What it proves.** The paper gives an efficient classical ML algorithm that, for a wide range of input-state distributions \(D\), learns to predict any local property of the output of an unknown \(n\)-qubit process, with small **average error over inputs drawn from \(D\)**. The algorithm can remain computationally efficient even when the unknown circuit has exponentially many gates; the theorem concerns the specified learning/prediction task and its sample/complexity analysis.

**Boundary.** “Small average error over \(D\)” is not a uniform guarantee over every input, adaptive design intervention, rare state, or distribution selected by an optimizer. The result concerns local output properties of quantum processes, not general material properties, reaction outcomes, preparation reachability, or physical fabrication. It therefore supports the mission’s warning about average-case prediction, but cannot be promoted to a guarantee of reliable invention under adaptive search.

## DFT complexity: the universal-functional barrier

**Primary source and date.** Norbert Schuch and Frank Verstraete, “Computational Complexity of interacting electrons and fundamental limitations of Density Functional Theory,” *Nature Physics* 5, 732–735 (2009); author record [arXiv:0712.0483](https://arxiv.org/abs/0712.0483) and DOI [10.1038/nphys1370](https://doi.org/10.1038/nphys1370). The arXiv record identifies the work and its 2010 revision; the journal article is 2009.

**What it proves.** They show that an efficient description/evaluation of the exact universal DFT functional would allow efficient solution of problems in QMA (in particular, the Hubbard-model ground-state problem in an external magnetic field), yielding a complexity-theoretic limitation. The result is a conditional hardness barrier: it does not say every practical DFT calculation is QMA-hard, nor that no useful approximation can ever be efficient.

**Boundary.** The theorem concerns ground-state electronic-structure functional complexity in the stated model, not the full mission: time-dependent response, finite temperature, reactions, processing histories, synthesis routes, and operating systems are outside the theorem. “Learned” does not by itself evade the barrier; any claimed efficient universal predictor needs a precisely stated promise, approximation notion, distribution, or restricted physical class.

## Closest rigorous inverse-design/active-experiment guarantee: SafeOpt

**Primary source and date.** Yanan Sui, Alkis Gotovos, Joel Burdick, and Andreas Krause, “Safe Exploration for Optimization with Gaussian Processes,” ICML 2015, *Proceedings of Machine Learning Research* 37:997–1005 (conference 2015), [PMLR paper](https://proceedings.mlr.press/v37/sui15.html) and [PDF](https://proceedings.mlr.press/v37/sui15.pdf).

**What it proves.** In a finite decision set, with noisy evaluations, a known safety threshold, a nonempty safe seed set, Lipschitz regularity, and a bounded RKHS norm/GP confidence condition, SAFEOPT samples only decisions certified safe with high probability. It converges to an \(\epsilon\)-near-optimal point in the safely reachable set (and gives a sample-complexity bound). This is a rigorous active-experiment/inverse-optimization guarantee under explicit assumptions.

**Boundary.** The guarantee applies to a static bandit-style decision set and a scalar unknown objective/safety function. It does not model multistep chemical transformations, stateful preparation histories, irreversible failures, changing environments, multiple coupled physical constraints, or unknown model misspecification. Its benchmark is the safely reachable optimum from the supplied seed under the assumed regularity; it is not a guarantee of global optimum or general physical synthesis.

## Exact open gap exposed by the boundary

Taken together, these four sources do not establish the combined capability required by the mission: from an arbitrary behavior specification and operating conditions, construct a physically executable composition/structure/preparation/operation policy, adapt experiments safely, and provide a calibrated guarantee that the resulting physical system satisfies the specification. MatterGen covers conditional candidate generation and one synthesis demonstration; Huang–Chen–Preskill cover average-case local quantum-process prediction; Schuch–Verstraete identify a universal-functional complexity barrier; SafeOpt covers safely reachable near-optimal active optimization under strong assumptions. The missing theorem/engineering capability exposed by this comparison is a validated bridge across **model applicability and uncertainty, adaptive search, multistep physical reachability, fabrication, and final behavior**. In particular, average prediction error or DFT agreement on a screened distribution cannot be silently converted into uniform guarantees on the novel candidates selected by an optimizer. This map makes no claim that a new combined theorem already exists, or that these four papers exhaust all relevant literature.

## Source-date notes

- MatterGen: journal DOI record, 2025 (the DOI and Nature article are the primary source).
- Huang–Chen–Preskill: PRX Quantum publication date 2023-12-06 in the Caltech primary repository.
- Schuch–Verstraete: Nature Physics publication 2009; arXiv record has later revisions.
- SafeOpt: ICML 2015 proceedings, published in PMLR volume 37.
