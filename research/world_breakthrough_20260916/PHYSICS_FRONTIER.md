# Physics frontier: a benchmark beyond certified small-molecule energies

## Finding

The clearest field-level benchmark is a **finite-temperature, thermodynamic-limit determination of the superconducting phase boundary of the doped repulsive square-lattice Hubbard model**, including the symmetry and competing-order context, with controlled uncertainty and independent validation. A convincing result would specify (U/t), density (including the near-1/8 regime), and (T/t), then establish by finite-size/bond-dimension/temperature extrapolation whether d-wave pair correlations have a genuine 2D finite-temperature transition (BKT behavior, or a rigorously supported absence of one), while simultaneously reporting spin/charge stripe correlations and thermodynamics.

This is a stronger benchmark than another accurate energy: it requires resolving an emergent phase in a sign-problematic, strongly correlated 2D many-body system, in the thermodynamic limit, with an error budget that survives method changes. The core requirement is one independently controlled physics result; agreement with a second method or microscope would strengthen it but is not logically required. The target observables are pair susceptibility/correlation-length scaling, helicity stiffness or an equivalent BKT diagnostic, and competing SDW/CDW structure factors.

## What current primary results establish—and where they stop

* Qin *et al.* combine constrained-path AFQMC and DMRG for the doped 2D Hubbard model and use cross-method checks, wide systems, and periodic geometries to reach the thermodynamic limit at **ground state**. In the pure nearest-neighbor model near optimal doping they report no superconductivity in the studied moderate-to-strong-coupling regime. This is an important benchmark for control, but it is not a finite-temperature phase-boundary result and depends on constrained-path AFQMC away from the sign-problem-free point. [Phys. Rev. X 10, 031016 (2020), arXiv:1910.08931](https://arxiv.org/abs/1910.08931)

* Xiao, He, Georges, and Zhang develop constrained-path AFQMC for **finite-temperature** doped 2D Hubbard systems. At dopings 1/5, 1/8, and 1/10 they resolve temperature-dependent spin and charge order, including stripe formation and a finite-temperature charge-order transition. The paper demonstrates strong control of competing order, but it does not establish a finite-temperature superconducting transition or a thermodynamic-limit d-wave condensate. [Phys. Rev. B 106, 075147 (2022), arXiv:2202.11741](https://arxiv.org/abs/2202.11741)

* Sinha *et al.* use finite-temperature iPEPS directly on the infinite square lattice, down to (T/t\approx0.17), with bond dimension up to 29, and obtain spin/charge correlators and specific heat in the slightly doped Hubbard model. This is a thermodynamic-limit tensor-network benchmark, but the reported regime is above the especially difficult low-temperature superconducting question and does not establish long-range d-wave order. [Phys. Rev. B 106, 195105 (2022), arXiv:2209.00985](https://arxiv.org/abs/2209.00985)

* Zhang *et al.* extend finite-temperature iPEPS to the doped 2D (t)-(J) model, reaching (T/t\sim0.1) and hole concentrations up to 1/4 in the thermodynamic limit. They find short-ranged d-wave pairing correlations over the studied window while resolving dopant-conditioned spin rearrangement. This is a recent and relevant limit-of-method result: thermodynamic-limit finite-(T) tensor networks can now expose local/intertwined correlations, but the published calculation does not demonstrate a superconducting phase transition. [arXiv:2510.04756 (submitted 6 Oct 2025)](https://arxiv.org/abs/2510.04756)

* Xu *et al.* report coexistence of superconductivity with partially filled stripes in a Hubbard model with next-nearest-neighbor hopping (`t-t'-U`), using tensor-network calculations. This is a major finite-temperature/ground-state superconductivity result in an extended model, but it does not settle the pure nearest-neighbor Hubbard benchmark proposed here. [Science, DOI:10.1126/science.adh7691](https://doi.org/10.1126/science.adh7691)

* Wang and Devereaux use numerically exact determinant QMC to identify finite-temperature d-wave signatures on the **electron-doped** side via imaginary-time-midpoint pairing correlations; they find no clear cooling signature for hole doping. The diagnostic is a finite-temperature advance, but the reported result is a signature rather than a thermodynamic-limit transition determination. [arXiv:2510.16616 (submitted 18 Oct 2025)](https://arxiv.org/abs/2510.16616)

* Zhang and von Delft introduce enhanced exponential TRG and reach `T ≈ 0.004t`, comparing finite-temperature data with zero-temperature iPEPS. This materially expands accessible temperatures, while the relevant finite-cylinder and thermodynamic-limit extrapolations remain distinct questions. [arXiv:2510.25022 (submitted 28 Oct 2025)](https://arxiv.org/abs/2510.25022)

As an experimental/theory control point, Feng *et al.* (2025) compare a doped, spin-imbalanced **attractive** Hubbard quantum gas with CP-AFQMC and find excellent agreement for short-range magnetic and charge correlations, then predict accessible FFLO precursors. This shows that microscope-to-AFQMC validation is becoming quantitative, but the attractive model and precursor regime do not answer the repulsive, d-wave thermodynamic-limit problem. [arXiv:2509.02688](https://arxiv.org/abs/2509.02688)

## Candidate discriminating test

Fix one public parameter line, for example the nearest-neighbor repulsive Hubbard model at (U/t=8), density near (n=0.875) (with density or chemical potential stated explicitly), and report a common grid of (T/t). The decisive submission would provide:

1. thermodynamic-limit extrapolations of the d-wave pair susceptibility, pair correlation length, and stiffness (or a finite-temperature scaling criterion that is equivalent to the BKT transition);
2. simultaneous spin and charge structure factors to distinguish superconductivity from stripe-driven finite-size enhancement;
3. a full systematic-error ledger: tensor-network bond dimension and environment error, AFQMC constraint/phase bias, Trotter or projection error, and finite-size extrapolation;
4. optionally, a blinded or preregistered comparison to site-resolved microscope correlators at the same dimensionless parameters, with no fitting of the answer to the experiment.

A field-level advance is achieved if one independently controlled calculation establishes (or rules out with a quantitative bound) the finite-temperature transition and its competing-order context in the thermodynamic limit. Agreement with another algorithm and experiment would make the claim substantially stronger. A result that only reports a variational energy, a short cylinder, a finite pairing amplitude without scaling, or a known formal identity would not meet this benchmark.

## Scope caveat

The benchmark is deliberately a **test**, not a prediction that superconductivity must occur. The scientifically valuable outcome can be a controlled upper bound or a controlled demonstration of no finite-(T) transition in the specified pure model. The novelty is the cross-validated thermodynamic-limit statement about an emergent phase and its competition with stripes.
