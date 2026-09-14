# V6 primary-source receipt

## Strong baselines for forced coupled systems

**Van Overschee and De Moor, “N4SID: Subspace Algorithms for the Identification of Combined Deterministic-Stochastic Systems,” Automatica 30 (1994), 75–93.** [DOI/primary record](https://doi.org/10.1016/0005-1098(94)90230-5). N4SID identifies stochastic/deterministic state-space models from input/output data using projections, QR, and SVD; the paper states convergence and numerical stability for its setting. This is the appropriate strong linear baseline for forced coupled systems. It learns a latent state realization, but does not by itself provide physical variable identification or causal mechanism discovery.

**Brunton, Proctor, and Kutz, “Discovering Governing Equations from Data by Sparse Identification of Nonlinear Dynamical Systems,” PNAS 113 (2016), 3932–3937.** [DOI/primary article](https://doi.org/10.1073/pnas.1517384113). SINDy assumes the governing vector field is sparse in a supplied candidate-function library and demonstrates noisy, parameterized, time-dependent, and externally forced examples. This is a useful weakly nonlinear baseline when state measurements and a suitable library are available. Delay coordinates can reconstruct useful coordinates under Takens-type assumptions, but they do not guarantee that the reconstructed coordinates are physically interpretable or that an arbitrary partially observed forced system becomes sparse. A fair test should include N4SID/ARX, delay-coordinate linear models, and SINDy/SINDYc with identical data and held-out interventions.

## Future comparison systems

**Ellis et al., “DreamCoder: Bootstrapping Inductive Program Synthesis with Wake-Sleep Library Learning,” PLDI 2021.** [ACM DOI](https://doi.org/10.1145/3453483.3454080). DreamCoder learns a symbolic program library and neural search policy from a corpus of synthesis problems in a supplied DSL; its reported gains are on related benchmark domains. It is a future comparison for reusable symbolic abstractions, not a baseline for physical identification unless the physical DSL, noisy observations, and certificates are specified.

**Bowers et al., “Top-Down Synthesis for Library Learning,” POPL 2023, Article 41.** [DOI](https://doi.org/10.1145/3571234). Stitch learns abstractions from a corpus of programs in a DSL by corpus-guided top-down synthesis and reports large speed/memory gains over prior library learning on those benchmarks. It should be compared only as a symbolic-library learner, with search cost and DSL coverage made explicit.

**Udrescu and Tegmark, “AI Feynman: A Physics-Inspired Method for Symbolic Regression,” Science Advances 6 (2020), eaay2631.** [DOI/primary article](https://doi.org/10.1126/sciadv.aay2631). AI Feynman recursively combines neural fitting with dimensional analysis, symmetry, separability, and symbolic search; its headline evaluation uses the Feynman equation corpus and related benchmark equations. It is an offline equation-fitting baseline, not evidence of active experiment selection or physical-law discovery from interventions.

## Measured dataset check

The closest accessible measured candidate found is **Willke, Rahm, and Kabelac, “Experimental Investigation of Coupled Transport Mechanisms in a PEM Based Thermoelectric Energy Converter,” LUH dataset (2024)**, [dataset DOI](https://doi.org/10.25835/ok3269b9). The repository explicitly provides CC BY 4.0 licensing, raw open-circuit measurements, electrochemical impedance data, code, and a 37.3 MB download. It is independently measured coupled thermal/electrical data and is suitable for offline model comparison.

However, it does **not** clearly provide a controlled input/output time-series protocol with interventions chosen by the learner. Therefore it cannot support an active experiment-savings claim without an added experimental protocol. A larger alternative is the Dryad smart-building dataset ([DOI](https://doi.org/10.5061/dryad.73n5tb363)); it contains six years of measured electrical, heating, cooling, and weather series, but the reduced download is about 320 MB and is observational facility data, not a controlled active-experiment benchmark.

## Claim boundary

N4SID/ARX and SINDy test offline prediction and identification under declared excitation and model classes. DreamCoder, Stitch, and AI Feynman test reusable symbolic compression or equation fitting. None establishes that a learner discovers a new intervention, chooses a cheaper physical experiment, or transfers a mechanism outside its supplied vocabulary. Those stronger claims require held-out active interventions, matched conventional planners, and complete acquisition/search/verification cost accounting.

## Measured TCLab anchor (frozen receipt)

The repository’s TCLab notebook documents step and sine tests on the Temperature Control Lab: [model notebook](https://dowlinglab.github.io/pyomo-doe/notebooks/tclab-model/). The two raw CSVs were downloaded on 2026-09-10 from commit `d250c5b0625afb35007c075df1e0f7125e74017d`:

- [step test at frozen commit](https://raw.githubusercontent.com/dowlinglab/pyomo-doe/d250c5b0625afb35007c075df1e0f7125e74017d/data/tclab_step_test.csv): 901 rows, 26,515 bytes, SHA-256 `0c5fc5917a3693ca73a9240f0f40aea711fe1458853b88dd9c9db54fded776ec`.
- [sine test at frozen commit](https://raw.githubusercontent.com/dowlinglab/pyomo-doe/d250c5b0625afb35007c075df1e0f7125e74017d/data/tclab_sine_test_5min_period.csv): 901 rows, 20,272 bytes, SHA-256 `12abe51379ee107ee593986e1923aa849939aa2254412c551886708c26aa7f6d`.

Both have columns `Time,T1,T2,Q1,Q2`, no blank/non-numeric cells in the downloaded rows, and approximately one-second sampling: time spans are 0–900.02 s (step) and 0–900.01 s (sine), with observed increments 0.98–1.02 s. The step file ranges are `T1=22.84–62.16`, `T2=22.84–34.76`, `Q1=50`, `Q2=0`; the sine file ranges are `T1=22.20–70.54`, `T2=21.87–35.09`, `Q1=0–100`, `Q2=0`. These are measured input/output records suitable for an offline baseline; the constant Q2 columns and the exact excitation semantics must be verified from the notebook before calling them a rich coupled-identification benchmark.

The repository license downloaded at the same commit is BSD 3-Clause (SHA-256 `147dfb9bad5ff01ae6dad17cbe313ed47076dc3e64211785d19ee865afe0acb5`). This freezes provenance and permits reproducible reuse under its terms. The records do not, by themselves, establish an active experiment-selection protocol or experimental savings. Those claims require new learner-controlled interventions and matched cost accounting.

## AI Feynman 2.0 clarification

The user's named comparison is the later **Udrescu et al., AI Feynman 2.0: Pareto-optimal symbolic regression exploiting graph modularity (NeurIPS 2020)**, not only the original Science Advances paper. The [primary proceedings page](https://papers.nips.cc/paper/2020/hash/33a854e247155d590883b93bca53848a-Abstract.html) describes neural-gradient modularity tests and a tradeoff between symbolic accuracy and complexity. It was researched but not executed in V6. A known equivalent ARX baseline suffices to reject this particular privileged linear-history compilation; a future nonlinear symbolic-method claim would need stronger applicable comparisons.
