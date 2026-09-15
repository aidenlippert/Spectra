# Additional library candidates for Spectra

Checked against primary project documentation on 2026-09-15. These complement
the earlier Block2, PySCF, GPU4PySCF, OpenFermion, QCElemental/QCEngine and RDKit
list. A listed capability is not a measured Spectra improvement.

| Candidate | Concrete place in this engine | Evidence needed before adopting |
|---|---|---|
| [SuiteSparseQR / PySPQR](https://github.com/yig/PySPQR) | Replace the dense rank-revealing QR used to eliminate redundant ideal coordinates in proposal construction. | Same projected constraints and recovered operator coefficients; preserved residual gates; exact final certificate replay; measured setup and whole-solve time. |
| [CHOLMOD](https://github.com/DrTimothyAldenDavis/SuiteSparse) | Factor the shifted sparse normal system used as the iterative solver's preconditioner. | Accuracy and setup/repeated-solve cost on the actual H8/H10 normal matrices, then whole-solve comparison. |
| [NVIDIA nvmath-python / cuDSS](https://docs.nvidia.com/cuda/nvmath-python/latest/host-apis/sparse/index.html) | Access GPU sparse direct solvers through a supported Python interface. | Count planning, factorization, host/device transfers and repeated solves; compare with optimized CPU sparse methods. |
| [python-flint / FLINT / Arb](https://python-flint.readthedocs.io/en/latest/general.html) | Compiled exact rational arithmetic, exact matrix operations, or rigorously enclosing real arithmetic. | Equal exact endpoints for a rational substitution, or a proved enclosure throughout a ball-arithmetic implementation. The previous GMP probe does not establish a FLINT speedup. |
| [Forte](https://forte.readthedocs.io/) | Selected configuration interaction and multireference chemistry for candidate states and physical-model comparisons. | Spin/state identification, charged determinant counts, matched model definitions and numerical comparisons. Perturbative corrections are not automatically certified upper bounds. |
| [ipie](https://github.com/JoonhoLee-Group/ipie) | CPU/GPU auxiliary-field quantum Monte Carlo as an independent correlation method. | Statistical uncertainty, approximation bias and comparable Hamiltonians. A phaseless AFQMC estimate is not an exact two-sided Spectra certificate. |
| [geomeTRIC](https://geometric.readthedocs.io/en/latest/) | Molecular geometry optimization, constrained scans and transition-state workflows using chemistry-program energies and gradients. | Consistent energy/gradient conventions, converged geometry and a clearly defined physical observable. |
| [Basis Set Exchange](https://molssi-bse.github.io/basis_set_exchange/) | Standardized basis data, formats and version provenance. | Bind exact basis definitions and transformations to integral-generation receipts. This improves reproducibility, not the difficulty of solving a fixed Hamiltonian. |

The immediate performance candidates are sparse QR and the sparse normal solver.
The immediate verification candidate is compiled exact arithmetic. Forte and
ipie supply complementary scientific comparisons; geomeTRIC and Basis Set
Exchange improve the molecular workflow. Installing all of them is not itself an
engine improvement: each needs a named caller and a measured or validated benefit.

The accompanying NVIDIA experiment tests selected numerical candidates in an
isolated environment. The standard-library rational checker remains the
authoritative accepting path. See REPORT.md for actual outcomes and failures.
