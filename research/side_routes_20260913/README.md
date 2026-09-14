# Complementary Spectra experiments — 13 September 2026 UTC

This side investigation completed three bounded experiments: an exact rational, fixed-particle-number chain certificate; a GPU benchmark on actual H10 Gram blocks; and a physical-model comparison against NIST. The chain route has a compact verifier but fails the desired accuracy scaling on the tested perturbed families. The A10 helps the larger tested blocks, with no demonstrated whole-solver speedup. The NIST comparison confirms that a highly converged finite-basis solver can still predict the wrong molecular geometry.

The main task's goal is to find an efficiently checkable condition under which accurate fermionic certificates can be **discovered and verified** without uncontrolled growth. The read-only [main-session snapshot](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/main_snapshot.json) shows that it is currently implementing SU(2) spin multiplets. It reports H10 Gram entries falling from 4.67 million to 1.22 million and the largest block from 725 to 500; those are dimension counts, not yet a measured runtime improvement. This side work occupies its own research/results directories and a separate directory on existing Lambda B.

**1. Exact positive-amplitude certificates with a small dynamic program.**

The implemented model is an open, nearest-neighbor, spinless fermion chain, restricted to exactly N particles:

\[
H=-\sum_{i=0}^{M-2}t_i(a_i^\dagger a_{i+1}+a_{i+1}^\dagger a_i)
  +\sum_{i=0}^{M-2}V_i n_i n_{i+1}+\sum_{i=0}^{M-1}h_i n_i+c,
\quad t_i\geq0.
\]

All inputs are rational. In the ordered occupation basis these nearest-neighbor hops have nonpositive matrix entries. Supply a strictly positive amplitude

\[
\psi(s)=\prod_i y_i^{n_i}\prod_i x_i^{n_i n_{i+1}},\qquad x_i,y_i>0.
\]

Define the local energy \(e(s)=(H\psi)(s)/\psi(s)\). The exact interval is

\[
L=\min_{|s|=N}e(s)\leq E_0\leq
U=\frac{\sum_{|s|=N}\psi(s)^2e(s)}{\sum_{|s|=N}\psi(s)^2}.
\]

The lower bound follows without a numerical eigenvalue computation. For each unordered pair of configurations with \(H_{st}<0\), form the matrix supported on s,t:

\[
K_{st}=(-H_{st})
\begin{pmatrix}\psi(t)/\psi(s)&-1\\-1&\psi(s)/\psi(t)\end{pmatrix}\succeq0.
\]

It has positive diagonal and zero determinant. Its off-diagonal entry is \(H_{st}\), so

\[
H-LI=\sum_{s<t:H_{st}<0}K_{st}+\operatorname{diag}(e(s)-L)\succeq0.
\]

The upper bound is the trial-state Rayleigh quotient. This argument is a standard positive-amplitude/stoquastic construction, not a claimed new general quantum theorem. The repo already proposed a related positive-amplitude route in `research/certificate_scaling/positive_amplitude_structural_route.md`; this implementation adds exact fixed-N evaluation for the declared fermionic chain. The broader stoquastic setting is discussed in [Bravyi et al., The Complexity of Stoquastic Local Hamiltonian Problems](https://arxiv.org/abs/quant-ph/0606140). Nonpositive off-diagonal entries alone do not make all ground-state problems easy.

The important computational step is that swapping neighboring occupations changes this amplitude by a ratio depending on at most four consecutive bits. Thus e(s) is a sum of local factors. The verifier computes its minimum with a dynamic program whose state stores the particle count and the last three bits. A parallel transfer recurrence accumulates the normalization and Rayleigh numerator, using exact fractions throughout.

This takes O(MN) rational arithmetic operations for fixed factor range, at most 8(N+1) states per layer, and O(M) local factor tables. It does not construct the exponentially large configuration matrix in the displayed proof. Rational bit lengths also matter: finite-range products and sums have polynomial bit growth in the input size. The receipt's `peak_rational_bits` records the observed minimum/normalization/numerator DP quantities, not every temporary allocation. The size guard limits the declared DP state budget; it is not a comprehensive memory limit for arbitrarily huge input integers.

The compact amplitude has O(M) rational parameters. This establishes cheap verification for a supplied amplitude. It does **not** establish an algorithm that discovers an accurate amplitude for a generic Hamiltonian, nor does it translate this certificate into the main cubic CAR/SOS file format.

**The experiment and its negative result.**

The [campaign](/Users/aidenlippert/Documents/Spectra/research/side_routes_20260913/chain_campaign.py) screened 40 rational parameter choices in each of five scenarios at M=12: 200 screens, then 25 evaluations at other sizes. Each scenario selected its smallest exact interval at M=12 and retained those same amplitude parameters at M=4,8,16,32,64. All 225 completed in 1.54 seconds of campaign wall time with four local CPU worker processes. These are parameter experiments across five scenarios, not 225 separate mathematical hypotheses.

| Scenario | Width, M=4 | M=8 | M=16 | M=32 | M=64 |
|---|---:|---:|---:|---:|---:|
| Interacting frustration-free control | 0 | 0 | 0 | 0 | 0 |
| Control + 0.0001 nearest-neighbor density interaction | 0.00005 | 0.00015 | 0.00035 | 0.00075 | 0.00155 |
| Repulsive interaction | 0.4412 | 1.1223 | 2.6419 | 5.8550 | 12.2993 |
| Varying bonds + interaction perturbation | 0.1192 | 0.3688 | 0.7612 | 1.6570 | 3.3394 |
| Varying site fields | 0.3385 | 1.3144 | 2.7812 | 6.0934 | 12.3274 |

Every width here is in the declared **abstract hopping energy unit**, not hartree. None is a claim of molecular chemical accuracy. See the exact fractions and work counts in the [campaign summary](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/positive_chain/summary.json).

At M=64 and N=32 the control has \(\binom{64}{32}=1,832,624,140,942,590,534\) occupation states. The dynamic program uses 246 states at its widest layer, 15,374 transitions, and 1,120 factor-table entries; its exact interval is [0,0]. This is a known, engineered frustration-free control, not a hard molecular example.

There is also an analytic explanation for the near-control failure. Take

\[
H_0=\sum_i[-(a_i^\dagger a_{i+1}+\text{h.c.})+n_i+n_{i+1}-2n_i n_{i+1}].
\]

Each bond is positive semidefinite and annihilates the uniform fixed-N amplitude. For \(H_\delta=H_0+\delta\sum_i n_i n_{i+1}\), \(\delta>0\), that same amplitude has local energy equal to delta times the number of occupied neighboring pairs. At half filling, its minimum is zero, while the uniform expectation of that number is N(N-1)/M. Therefore this witness's width is exactly

\[
U-L=\delta\frac{N(N-1)}{M}=\delta\left(\frac{M}{4}-\frac12\right).
\]

That width grows linearly despite a path interaction graph and cheap verification. This is a rigorous limitation of the specified witness, not a lower bound on the best possible certificate or a proof that other amplitude families fail. The 40-point screen likewise does not globally optimize all positive amplitudes.

The next useful mathematical target for this route is **control of the total local-energy variation as factor range grows**. For example, a proved error bound of the form M exp(-r/xi), combined with a verifier costing M N exp(O(r)), could justify choosing r proportional to log(M/epsilon) when xi is bounded. That estimate is a hypothesis, not something implied by correlation decay or proved by this experiment. Discovery of the factors would still need its own cost bound. Optimizing the certified width directly gives this route a falsifiable objective.

Verification: five focused tests passed, including exact comparison against an independently assembled CAR Hamiltonian in every particle sector of a nonuniform six-mode problem, numerical spectral containment, tamper/refusal tests, and the analytic extensive-width formula. A fresh `python -S` process replayed all 225 stored bounds in 1.99 seconds. This replay uses the same exact algorithm; the independent implementation is the small CAR test. See [tests](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/chain_tests.log) and [replay receipt](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/chain_replay.json).

**2. Actual H10 blocks on the existing Lambda A10.**

The [benchmark](/Users/aidenlippert/Documents/Spectra/research/side_routes_20260913/gpu_projection_gate.py) compared float64 positive-semidefinite projection on the eight actual H10 blocks of each tested size. The CPU used SciPy divide-and-conquer `eigh` with one BLAS thread. GPU methods used CuPy 13.6.0, either one matrix at a time or its batched interface. Timings include reconstruction and GPU input/output transfer, with two warmups and three measured trials. File loading and package installation are excluded.

| Eight blocks per call | CPU median | A10, one matrix at a time | A10, batched |
|---|---:|---:|---:|
| 225 by 225 | 0.03250 s | 0.06284 s | 0.11833 s |
| 725 by 725 | 0.58095 s | 0.25556 s | 2.13279 s |

The larger sequential GPU case is 2.27 times faster than this CPU baseline. The small sequential GPU case is 1.93 times slower, and the large batched GPU case is 3.67 times slower than CPU. Maximum absolute differences from the CPU projections were at most 3.25e-14. Input matrix entries reached approximately 0.00187 and 0.01001 respectively. These are numerical comparisons, not exact PSD certification.

This supports testing a hybrid dispatch for small versus large blocks. It does not yet support a whole-solver speedup claim or an assumed benefit from a stronger GPU. The main task's new spin decomposition changes the block sizes again, so benchmark its actual blocks before choosing a threshold. These kernel tests used roughly 293 MB of CuPy's memory pool; memory capacity was not the limiting issue here. The [GPU receipt](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/gpu/gpu_receipt.json) and [source metadata](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/gpu/source_metadata.json) contain source hashes and timings.

CuPy was installed only into `/home/ubuntu/spectra-side-20260913/vendor` on existing Lambda B, and the benchmark read the main run's raw matrices. The main solver and Python environment were not edited. No new instance was launched or existing instance terminated. The user has authorized considering stronger GPUs; a larger hardware run remains a next experiment, not an outcome of this one.

**3. A real NIST physical-model check.**

The [NIST Chemistry WebBook ground-state H2 row](https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1000) lists an equilibrium internuclear distance of 0.74144 angstrom. The frozen [reference record](/Users/aidenlippert/Documents/Spectra/research/side_routes_20260913/nist_h2_reference.json) preserves species, charge, state, quantity, units, source and compilation date. The displayed row supplies no uncertainty for that value, so the record leaves uncertainty unspecified.

The experiment optimized the neutral singlet H2 geometry with nonrelativistic Born-Oppenheimer FCI in three finite orbital bases, using PySCF 2.14.0 and SciPy 1.18.1. It includes nuclear repulsion in each geometry objective. The optimizer searched 0.60–0.95 angstrom; the final point was lower in energy than both tested neighbors at plus/minus 0.001 angstrom. Numerical state normalization, singlet spin and Hamiltonian residuals were checked. These checks do not constitute a formal global geometry certificate. [PySCF's FCI documentation](https://pyscf.org/user/ci.html) describes the numerical method.

| Basis | Calculated r_e, angstrom | Difference from NIST, angstrom |
|---|---:|---:|
| STO-3G | 0.73486523 | -0.00657477 |
| cc-pVDZ | 0.76089336 | +0.01945336 |
| cc-pVTZ | 0.74262280 | +0.00118280 |

The accepted run finished in 4.39 seconds. All evaluated FCI residual norms were below 9.01e-13 hartree. The cc-pVDZ geometry is farther from the reference than STO-3G; this small ladder must not be presented as monotonic convergence of every observable. Finite-basis error and omitted physical effects remain even after the finite electronic problem is solved very accurately. This measures one physical property of one species, not a general NIST validation suite or a force certificate.

An initial cc-pVTZ iterative solve failed its convergence check and was rejected. The accepted run uses full diagonalization of at most 784 determinants and retains the residual gate. The [initial log](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/nist_h2.log) is preserved; the [accepted receipt](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/nist_h2_dense/summary.json) contains every evaluated geometry. No pass threshold for agreement with experiment was invented.

The practical validation structure is now concrete: compare certificates against an independent solver **for the identical Hamiltonian**, and compare predicted observables against experiment **with their physical definitions and model errors recorded**. A small certified electronic-energy interval cannot by itself settle basis error, reaction kinetics, finite-temperature behavior, FeMoco, or superconductivity.

Reproduce the exact chain checks from `/Users/aidenlippert/Documents/Spectra`:

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -m unittest research.side_routes_20260913.test_positive_chain -v
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.side_routes_20260913.replay_campaign --directory results/side_routes_20260913/positive_chain --out /tmp/spectra-side-chain-replay.json
```

Run the NIST experiment with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv-molecule/bin/python -m research.side_routes_20260913.nist_h2_geometry --out <new-output-directory>`. The experiment refuses to overwrite an existing output directory. [Artifact hashes](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/manifest.json) identify the files delivered with this investigation.
