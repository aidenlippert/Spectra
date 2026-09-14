# Retaining the coupling states: exact adaptive Schur bounds

The larger-perturbation precision problem is now resolved on both tested ten-mode fixtures. At strength 1/100, the connected-cycle interval has width **9.197032362665529e-8** and the mixed-interaction interval has width **8.325064374031218e-8**. These improve the preceding six-dimensional bounds, whose widths were approximately 9.61e-4 and 1.08e-3. The new construction retains the states reached by the perturbation and treats the actual Hamiltonian on that enlarged space exactly.

This result establishes finite-instance compression. It does not establish a compact boundary for the entire physical-marginal cone, quartic-relaxation optimality, or scalable reference discovery.

## Construction and proof

Use the same M=10, N=5 reference H0 and fixtures as [the reference-resolvent experiment](marginal_general_schur_results.md). Let Z span its six-dimensional symmetric space P. The verified reference decomposition supplies QH0Q >= cQ on Q=I-P. A rational Hermitian-pair estimate supplies eta >= ||H-H0||.

Start with the coupling W=Q(H-H0)Z. Exact elimination constructs its H0-Krylov closure R and checks closure by substituting H0R back into the resulting basis. Set U=[Z,R]. Additional rounds add the H0-closure of leakage from the current U. The basis need not be orthonormal. Define

\[
G=U^TU,\qquad A=U^THU,\qquad
L=HU-UG^{-1}A.
\]

The checker verifies full rank through exact positive metric pivots and checks U^TL=0. Since U contains P and is H0-invariant, its orthogonal complement S lies inside Q. For q=c-eta-b>0,

\[
S(H-bI)S\succeq qS.
\]

The Schur complement therefore gives the sufficient condition

\[
\boxed{A-bG-\frac{L^TL}{c-\eta-b}\succ0
\quad\Longrightarrow\quad H\succeq bI.}
\]

Every acceptance pivot is rational. A separately evaluated integer-amplitude Rayleigh witness supplies the upper bound. Optional reference-resolvent refinement replaces L^TL/q with the exact reference response on the leakage, using an automatically discovered and exactly verified annihilating polynomial. The strongest results below use enlarged spaces with the scalar complement bound.

## Accepted intervals

| Fixture / strength | Retained dimension | Certified width |
|---|---:|---:|
| Cycle, 1/1000, one round | 14 | 8.815006208697191e-10 |
| Mixed, 1/1000, one round | 36 | 1.4971624774245032e-9 |
| Cycle, 1/100, one round | 14 | 2.116834216764047e-5 |
| Mixed, 1/100, one round | 36 | 2.1972817593562723e-5 |
| Cycle, 1/100, one round plus resolvent | 14 | 4.2426812743194344e-6 |
| Mixed, 1/100, one round plus resolvent | 36 | 4.55615412567866e-6 |
| Cycle, 1/100, two full rounds | 32 | **9.197032362665529e-8** |
| Mixed, 1/100, targeted second round | 42 | **8.325064374031218e-8** |

The strongest intervals, shown as rounded decimal endpoints, are

\[
3.279699955076861\le E_{\rm cycle}\le3.279700047047185,
\]
\[
3.2805774206\le E_{\rm mixed}\le3.280577503850644.
\]

Exact endpoints are stored in the certificates. Energies use the fixtures' interaction units. Both Hamiltonians break all five individual pair charges and global left/right exchange.

Full enrichment gives dimensions 6→14→32 for the cycle and **6→36→130** for the mixed fixture. The latter is an exact dimension diagnostic only; no 130-dimensional positivity solve was attempted. Targeting one physical direction gives **6→36→42** for the mixed case.

## Selecting a direction without trusting a floating certificate

The mixed second round uses the leakage Lx of a proposed low-energy state, then closes this one vector under H0. Numerical diagonalization proposes the state; it does not establish the bound.

Direct rounding of coefficients in the elimination basis failed badly because its metric condition number was approximately 6.24e20. A retained Ritz value near 3.28058444 became an exact rounded-state Rayleigh value near 3.98334. This run still produced a valid lower certificate, but its interval width was a weak 0.0373239. It remains recorded as an accepted weak experiment, not a false certificate.

The successful selector instead rounds the state in physical Fock coordinates, projects that rational state exactly onto U using G^-1 U^T, and clears denominators to obtain an integer coefficient recipe. Its numerical Ritz value is 3.2805844370809556; its exact projected-state Rayleigh value is approximately 3.2805844370812784. The largest targeting integer uses 95 bits. The final bound is accepted by one exact positivity check after a numerical lower-bound proposal.

An alternative exact Gram-Schmidt construction improved rounding quality but generated very large rational numbers. It was stopped after an observed 300 seconds and produced no accepted certificate. This is a resource-limited experiment, not an infeasibility result. The default targeted path uses physical-coordinate rounding and exact projection.

The strongest cycle run took 8.73 seconds and the successful mixed run took 24.22 seconds in individual observed runs. These are not controlled scaling benchmarks. A subsequent removal of repeated matrix-action construction has not been separately timed.

## What remains

A subsequent exact first-round growth diagnostic extends the declared fixtures to 12, 14, and 16 modes. For m=M/2, it uses an (m+1)-dimensional signed symmetric embedding, centered matched perturbations s(i-(m-1)/2), cross edges weighted -s(1+i/m), and the same mixed local terms with right-hand mode indices shifted with m. The ten-mode construction is explicitly checked against the original fixture.

| Modes | Full fixed-particle sector | Cycle retained dimension / support | Mixed retained dimension / support |
|---|---:|---:|---:|
| 10 | 252 | 14 / 112 | 36 / 128 |
| 12 | 924 | 17 / 256 | 46 / 288 |
| 14 | 3,432 | 20 / 576 | 56 / 640 |
| 16 | 12,870 | 23 / 1,280 | 66 / 1,408 |

These are exact closure ranks, not energy certificates at the larger sizes. The observed dimensions follow 3m-1 and 10m-14 over these four sizes; no all-size rank formula is proved. Independent substitution checks H0-invariance and orthogonality to the reference space. The 16-mode mixed diagnostic took 8.01 seconds in one run. The script and all eight receipts are saved under `results/marginal_enlarged_schur/dimension_probe.py` and `dimension_m*_*.json`.

A subsequent [defect-Dicke implementation](marginal_defect_dicke_results.md) now supplies an exact implicit encoding for the first retained space and replays both first-round enlarged lower proofs. Its atoms fix local pair occupations and retain a symmetric spectator right count. The [complete implicit extension](marginal_implicit_certificate_results.md) also now ports targeted second enrichment and upper witnesses, reaching intervals below 1e-7 on both fixtures. Accurate, efficient higher-order compression at larger sizes remains a target. A small matrix dimension alone does not supply efficient construction.

The retained positivity matrices are small relative to the 252-dimensional particle sector, but their vectors are not uniformly sparse: the strongest cycle basis touches 232 configurations and the mixed basis 236. Their upper witnesses contain 252 amplitudes. The current numerical selector constructs the full particle-sector Hamiltonian. Thus these experiments compress the proof matrix while retaining substantial full-sector construction cost.

The remaining central tasks are to discover usable reference spaces and complement bounds beyond the specially structured H0; control growth of retained dimension, configuration support, and rational bit size; and preserve accuracy when the perturbation norm consumes the available complement gap. General chemistry also requires controlling finite-basis and physical-model error separately. None of these follows from the present finite certificates.

## Verification and artifacts

Implementation: `experiments/marginal_enlarged_schur.py`. All nine completed certificates under `results/marginal_enlarged_schur/` passed independent standard-library exact replay, including binding the stored Hamiltonian to its expected fixture. The ninth is the weak legacy targeting experiment described above. The partial orthogonal experiment and the mixed 130-dimensional diagnostic are not counted as certificates.

```sh
python3 -S -m experiments.marginal_enlarged_schur --verify results/marginal_enlarged_schur/cycle_1_100_rounds2/certificate.json
python3 -S -m experiments.marginal_enlarged_schur --verify results/marginal_enlarged_schur/mixed_1_100_rounds2_targeted_physical/certificate.json
```

The complete marginal test suite passed **120 tests in 30.510 seconds**. New checks cover exact closure and dimension limits, nonorthogonal metric projection, positivity gates, second-round targeting, fixture binding, and rejection of false lower bounds and corrupted upper witnesses. Exact residues and pivots are retained in the replay receipts.
