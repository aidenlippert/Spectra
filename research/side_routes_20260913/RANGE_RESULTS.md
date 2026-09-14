# Longer-range search and an algebraic alternative

This continuation produced an exact certificate for a known interacting chain through 10,000 sites, after a broader numerical search exposed the limitations of fitting local amplitudes. The useful new research direction is recognition of compact algebraic identities. The solved family is special; scalable discovery for unfamiliar molecular Hamiltonians remains open.

**The strongest result: a recognized algebraic family.**

Independent reference calculations returned ground energies -M/2 for the repulsive t=V=1 benchmark at M=8,12,16. A second, dense CAR implementation confirmed the small cases. Literature checking identified this as a known supersymmetric XXZ point, associated with [Yang–Fendley](https://arxiv.org/abs/cond-mat/0404682) and developed through local supercharges by [Hagendorf–Liénardy](https://arxiv.org/html/1612.02951).

The [implemented recognizer and verifier](/Users/aidenlippert/Documents/Spectra/research/side_routes_20260913/supersymmetric_chain.py) accept uniform positive hopping t=a, interaction V=a, a uniform field mu, and an additive constant c. For even M and N=M/2 they certify the exact energy c+mu N-a N. They check a local positive-factor identity and a compact witness that the lower bound is attained. The [short proof](/Users/aidenlippert/Documents/Spectra/research/side_routes_20260913/SUPERSYMMETRY_PROOF.md) specifies the algebra.

For M=10,000, N=5,000, a=1, mu=c=0, the accepted interval is **[-5000,-5000]**, in the abstract model's energy units. The witness JSON is 251 bytes excluding the separately supplied Hamiltonian. Recognition plus verification took 0.260 seconds in the recorded single local run. It checks 29,998 Hamiltonian coefficients and constant-size tensors, without enumerating occupation states. This is an input-length complexity argument backed by timings, not a broad hardware benchmark.

The [nine-size ladder](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/supersymmetry_ladder.json) covers M=4,8,12,16,32,64,256,1024,10,000. All nine certificates passed a fresh `python -S` [replay](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/supersymmetry_replay.json). Tests independently construct the full length-changing matrices for M=2 through 6 and compare their factorization exactly with CAR Hamiltonians in every particle sector. They also test scaled/shifted models and reject changed interactions, boundary fields, false saturation, and tensor tampering. This is a new implementation of known mathematics. It is not yet part of the main fixed-M CAR/SOS certificate format.

**What the larger search established.**

The [finite-range verifier](/Users/aidenlippert/Documents/Spectra/research/side_routes_20260913/finite_range.py) extends the previous pair-amplitude family to

\[
\psi(s)=\prod_i y_i^{n_i}\prod_{d=1}^{r}\prod_i x_{d,i}^{n_i n_{i+d}},\qquad r=0,1,2,3.
\]

A nearest-neighbor hop changes amplitudes within at most 2r+2 consecutive sites. The exact dynamic program retains the particle count and at most 2r+1 trailing bits. For fixed range its arithmetic-operation count is O(M N 2^(2r+1)); exact rational bit growth is also charged in the receipts. Increasing the range raises verification cost substantially.

The numerical discovery family uses one pair weight per distance and three fixed site features: alternating sites, endpoints, and integer-scaled local fields. These are a limited parameterization of the general verifier input. Seventy-two budgeted searches cover three scenarios, four ranges, three objectives, and two seeds. They evaluate **77,540 numerical proposals**. Each optimized proposal is rounded to positive rational weights and checked exactly, along with its starting point. The three objectives are minimum interval width, minimum variational upper energy, and maximum lower bound.

All parameters are fitted at M=12. Selection at each range also considers smaller-range solutions embedded by adding unit pair weights. The selected parameters are then fixed for evaluations at M=8,12,16,32,64: **180 exact transfer certificates**, including the training size. These are parameter studies in three abstract chains, not 72 separate physical hypotheses.

Both numerical campaigns ran on four CPU workers on existing Lambda B, at reduced scheduling priority. Their measured search phases took 4.34 and 2.56 seconds; their full runs, including exact transfer calculations, took 42.81 and 27.21 seconds. These times exclude upload/setup and subsequent replays/reference work. Numerical screening was cheap enough that a larger GPU was unnecessary for this experiment. Only 2 of 72 differential-evolution runs declared convergence before their budget; 52 of the subsequent Powell runs did. No global optimality conclusion is drawn from these searches.

The [first campaign](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/finite_range_campaign/summary.json) and [lower-bound campaign](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/finite_range_lower/summary.json) preserve costs, seed outcomes, parameters and exact fractions. Source snapshots preserve both runner versions.

**Separately optimizing the two bounds helps.**

The best lower and upper bounds need not use the same amplitude. The [combiner](/Users/aidenlippert/Documents/Spectra/research/side_routes_20260913/combine_bounds.py) selects the maximum lower and minimum upper claims, checks identical declared Hamiltonians, and independently replays each selected certificate. Its tests reject Hamiltonian mismatches and tampered bounds. Candidate selection uses certified values at the target size; it does not refit the amplitude parameters there, and all candidate evaluation work is recorded.

| Scenario | Previous M64 width | New combined M64 width | Reduction factor |
|---|---:|---:|---:|
| Repulsive t=V=1 | 12.29928 | 3.86842 | 3.18x |
| Varying bonds | 3.33939 | 1.03584 | 3.22x |
| Varying fields | 12.32743 | 3.73426 | 3.30x |

These values describe the fitted-amplitude route, before applying the exact algebraic solution to the repulsive case. All widths use the abstract hopping energy unit, not hartree.

| Scenario | New combined width, M12 | M16 | M32 | M64 |
|---|---:|---:|---:|---:|
| Repulsive | 0.63337 | 0.88508 | 1.89387 | 3.86842 |
| Varying bonds | 0.17369 | 0.25648 | 0.50745 | 1.03584 |
| Varying fields | 0.62786 | 0.96683 | 1.85153 | 3.73426 |

The intervals still grow with size at these fixed ranges. This is finite evidence about the searched family, not a proof that every finite-range amplitude must fail or that no compact certificate exists. Indeed, the exact algebraic certificate solves the repulsive control that this amplitude search leaves loose.

**The main remaining error is in the lower bound.**

Nine [independent sparse CAR references](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/ed_reference.json), at M=8,12,16 in each scenario, passed numerical eigenvector residual checks below 2.2e-12. They enumerate the small fixed-N space only for validation and were not used to fit parameters. Every combined certificate contains its numerical reference within the stated comparison tolerance.

At M=16, the fraction of interval width attributable to the lower bound's distance from the reference is 69.7% for the repulsive case, 90.7% for varying bonds, and 76.4% for varying fields. Energy-only fitting generally gives a good upper bound while leaving a poor lower one. The [combined analysis](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/range_analysis/summary.json) records both errors separately and replays the selected exact certificates for all 60 combinations of case, maximum range, and size. This analysis took 71.18 seconds locally with three workers; its replay cost is additional to the search timings above.

Fresh stdlib-only replay accepted **324 finite-range certificates**: 144 training candidates and 180 transfer witnesses. The first [replay receipt](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/finite_range_campaign/replay.json) and second [replay receipt](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/finite_range_lower/lower_replay.json) distinguish replay from independent small-matrix testing. Ten focused tests cover the finite-range extension, objective/combination logic, and algebraic certificate.

This round supports adding **local algebraic factorization** to the structural candidates. The concrete next gate is finding related identities or controlled residual certificates away from the exactly matched family. The current recognizer correctly rejects such perturbations; it provides no claim about general molecular chemistry, FeMoco, or finite-temperature superconductivity.

Reproduce the structural tests and a fresh ladder from `/Users/aidenlippert/Documents/Spectra`:

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m unittest research.side_routes_20260913.test_supersymmetric_chain -v
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.side_routes_20260913.supersymmetric_chain --out /tmp/spectra-supersymmetry-ladder.json
```

The [continuation manifest](/Users/aidenlippert/Documents/Spectra/results/side_routes_20260913/range_manifest.json) records delivered source and result hashes. The main session's solver and goal files were not modified by this side investigation.
