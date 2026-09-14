# Joint pure-coherence energy certificates and a transfer limitation

ENERGYv19 and FAMILYv15 integrate the two independent pure-coherence constraints identified after complete diagonal spin-word closure. Reoptimizing older coefficients jointly with the new terms produces stronger exact energy lower bounds. Both new projector moments close, but full quantum overlap still requires a separate test. The frozen held-out comparison reveals a transfer limitation.

## Accepted matched results

| Case | Periodic lower/site | Million-site open lower/site | Enlarged family ceiling | Family gap |
|---|---:|---:|---:|---:|---:|
| W_zero | -0.64290658388107558 | -0.64290908388107559 | -0.64286167896845148 | 4.49049126241e-05 |
| W_plus_1 | -0.66053609624452647 | -0.66054059624452655 | -0.66036907503281894 | 0.000167021211708 |

The family ceilings limit attainable lower certificates in this specified relaxation. They are not physical ground-energy upper bounds. Physical uppers remain -0.6106763470511881 for W=0 and -0.6184244823693281 for W=1.

- W_zero: the lower improves by 5.53523367336e-06/site over the preceding certificate. Signed separation from the preceding full-family ceiling is -5.19779969776e-05; strict full-family separation proved: False. The selected enlarged ceiling uses 183 positive sources, maximum weight length 1859 characters. Selected evidence: `results/marginal_graded_hubbard8/pure_coherence/W_zero/refined`.
- W_plus_1: the lower improves by 2.65871180638e-05/site over the preceding certificate. Signed separation from the preceding full-family ceiling is -0.000353425424981; strict full-family separation proved: False. The selected enlarged ceiling uses 194 positive sources, maximum weight length 1864 characters. Selected evidence: `results/marginal_graded_hubbard8/pure_coherence/W_plus_1/scaled`.

The preceding ENERGYv18 and FAMILYv14 proofs were freshly replayed under the current implementation. A stronger certificate alone does not prove the new directions outperform every reoptimized recipe in the preceding family.

## What joint reoptimization establishes

The earlier exact one-source and three-source caps allowed arbitrary real coefficients on the two new terms while keeping all older fields frozen. Those physical caps were reconstructed again here under current code. Both jointly optimized lower bounds exceed those caps:

| Case | Joint lower minus frozen-old-recipe ceiling | Removing only new terms from current fixed recipe loses |
|---|---:|---:|
| W_zero | 5.49754196733e-06 | 0.00113616 |
| W_plus_1 | 2.64143211116e-05 | 0.0010779 |

These comparisons have different scopes. Exceeding the frozen-old-recipe ceiling proves that allowing older coefficients to respond mattered. The matched ablation measures removal of the new terms from the resulting fixed recipe. Neither comparison establishes separation from the reoptimized preceding full family. Certificate byte identity binds the comparisons to the selected polished energy when a different mixture directory supplies the best ceiling.

## Exact closure and residual quantum consistency

The new sources are |346>+|409> and |314>+|614>. Removing the diagonal from their symmetry-averaged reflection-odd projector telescopes leaves five-site operators with 8 and 16 entries and six-site telescopes with 60 and 120 entries. Their norms are bounded by 1/2. The first carries four-bit occupation-dependent spin exchange, the second six-bit coherence. Prior exact annihilator functionals establish independence from all older affine terms.

Each selected positive mixture has both new moments exactly zero, all 120 spin-word moments zero, and the entire preceding constraint hierarchy. An independent complete five-site density-matrix reconstruction also checks the two original positive-projector differences, rather than relying only on the production moment fields. Full diagonal prefix/suffix laws still agree and retain their stationary classical order-five extension.

- W_zero: complete overlap has 0 diagonal mismatches and 29340 nonzero upper-triangle differences. Quantum extension of this particular symmetry-averaged mixture refuted: True. Positive-projector vector {'358': 1, '409': 1} detects mismatch 0.000162442864213.
- W_plus_1: complete overlap has 0 diagonal mismatches and 30108 nonzero upper-triangle differences. Quantum extension of this particular symmetry-averaged mixture refuted: True. Positive-projector vector {'103': 1, '358': 1} detects mismatch -0.000137849657358.

Local positivity remains intact. A residual overlap refutes extension of the particular mixture, not all mixtures at its energy or the accepted ceiling. Classical occupation-law consistency does not extend the quantum coherences or impose fixed global particle number. The complete symmetry-compatible real test space remains 3960-dimensional: 120 diagonal and 3840 offdiagonal directions.

## Frozen transfer

The initial W=0 recipe before polishing was frozen and transferred to U=5, t=1, V=1/4, W=-1/5. All projector sources, penalties and auxiliary coefficients were preserved; the recorded physical profiles were rescaled or shifted and only the scalar threshold was recomputed. Geometry and filling were unchanged.

| Million-site open-chain recipe | Lower/site |
|---|---:|
| New recipe with pure-coherence terms | -0.5223150350391923 |
| Only those two terms removed | -0.52285071503919223 |
| Previous frozen spin-word recipe | -0.52200193331684708 |

The signed contribution of the new terms is 0.00053568/site. The signed change of the complete new recipe relative to the previous frozen recipe is -0.000313101722345/site: a negative value means a worse lower bound. Thus helpful new terms do not establish robustness of the full jointly adapted recipe. All three energy replays and exact frozen-field comparisons pass. The physical upper remains -0.48856168029895469. No broader transfer claim is made.

## Implementation and validation

The new production module reuses the existing signed projector builder and removes diagonal entries, which commutes with its checked signed permutations and embeddings. Both new operators enter every local matrix before PSD replay. ENERGYv19 requires the full previous hierarchy; older versions reject the new field. FAMILYv15 requires both new zero moments and rejects fixed-coefficient fields. Its conservative source cap is 209, while the numerical system has 200 independent rows. The 4096-character weight limit, 4096-state coverage, 94 blocks and maximum local PSD dimension 200 are unchanged.

The joint search varies 199 coefficients: 120 full spin-word directions, 75 offdiagonal directions, two hopping profiles and two penalties. Each run is bounded by 500 full-spectrum evaluations, with gradient/Hessian checks and fresh physical reconstruction. Initial W0/W1 searches used 460/500 evaluations; polishing used 461/500. Family pricing uses bounded physical integer vectors, then exact fraction-free reconstruction and independent physical replay. Initial pricing has 40 two-eigenvector rounds and refinement has up to 80 one-eigenvector rounds. Iteration limits and remaining negative reduced eigenvalues do not prove convergence. Timing records are not controlled performance benchmarks.

- Rejected exact basis in `results/marginal_graded_hubbard8/pure_coherence/W_plus_1/refined/diagonal_family_limit_diagnostic.json`: Exact proposed basis is inconsistent. This attempt contributes no accepted ceiling.
- Zero-pricing export retry in `results/marginal_graded_hubbard8/pure_coherence/W_plus_1/scaled/diagonal_family_limit_diagnostic.json` scaled the numerical constraint rows by 1000 and used highs-ds_without_presolve. It produced an exact proposal; independent physical replay and the original 4096-character weight gate decide acceptance. No additional pricing or gate relaxation was used.

All 74 focused integration tests, 3 fraction-free tests, and the full 1128 tests plus 102 subtests pass. The existing calibration return-value warning remains. A copied fraction-free test initially retained the old dimension boundary; it was corrected to exercise 200 accepted rows and 201 refused rows, then passed. The failed test log is preserved. No production source changed after full-suite collection. The collector audits 38 current receipts and 1591 proof hash entries, plus 291 historical receipts against unchanged sources or preserved snapshots. Both previous provenance manifests are accounted for. No agents, GPU or paid resources were used.

Full-family numerical attainment, general quantum representability, generic molecular/long-range/higher-dimensional transfer, and scalability at requested accuracy remain unproved. The goal remains active.
