# Full spin-word constraints: energy certificates and classical closure

ENERGYv18 and FAMILYv14 implement the complete 120-dimensional diagonal five-site space that is even under particle-hole and spin flip and odd under reflection. Both matched energy bounds improve after exact PSD replay. The selected positive local mixtures now have identical full spin-resolved five-site occupation laws, with an explicit stationary classical extension. Full quantum extension remains a separate requirement.

## Accepted matched results

| Case | Periodic lower/site | Million-site open lower/site | Family ceiling | Family gap |
|---|---:|---:|---:|---:|
| W_zero | -0.64291211911474899 | -0.64291461911474901 | -0.64285460588409793 | 5.7513230651e-05 |
| W_plus_1 | -0.66056268336259027 | -0.66056718336259035 | -0.66018267081954596 | 0.000380012543044 |

Family ceilings bound attainable lower certificates in this specified family; they are not physical ground-energy upper bounds. Physical upper bounds remain -0.6106763470511881 (W=0) and -0.6184244823693281 (W=1).

- W_zero: lower improvement 1.35318934685e-05/site; signed separation from the preceding whole-family ceiling -2.1600537065e-05; strict separation proved: False. Removing only the new spin-word field from the adapted fixed recipe loses 0.03243578/site. The selected ceiling uses 192 positive sources, with maximum rational weight length 2103 characters. Refinement lowers the initial ceiling by 5.39849675091e-05/site. Selected evidence: `results/marginal_graded_hubbard8/spin_word/W_zero/scaled`.
- W_plus_1: lower improvement 4.35795087681e-05/site; signed separation from the preceding whole-family ceiling -0.000165402962809; strict separation proved: False. Removing only the new spin-word field from the adapted fixed recipe loses 0.01586254/site. The selected ceiling uses 198 positive sources, with maximum rational weight length 2768 characters. No tighter refinement passed exact acceptance. Selected evidence: `results/marginal_graded_hubbard8/spin_word/W_plus_1/final`.

The preceding energy and family proofs were freshly replayed under current code. Matched ablations are bound to the selected energy by certificate byte identity. Their larger effects concern fixed, jointly adapted coefficients and do not establish improvement over a reoptimized older family.

## Exact diagonal closure and remaining quantum obstruction

The independent diagnostic reconstructs all 4096 six-site occupation probabilities from the positive source mixture and its eight symmetry images. Every one of the 1024 five-site prefix probabilities equals the corresponding suffix probability exactly. All 120 new family moments vanish, as do the preceding constraints.

The closure also has a general finite-window explanation. For any mixture averaged over these symmetries, the prefix-minus-suffix diagonal vector is particle-hole even, spin-flip even and reflection odd. It therefore lies in the complete 120-dimensional basis. The zero telescope moments make it orthogonal to every basis vector, so it must vanish. Disjoint signed orbits establish independence. This implication concerns the diagonal of the symmetry-averaged mixture; the independent integer reconstruction verifies it for each accepted instance.

For a five-word prefix u with mass q(u)>0, append symbol a with probability p6(u,a)/q(u). Prefix/suffix equality proves normalization and stationary inflow. Zero-mass prefixes may append empty deterministically. This constructs an order-five stationary classical process reproducing the full six-site occupation law over empty, up, down and double. It does not preserve the offdiagonal matrix elements of the supplied quantum mixture or impose a fixed global particle number.

- W_zero: 960 positive-mass Markov states and 2566 positive flows. Independent complete density-matrix replay finds 0 diagonal differences and 30136 nonzero upper-triangle differences. Stationary quantum extension of this particular symmetry-averaged mixture refuted: True. The normalized positive projector from vector {'346': 1, '409': 1} has expectation mismatch -0.000219649124863.
- W_plus_1: 1024 positive-mass Markov states and 3228 positive flows. Independent complete density-matrix replay finds 0 diagonal differences and 30324 nonzero upper-triangle differences. Stationary quantum extension of this particular symmetry-averaged mixture refuted: True. The normalized positive projector from vector {'314': 1, '614': 1} has expectation mismatch -0.000270921208986.

A nonzero offdiagonal overlap difference refutes extension of the particular mixture, not all mixtures at its energy or the family ceiling. Local positivity is intact. The exact symmetry census gives 3960 allowed real overlap directions: 120 diagonal and 3840 offdiagonal. Closing the diagonal part does not close the full space.

## Frozen transfer

The initial W=0 recipe before polishing was frozen and transferred to U=5, t=1, V=1/4, W=-1/5. All projector sources, penalties and auxiliary coefficients remain fixed; the recorded physical profiles are rescaled or shifted and only the scalar threshold is recomputed. Geometry, filling and interaction range are unchanged.

| Million-site open-chain recipe | Lower/site |
|---|---:|
| With new spin-word field | -0.52200193331684708 |
| Only that field removed | -0.53905959331684705 |
| Previous frozen coherent-projector recipe | -0.52245353505179537 |

The signed frozen contribution is 0.01705766/site; the full recipe improves on the previous frozen recipe by 0.000451601734948/site. The physical upper remains -0.4885616802989547. All three energy proofs and the exact frozen-field comparison pass. This single parameter transfer does not establish broader transfer robustness. The preceding stage's adverse two-term transfer remains historical evidence; this result does not erase it.

## Implementation, numerical limitations and validation

Canonical signed symmetry orbits partition all 1024 five-site spin words and produce exactly 120 independent diagonal directions. Independent full-group projection tests establish completeness and containment of the previous 52 charge directions and nine selected shapes. The new field enters the local matrix before PSD replay. Coverage remains 4096 states and 94 blocks, with maximum local PSD dimension 200. Older energy versions reject the new field. FAMILYv14 requires the preceding hierarchy and every new zero moment. Its conservative source cap is 207; the actual numerical system has 198 independent rows because the selected nine shapes are redundant. The 4096-character rational-weight limit is unchanged.

The numerical energy search has 197 coordinates and a 500 spectral-evaluation cap. Initial W0/W1 runs used 410/500 evaluations and polishing used 417/500. Family pricing used 40 rounds with two eigenvectors, then 80 rounds with one, preserving accepted physical vectors and determinant anchors. The first W0 family search failed with a HiGHS Unknown status. Its log is retained. A separate bounded fallback tries dual simplex and interior point without presolve, retains explicit candidate-feasibility checks, and never promotes numerical output to acceptance. Both first refinement exports failed exact consistency despite successful numerical LP status. A bounded zero-pricing retry scales all numerical constraint rows by 1000 and uses dual simplex without presolve; it retains the original rational equations and acceptance gates. Each selected mixture is reconstructed by bounded 198-row fraction-free arithmetic and replayed from physical vectors. The reported directory identifies the best accepted ceiling; failed proposals do not count as refinements. Negative remaining reduced eigenvalues and finite iteration limits do not prove convergence; recorded durations are not controlled benchmarks.

All 63 focused integration tests, 3 focused fraction-free tests, and the full 1074 tests plus 102 subtests pass. The existing calibration return-value warning remains. No production source changed after suite collection. The collector checks 32 current accepted receipts with 1292 source-hash entries and 254 historical receipts against unchanged sources or preserved snapshots. The preceding 481-file manifest is accounted for. Numerical attempts and failure logs remain in provenance. No agents, GPU or paid resources were used.

General quantum representability, exact numerical-limit attainment, generic molecular or long-range/higher-dimensional transfer, and scalability at requested accuracy remain unproved. The goal remains active.
