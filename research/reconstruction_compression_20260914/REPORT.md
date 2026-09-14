# Reconstruction-preserving proof compression — completed bounded pass

The main target remains open. The strongest independently constructed compressed intervals are **2.741879 mHa on H6** and **12.153150 mHa on corrected H8**. Both exceed 1.6 mHa. The preserved 0.158560 and 1.176860 mHa milestones and the CH2 model result remain unchanged.

The useful accepting-path result is a rigorous interval MPS verifier: H8 takes **3.957 s**, compared with **33.351 s** for the integer checker, an observed **8.43×** ratio on the same rounded state and declared endpoint. The endpoint allowance is **0.00005067 mHa**. Exact equality remains ambiguous and requires refusal or integer fallback.

## Complete accepting results

| Case | Preserved exact interval (mHa) | Fast-upper + preserved-lower interval (mHa) | New complete accepting time (s) |
|---|---:|---:|---:|
| H6 | 0.158560146 | 0.158617841 | 6.056 |
| H8 | 1.176859796 | 1.176910461 | 40.515 |

These new complete bundles recheck the inherited full-cubic lower proof. They establish a cheaper verifier with a charged, tiny upper allowance; they do not establish a cheaper lower constructor. Both upper paths were compared with the same endpoint ceil(2^24 U)/2^24. The main lower searches used the original exact U throughout.

## Full lower-plus-frozen-upper comparisons

Every row below has an independently replayed rational lower bound. “Passed” in a process receipt only means the program exited successfully; the interval target is assessed separately. Search time includes numerical construction/solve/export and excludes the shared preparation listed below. Rank is a cap per odd block; quadratic repair blocks stay full.

| Case / rule | Rank | Gram entries | Certified interval (mHa) | Search/export (s) | Accepting lower replay (s) | Peak solver MiB |
|---|---:|---:|---:|---:|---:|---:|
| H6 / nested cuts, paired | 8 | 6,410 | 11.200718 | 4.674 | 4.299 | 220.7 |
| H6 / MPS, paired | 32 | 14,244 | 2.741879 | 68.185 | 5.307 | 726.4 |
| H6 / MPS, paired; exact elimination; indirect solve failed numerically | 32 | 14,244 | 80368.842881 | 76.405 | 4.274 | 668.7 |
| H6 / MPS, paired; exact elimination | 32 | 14,244 | 2.774471 | 68.297 | 4.454 | 749.7 |
| H6 / MPS, paired | 8 | 6,496 | 8.538141 | 32.682 | 4.252 | 170.1 |
| H6 / simple, paired | 8 | 6,496 | 13.109728 | 2.739 | 0.226 | 156.6 |
| H6 / teacher, unpaired | 8 | 7,268 | 13.109958 | 20.134 | 5.177 | 676.0 |
| H6 / MPS, unpaired | 8 | 7,268 | 13.110076 | 47.818 | 4.888 | 608.7 |
| H8 / nested cuts, paired | 8 | 18,944 | 23.428920 | 13.979 | 24.376 | 887.8 |
| H8 / MPS, paired; exact elimination | 16 | 21,248 | 14.304603 | 28.712 | 26.520 | 704.2 |
| H8 / MPS, paired | 24 | 25,088 | 12.153150 | 80.236 | 29.184 | 1145.2 |
| H8 / MPS, paired | 8 | 18,944 | 20.927569 | 20.858 | 27.120 | 421.9 |
| H8 / simple, paired | 8 | 18,944 | 24.166273 | 10.387 | 0.793 | 293.6 |
| H8 / teacher, unpaired | 8 | 19,776 | 24.165545 | 64.419 | 26.726 | 1692.4 |
| H8 / MPS, unpaired | 4 | 18,624 | 24.165336 | 41.314 | 3.085 | 463.1 |

Shared Hamiltonian-map and MPS-moment preparation costs 7.447 internal seconds for H6 and 47.874 for H8, plus interpreter/import overhead in the run ledger. Every standalone candidate must pay that preparation. The initial H6 quotient-frame trial, replaced to match the teacher convention, is also charged in the campaign total.

The two useful rank-32 H6 rows use the existing wedge spectral residual checker; other rows use the full CAR coefficient remainder. Their separate spectral-proposal costs remain in the ledger, and the accepting column includes the selected checker’s complete lower replay. The result index includes proof bytes and exact endpoints. No omitted degree-six terms or floating-point feasibility claims substitute for acceptance.

The best searches use approximately 15.3× fewer Gram entries on H6 and 48.1× fewer on H8 than their uncompressed dictionaries. That does not give a matched-accuracy speedup. In particular, H8’s rank-24 projected coefficient map has 3,443,456 nonzeros and exceeded the memory target.

At rank eight, MPS-guided paired directions improve the interval relative to simple paired directions on both cases. The cut maps improve less. The tiny unpaired subspaces return almost the same bounds as the quadratic/simple control. This supports the relevance of state information, while showing that the tested small spaces lose too much reconstruction capability. It does not prove that larger or differently structured subspaces must fail.

## Exact equality reduction and the cut experiment

The paired anticommutator construction cancels degree-six terms before dense projection. An exact modular full-column-rank calculation proves that the quartic sector multipliers must then be zero in an exact reconstruction. H6’s active equations drop from 4,859 to 499 and its multiplier dimension from 499 to 25; H8 drops from 28,461 to 1,525 equations and from 1,525 to 41 multipliers. Equivalence applies to the paired model’s zero-remainder feasible set, not every residual-corrected certificate. It did not resolve its accuracy or convergence limitations.

The cut experiment constructs nested fixed maps at two orbital partitions, with a per-node dimension cap of two and an independent root rank. All root positivity constraints share global CAR moments and the full molecular Hamiltonian. Unknown states are never restricted to the MPS. However, this prototype still constructs full dictionary moments/maps and exports expanded factors. It is a finite operator-map restriction, not an efficient implementation of a molecular RDM consistency hierarchy. That higher-upside part remains unresolved.

## Frozen transfer

The generation and search rules were hashed before a new H6 geometry at **2.07 Å** was generated. The independent bond-48 MPS, paired rank-32 lower, exact equality reduction and 65-second solver budget produced a complete interval of **4.124284 mHa**. The target was not met. The frozen transfer uses coefficient-L1 acceptance, which gave 4.011289 mHa for the corresponding original-fixture construction; it does not add the later spectral-residual improvement. No FCI trial state or full-cubic teacher solve was used on this geometry. The exact interval concerns the rational electronic Hamiltonian; nuclear/geometry/model errors are separate.

## Accounting, validation, and preservation

1076.784 seconds of bounded process wall time are recorded against the initial 2,400-second campaign budget, including failures, preparation, tests and replay. Independent jobs sometimes overlapped, so this sum is not calendar elapsed time. Editing, file inspection, final inventory and inherited discovery are outside that sum. The ledger contains 63 completed process receipts.

The teacher rows also inherit the old full-cubic discovery (259.309 s for H6 and 721.563 s for H8). Frozen MPS construction is another inherited dependency: the retained H8 warm path took 414.36 recorded seconds. These historical figures do not account for every earlier failed attempt. None is silently counted as ordinary Hamiltonian input or free preprocessing.

Two numerical builds exceeded the 1 GiB working target: h8_low_rank24 (1.118 GiB), h8_teacher_rank8 (1.653 GiB). Their costs are retained; subsequent smaller cases and the frozen transfer stayed below that target. The indirect linear solver trial failed to converge to useful positivity and is retained as an 80,368 mHa failure, not a successful compression.

Twelve focused tests pass. They cover exact CAR reconstruction and degree-six retention, adjoint pairing, numerical moments against an explicit tiny oracle, positive congruences, modular rank and its refusal case, MPS/sector/Hamiltonian mutations, rigorous primitive enclosures, bad upper endpoints, and equality fallback. Initial implementation/test and report-assembly failures remain in the prospective logs.

All **9,864 inherited files** match their starting hashes. The fresh rule hashes are unchanged. No external resources were provisioned or published, and the CH2 model result was preserved.

The result supports keeping the faster rigorous upper checker. The lower-bound research still needs a representation that carries enough cross terms and constructs the coefficient map compactly. Increasing dense reduced rank already grows memory before H8 approaches the accuracy target. No exact energy-family obstruction was proved in this pass.

[Mathematical derivation and limitations](/Users/aidenlippert/Documents/Spectra/research/reconstruction_compression_20260914/MATHEMATICS.md) · [Exact result index](/Users/aidenlippert/Documents/Spectra/results/reconstruction_compression_20260914/final_result.json) · [Accounting audit](/Users/aidenlippert/Documents/Spectra/results/reconstruction_compression_20260914/audit.json)
