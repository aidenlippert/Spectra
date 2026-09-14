# Coherent-projector energy certificates and a transfer limitation

ENERGYv17 and FAMILYv13 integrate the two exact five-site projector-overlap constraints found in the preceding diagnostic. Both matched-model lower bounds improve after fresh exact PSD replay. The new terms help the matched fixed recipes, but their frozen contribution is harmful at the held-out couplings. Neither enlarged numerical limit is resolved.

## Accepted matched results

| Case | Periodic lower/site | Million-site open lower/site | Family ceiling | Remaining family gap |
|---|---:|---:|---:|---:|
| W_zero | -0.64292565100821752 | -0.64292815100821743 | -0.64289051857768398 | 3.51324305335e-05 |
| W_plus_1 | -0.66060626287135837 | -0.66061076287135845 | -0.66039728039978174 | 0.000208982471577 |

Family ceilings bound attainable LOWER certificates in the specified family; they are not physical ground-energy uppers. Physical upper bounds remain -0.6106763470511881 (W=0) and -0.6184244823693281 (W=1). Each selected family mixture has 139 positive physical sources and satisfies both new moments and the preceding hierarchy exactly. Independent projector contractions also vanish.

- W_zero: lower improvement 6.1237590234e-06/site; signed separation from the preceding whole-family ceiling -2.10050762005e-05; strict whole-family separation proved: False. Removing only the two new terms from this fixed recipe loses 0.0005827/site.
- W_plus_1: lower improvement 2.57706019341e-05/site; signed separation from the preceding whole-family ceiling -0.00057164758007; strict whole-family separation proved: False. Removing only the two new terms from this fixed recipe loses 0.00020318/site.

The larger matched ablation effects concern jointly adapted fixed coefficients. They do not imply the same improvement over a reoptimized older family. The preceding energy and family certificates were freshly replayed under the current implementation for these comparisons.

## Frozen transfer and signed ablation

The INITIAL W=0 recipe, before polishing, was frozen and transferred to U=5, t=1, V=1/4, W=-1/5. Projector sources, penalties and all auxiliary coefficients were preserved. Physical profiles were rescaled or shifted by the recorded recipe; only the scalar spectral threshold was recomputed. Geometry, filling and interaction range stayed fixed.

| Million-site open-chain recipe | Lower/site |
|---|---:|
| With both coherent-projector terms | -0.52245353505179537 |
| Only those two terms removed | -0.52164979505179543 |
| Previous frozen three-spectator recipe | -0.52258411426757478 |

The signed contribution of the frozen new terms is **-0.00080374/site**. Removing them improves the bound: their matched benefit does not transfer with these coefficients. The full new frozen recipe nevertheless improves on the previous frozen recipe by 0.000130579215779/site. Its other jointly adapted coefficients matter. The physical upper is -0.4885616802989547. All three energy replays and exact frozen-field comparisons are accepted.

An initial use of the matched-target upper-state driver correctly refused all three held-out cases because that target has no matching stored filter recipe. Those refusal logs are retained. The dedicated transfer driver then evaluated the fixed W=0 physical filter at the actual held-out Hamiltonian and passed exact lower/upper replay; the source-selection gate was not weakened.

## What remains outside the implemented overlap constraints

An exact signed-permutation orbit census finds 32,264 real Hermitian five-site matrix directions preserving both spin numbers. Imposing particle-hole evenness, spin-flip evenness and reflection oddness leaves **3,960 directions: 120 diagonal and 3,840 offdiagonal**. Every matrix unit belongs to one checked orbit; disjoint consistent orbits give independent basis vectors and sign-conflicted orbits give none. This is an operator-space dimension statement, not a performance or general representability result.

- W_zero: stationary extension of the new selected symmetry-averaged mixture refuted by complete overlap replay: True. Nonzero upper-triangle differences: 30760. Positive-projector witness {'102': 1, '153': 1} has exact mismatch summarized by -0.000245148573397.
  The spin-resolved diagonal word laws have exact total variation 0.000832856009136, detected by a positive diagonal projector onto 308 states. Every coarser charge-word overlap is exactly zero. Thus even full diagonal spin-word consistency remains open, despite charge-law closure; the signed indicator telescope has norm at most one.
- W_plus_1: stationary extension of the new selected symmetry-averaged mixture refuted by complete overlap replay: True. Nonzero upper-triangle differences: 31000. Positive-projector witness {'103': 1, '358': 1} has exact mismatch summarized by -0.00021321913098.
  The spin-resolved diagonal word laws have exact total variation 0.00194154048642, detected by a positive diagonal projector onto 332 states. Every coarser charge-word overlap is exactly zero. Thus even full diagonal spin-word consistency remains open, despite charge-law closure; the signed indicator telescope has norm at most one.

Closing the two added moments does not close full quantum consistency. A failed overlap test refutes the extension of that particular mixture, not all mixtures at its energy or the accepted family ceiling. The census motivates a broader operator-space treatment rather than assuming two more scalar tests establish representability.

## Implementation and validation

The production implementation fixes two canonical projector sources. It checks signed orbit closure, Hermiticity, spin numbers, reflection/particle-hole/spin-flip invariance and fermionic telescoping on the full 4096-state space. Both operators enter the local matrix before exact PSD acceptance. Local coverage remains 94 blocks, with maximum local PSD dimension 200. Older energy versions reject the new field. FAMILYv13 requires the full preceding hierarchy and two exact zero moments, allows at most 139 sources, and forbids fixed-coefficient fields. Older source caps and the 4096-character rational-weight limit are preserved.

Numerical proposals optimize 138 coefficients with at most 500 full-spectrum evaluations per run. Initial W0/W1 runs used 412/500 evaluations; polish runs used 411/500. Fresh physical reconstruction and derivative checks accompany the nonaccepting proposals. Initial family pricing ran 40 rounds with two eigenvectors per block; refinement ran 80 rounds with one. Exact bounded 139-row fraction-free reconstruction is separate from numerical selection. Remaining negative reduced eigenvalues and iteration/evaluation limits do not establish convergence.

All **81 focused integration tests**, **3 focused fraction-free tests**, and the full **1039 tests plus 102 subtests** pass. The existing calibration return-value warning remains. No production source changed after full-suite collection. 29 current accepted receipts and 1119 source-hash entries were checked; 225 historical receipts are audited against unchanged sources or preserved snapshots. Failed precondition attempts and untrusted proposals are retained. No GPU, paid resources or new agents were used.

General quantum representability, exact numerical-limit attainment, generic molecular or long-range/higher-dimensional transfer, and scalability at requested accuracy remain unproved. The goal stays active.
