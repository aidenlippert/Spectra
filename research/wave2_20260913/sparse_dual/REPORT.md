# Wave 2: adaptive dual guided factor pricing

Discovery started from the original rational Hamiltonians only. No source
factors, FCI space, or upper state entered the optimization. The exact CAR
replay was run on every exported incumbent.

H4 fixture SHA-256: `8e64505deb4cc06a87f75e669b0f3073d15ebca72daa6d086e06eeea0b0f5120`.
The no-pruning run reached the 512 atom cap in one solve: 512 retained atoms,
5,772 complete pricing-map nonzeros, 849 coefficient rows, 0.031 s LP,
0.009 s export/replay, and certified lower `-8.40836393255014 Ha`.

The pruning comparison retained 120 atoms after 9 solves, added 640 and
pruned 520 candidates cumulatively, and gave lower `-7.06726696579605 Ha`.
That is a stronger lower bound than `-8.40836393255014 Ha` because larger
values are better lower bounds. It is not matched: no-pruning started at 512
atoms and made zero adaptive additions, while pruning opened slots. This was
aggressive inactive-atom pruning, not conservative incumbent preservation.

The matched larger-cap H4 run used width 4, batch 64, 20 iterations, and cap
2048. It started with 512 seed atoms, made 1,792 adaptive additions, pruned
none, reached 1,792 retained atoms after 21 LP rounds, and exact-replayed to
lower `-4.90788083670932 Ha`. Pricing-map nonzeros were 5,772; LP time was
3.281 s, pricing 0.525 s, export/replay 0.716 s, and total wall 4.748 s.

H6 fixture SHA-256: `ea50a4ab19a353bc89691bb344fcf8d610f74e1c2c76bf07c6869be2b24dbec8`.
The no-pruning run reached 512 atoms with 29,682 complete pricing-map
nonzeros, 4,501 coefficient rows, 0.049 s LP, 0.063 s export/replay, and
certified lower `-40.27557013159622 Ha`. This is a negative accuracy result;
the sparse family is not competitive at this budget.

The paired selected-refinement upper witnesses were not available in the
Wave 2 tree, so no upper interval claim is made here. They may be attached as
external state-assisted context later, with provenance kept separate from
H-only discovery. No asymptotic scaling or complete treewidth claim was tested.

Outputs:

- `results/wave2_20260913/sparse_dual/h4_noprune/`
- `results/wave2_20260913/sparse_dual/h4_conservative/`
- `results/wave2_20260913/sparse_dual/h6_noprune/`
- `results/wave2_20260913/sparse_dual/h4_width4_cap2048/`

Next promising parameter: H4 width `8` with no pruning and atom cap `2048`,
then compare exact lower and pricing cost against the width-4 run. The H6
result suggests improving atom ranking before scaling H6 caps.

## Width-0 full-pricing follow-up

The requested H4 width-0 run used cap 4096, batch 64, 60 iterations, no
pruning, and a 120-second wall budget. It started with 512 seed atoms, made
1,648 accepted additions, reached 2,160 retained atoms, completed 29 LP
rounds, and stopped on the wall budget. Exact replay gave lower
`-4.19353851257858 Ha` and residual L1 `0.66865951257858`. The run consumed
72.685 s LP, 0.504 s pricing, and 120.873 s wall time. It therefore improved
the width-4/cap-2048 lower (`-4.90788083670932`) in the expected direction
(larger lower values are stronger), while not reaching the requested cap.
The machine receipt records the original fixture path, sector, and canonical
H hash.
