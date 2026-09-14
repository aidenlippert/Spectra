# Optional adaptive-pricing modes audit

The full-row path with pruning passes an independent M4 control. With one
requested iteration it produces two LP checkpoints, preserves a replayed
best certificate, prunes inactive atoms, and records cumulative additions and
prunes. The hopping toy's lower remains below the exact variational
upper -7/50 from symmetric trial amplitudes (1,3,1). Pruning is therefore a valid operational mode in this bounded test;
the atom cap is still a hard retained-pool cap while cumulative additions and
re-pricing are reported separately.

The corrected `half_rows=True` path now filters ideal terms to the retained
rows. On the M4 hopping control, the initial full- and half-row numeric LP
objectives agree to better than `1e-7`, and the half-row exported certificate
replays below the exact variational upper -7/50. The row weights are 1 on
diagonal entries and 2 off diagonal, while pricing symmetrizes the reduced
dual moment matrix. This is a valid finite control, though it is not yet a
measured speedup claim.

The control also confirms that all 32 M4 body-2 multiplier basis elements are
Hermitian. This supports the reduced-row treatment for the current basis; a
future basis extension should retain an explicit Hermiticity check.

The full-row initial checkpoint and pruned replay remain valid finite
certificates after exact `verify()`. As before, pruning and greedy pricing do
not prove global optimality; all pricing maps and cumulative candidate costs
must remain in the scaling report.

## Square width-4 diagnostic

The later square ablation separates the effects more clearly. Full rows without
pruning reached `-3.47147133` in about 71.97 s; half rows without pruning
reached `-3.47222092` in about 69.02 s. Thus half-row reduction alone was
comparable on this run. Full rows with pruning stopped at `-4.45900878` in
about 7.88 s after pruning 5,251 atoms, while half rows with pruning reached
`-4.38075138` in about 23.4 s after pruning 12,596 atoms. These are exact
replayed incumbent certificates, so the difference is not an exact-validity
failure. Aggressive pruning is the leading measured suspect.

Half-row reduction can still return a different dual representative because LP
duals are nonunique, so dual degeneracy remains plausible but is not proven by
these measurements. The implementation does preserve every current LP column
with weight above `1e-9`, but it does not preserve every atom or constraint from
an older exact-best checkpoint; released keys can be repriced. Historical pool
churn is therefore the first mechanism to address.

The recommended next change is to disable aggressive deletion or retain
historical pricing constraints/atoms in a cut-management pool, while keeping
the incumbent exact certificate as a fallback. Only after that ablation should
dual canonicalization or a secondary minimum-norm objective be tested. This is
a heuristic-search change and leaves exact certificate validity intact.
