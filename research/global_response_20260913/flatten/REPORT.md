# Flattened two-level response schedule

This pass implements a semantics-preserving flattened execution wrapper for
the existing nested `recursive.retained_action`. It memoizes only identical
base-H right-hand sides. Projections, polynomial recurrences, and the order of
the two congruences are unchanged. The wrapper therefore tests whether the
current nested schedule contains common subexpressions; it does not introduce
a new bound or imply an enumeration-free implementation.

For first response degree `k1` and second degree `k2`, the existing schedule
requires `2*k1+1` parent-H actions per first retained action. The outer
polynomial has `(k2-1)*(2*k1+1)` actions on its right-hand side, and a complete
second retained action has `(2*k2+1)*(2*k1+1)`. With the current H6 parameters
(`k1=4`, diagnostic `k2=3`) this is 63 actions. For the production values
(`k1=26`, `k2=119`) it is 12,667 actions.

The exact rational tests use a noncommuting four-dimensional symmetric block.
The flattened and recursive schedules agree exactly. The zero-vector case
shares repeated right-hand sides, demonstrating that the cache is real rather
than a relabeling of the count.

The matched production-order H6 diagnostic uses the largest 200-state conserved
block generated from the 924-state fixture, with the physical Q1/R2/P2
occupation masks, `k1=26`, `k2=119`, and the exact `eta1=5.207645551389818e-7`
shift included in the inner operator. The final `eta2` shift is also applied.
The formal and measured counts are exactly 12,667 base-H right-hand sides. On
an unscaled float64, one-BLAS-thread run, both schedules produced maximum
absolute difference `0.0`; the nested run took 12.834 s and the memoized run
13.641 s. All 12,667 right-hand sides were distinct, so generic caching saved
no work. The block was explicitly enumerated for this diagnostic (924 labels
and 132,004 matrix entries); it is evidence about schedule cost only.

An earlier scaled-H diagnostic is superseded and must not be used as evidence;
the receipt now points to this unscaled physical-mask run.

The algebraic identity `F_k(d)=c_k p_{2k}(d)`, with
`c_k=T_{2k}(z0)/(T_{2k}(z0)+1)`, follows from
`T_k(z)^2=(T_{2k}(z)+1)/2`. It can reduce scalar recurrence bookkeeping, but it
does not reduce the leading H-action degree: the distinct outer right-hand
sides still require the same 12,667 oracle applications. A 10x reduction
therefore requires a genuinely new operator oracle or exploitable structure in
the physical vectors.
