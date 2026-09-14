# V9 trajectory wide Galerkin baseline

This preflight implements a conventional exact adaptive Galerkin proposer. A
finite Pauli basis is seeded by the initial support. For each fixed basis it
forms the projected generator by applying every required column and propagates
the projected Taylor recurrence exactly. At every time coefficient it computes
the full omitted residual `(I-P)G c_k`, including collisions and signs.

Labels are ranked by accumulated power integration weight
`sum_k |r_k(p)| T^(k+1)/(k+1)`. A deterministic batch is inserted and the
projected recurrence is restarted. The expansion count is bounded at 16 and
the polynomial degree at 24. Failed rounds, residual scoring, generator cache
work, and the final fresh v7 checker replay remain in the report.

The exported candidate is an ordinary exact `Piece`, so the unchanged checker
recomputes jumps and every residual. A result is either independently
certified at `1/1000` for the requested horizon or reported as a bounded
refusal. This is a baseline construction, not learner training and not a
held out evaluation. No global magnitude drop rule or nilpotent high gain
contractivity claim is used.
