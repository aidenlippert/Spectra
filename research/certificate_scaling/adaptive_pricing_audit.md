# Adaptive factor pricing audit

The current implementation was checked with three independent controls in
`test_adaptive_factor_pricing.py`. The (2\times2) oracle
`[[0,1/8],[1/8,1/4]]` produces a negative width-2 pricing direction. A small
H-only run exports a rational certificate that passes the production checker;
the replayed lower equals the receipt lower, source factors and upper witnesses
are both excluded, and the complete pricing map is nonempty. A body-2 ideal
run constructs a larger coefficient basis (more than the 47 rows in the
body-1 M4 control) and replays successfully with residual degree at most six.
All three controls pass.

The iterative snapshot is internally coherent on the recorded square run:
each history entry stores the LP lower and exact replay lower, while the final
receipt retains the best exact replayed certificate even though later rounds
may have more atoms. The receipt records the full pricing Gram dimensions,
map nonzeros, LP time, pricing time, export time, and the 2,828-word-pair seed
scan. It explicitly marks omitted-family optimality as unproved. A full
residual charge makes that omission irrelevant to lower-bound validity; it is
still relevant to discovery-scaling claims.

The exact moment maps use the production `gram_map`, including the dagger on
the left word. Candidate directions are quantized before being added to the
LP, and accepted factors are quantized again and passed through `verify()`.
The method therefore has a valid finite-certificate path. Numerical LP
objectives and dual pricing eigenvalues remain proposal diagnostics, not exact
proofs.

The current controls do not establish that the greedy pricing rule is globally
optimal, that the full candidate pool can be searched subquadratically, or that
the method transfers to larger active spaces. The recorded complete pricing
cost must remain part of every scaling comparison.
