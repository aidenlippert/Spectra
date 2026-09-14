# V8 time reduction probe

This development experiment tests fixed-grid time stepping as a conventional
headroom candidate.  A segment of duration `h` uses the exact sparse Taylor
recurrence

\[
 P_j = G^j O_0/j!,\qquad R(t)=P'_m(t)-G P_m(t),
\]

and the endpoint polynomial becomes the initial observable for the next
segment.  The unchanged v7 checker receives all segments together.  Thus each
endpoint mismatch is an explicit `jump:j` record, and the claimed error is the
sum of every jump and integrated residual over the complete horizon.  There is
no resetting of the error budget at segment boundaries.

The algebraic preflight does not give a general advantage.  Splitting a Taylor
polynomial reduces the local factor from `T^(m+1)/(m+1)` to `h^(m+1)/(m+1)`,
but repeats support generation and exact rational arithmetic. With the next
segment initialized to the previous polynomial endpoint, interface jumps are
exactly zero. The previous approximation error is already accumulated through
the earlier residual integral; it must not be added again as a jump.
The method remains ordinary time stepping, not a new mathematical operation.

## Probe and accounting

`experiments/v8_time_reduction_probe.py` uses the fixed v7 task (XXZ and mixed
families, widths 3/4/6, gamma `0`, `1/5`, `2`, horizons `1/5`, `1/2`,
tolerance `1/1000`).  It tries segment counts 1, 2, and 4 independently and
orders 4, 8, 12, 16, 24 for each arm, with the v7 caps (support 512, order
24).  Every candidate constructs witnesses with power and Bernstein integration
and all three existing norm groupings, then performs a fresh independent check.
Each case also gets a fresh adaptive-Taylor baseline. Failed search work is
retained per arm; a structural budget refusal closes only that fixed-grid arm.
Run from the repository root with:

```text
/Users/aidenlippert/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 -m experiments.v8_time_reduction_probe
```

The machine-readable output is
[`results/v8/time_reduction_probe.json`](../../results/v8/time_reduction_probe.json).

## Result

Across the 36 fixed physical cases and their three independently searched
segment arms, there were 172 attempted order settings: 72 certified, 64 over
tolerance, and 36 refused. There are 24 certified arms for each segment count
(12 cases per family); the mean elapsed attempt time for certified arms was
about 0.069 s for one segment, 0.112 s for two, and 0.211 s for four on this
run. Segmentation therefore increased measured complete attempt cost while
offering no new feasibility case. Width-6 mixed cases refused at the existing
support budgets, and width-6 XXZ cases only reached order 4 before refusal.
This targeted diagnostic does not establish a failure theorem or complete-cost
headroom, and supplies no reason to replace the strongest v7 conventional arms.

The output includes construction, witness, checker, support, and generator
counters for certified attempts, while refusal rows retain the work counters
at failure.  All phases therefore count toward the development comparison.

## Checker conclusion

No new functional checker is needed.  The existing piecewise checker already
verifies total horizon equality, every jump, every residual record, exact norm
witness coverage, and one accumulated bound.  A future time-stepping method
would need only this same interface; any proposed optimized checker would need
to preserve explicit jump records or provide an independent proof that their
errors are included in the final-horizon bound.

This is a conventional baseline measurement, not a novelty claim, learner
result, or held-out evaluation.
