# V7 grouping probe

This is a preliminary norm-grouping experiment, not evidence of full dynamics or autonomous discovery headroom.

The current probe constructs exact residual coefficients `G^m O / m!` by repeated `experiments.pauli.commutator_i` application for a local XX+ZZ+Z chain, with no truncation. It compares first fit, coefficient-weighted greedy packing, and singleton l1. It records separate exact norm intervals, executed anticommutation comparisons during construction, and independent checker pair checks. Every partition is independently verified for coverage, pairwise anticommutation, and rational square-root bounds.

The earlier artifact is preserved as `results/v7/grouping_probe_invalidated.json`. It must not be used as evidence: it summed Taylor coefficients at `t=1` rather than measuring a residual, rounded complex generator coefficients, silently sliced terms, used an invented construction proxy, combined unlike costs with an arbitrary weight, and never invoked independent checking.

Run `python experiments/v7_grouping_probe.py`; output is `results/v7/grouping_probe.json`. Results support only bounded exploratory norm grouping within this declared chain family.
