# V7 query identifiability gate

For an exact linear parameter law `y = X theta`, a scalar query `r theta` is
identified by the observed design exactly when `r` is in the row space of
`X`. The implementation returns rational proof weights `w` satisfying
`X.T w = r`, and an independent verifier checks that identity using
`Fraction` arithmetic. When the query is outside the row space it returns a
nonzero exact nullspace vector `v` with `Xv = 0` and `rv != 0`. Given compatible
data, `theta0 +/- scale*v` are two compatible models that disagree by the
requested query margin.

The bounded implementation supports at most 32 columns and uses exact
Gauss-Jordan elimination. It does not use a numerical tolerance as proof and
does not infer parameter constraints. Constraints require a separately
verified feasible-set certificate.

For noisy data, the narrow supported certificate takes a supplied `theta0`,
declared `eta`, and exact nullspace witness. It checks
`||X theta0-y||_infinity <= eta`; the two perturbed models have the same
residual and disagree on the query. An identified row-space query receives
only the error bound `eta * ||w||_1`, not an exact point claim. In the TCLab
ARX parameterization, the free intercept and input lag coefficients permit a
constant-heater-column nullspace witness unless additional constraints are
verified.
