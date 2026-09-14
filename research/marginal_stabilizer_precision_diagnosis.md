# Reduced quartic precision diagnosis

The tightened replay in `results/marginal_stabilizer_refine` does not improve the
accepted bound: Clarabel tolerance `1e-11` and factor denominator `10^12` give
lower `3.2758655432377815`, versus `3.275865704229399` previously. The change is
only `1.6e-7`, so neither solver tolerance nor Gram-factor decimal rounding
explains the approximately `0.00514` gap.

The multiplier serialization error is bounded independently. For the full
degree-three multiplier basis, the coefficient map has 8,351 columns; their
total absolute column l1 mass is 60,670 and the largest is 15. Rounding each
multiplier coordinate to `10^-10` therefore contributes at most

`(60,670)(0.5)(10^-10) = 3.0335e-6`

to the canonical residual l1 norm. The actual symmetry-filtered 541-column
basis is a subset, so this is a conservative bound. It is over three orders
of magnitude too small to explain `0.00514`.

The loss is therefore attributable to the numerical reduced certificate's
residual / PSD reconstruction and its interaction with the original export,
not multiplier rounding. A separate full Reynolds reduction subsequently
reached numerical objective `3.281006695866` and exact lower `3.281005630613`
(width about `1.065e-6`), showing that quartic certificates can reach that narrower interval. This
is consistent with a conditioning/reconstruction problem in the earlier run;
these measurements do not isolate the exact numerical contribution responsible.
