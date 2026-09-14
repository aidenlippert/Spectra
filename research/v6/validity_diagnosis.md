# Post-run diagnosis: constant input does not identify new input effects

This note was written after the frozen V6 test. It does not revise its selection rules, coverage target, data split, or reported forecasts.

The measured TCLab step data have Q1=50 and Q2=0 throughout the fitting record. Selection and calibration use later portions of the same constant-input record. The sine experiment varies Q1 between 0 and 100. See the [primary experiment description](https://dowlinglab.github.io/pyomo-doe/notebooks/tclab-model/) and [frozen source receipt](primary_sources.md).

Consider even the affine one-step family

\[
y_{t+1}=a^T H_t+b^Tu_t+d.
\]

If all training inputs equal u₀, then for every vector v the coefficients

\[
b'=b+v,\qquad d'=d-v^Tu_0
\]

give exactly the same training predictions. Under a new input u they differ by vᵀ(u−u₀). Consequently these observations alone cannot identify the response to a new input contrast. More history or an invertible change of coordinates cannot restore information absent from the experiment. Additional physical restrictions can constrain this ambiguity; such restrictions would need to be stated and validated.

In V6, centering the measured training inputs gives rank zero. All 55 sine-test forecast blocks contain a contrast outside that training span. The diagnostic in `experiments/v6_audit.py` uses training and planned **input values only**; it does not inspect the test outputs to decide which blocks are unsupported. A unit test verifies this with unavailable/nonfinite target outputs.

This is a necessary excitation test, not a sufficient identification theorem. Rank-complete inputs alone do not prove hidden-state observability, independent residuals, correct model order, or calibrated uncertainty. The synthetic training inputs have rank two, yet their empirical forecast-block coverage still falls below the requested average threshold under withheld forcing. Serial dependence, distribution change and fitted-model uncertainty remain relevant.

The appropriate conclusion is a validity failure: low error on the constant-input validation segment did not justify narrow intervals on a different forcing regime. The observed measured coverage was 2/55≈3.64%, not the requested 90%. V6 does not retrospectively count an input-span warning as successful prospective abstention; it records the diagnosis separately. A new refusal policy would need a new frozen evaluation.

This failure is distinct from the headroom null. The privileged companion representation and stateful ARX are algebraically equivalent regardless of whether their common prediction law is well calibrated. The test therefore exposes both an unsuitable candidate for a method-learning investment and a limitation of the physical data for intervention claims.
