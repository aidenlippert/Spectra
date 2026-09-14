# Exact obstruction to independently bounded transition ratios

The proposed shared D≥2 tail bound fails even when every individual metric-ratio maximum is computed exactly. Its conservative row polynomial has minimum in

\[
[-10.204432073439493,\;-10.204432073014694]\ \mathrm{Ha},
\]

an exact rational bracket of width **4.2479938036440263e−10 Ha**. The required threshold is−6.264 Ha. Thus stronger positivity solvers cannot rescue this particular envelope. The existing complete charge-spin proof remains valid; this is an obstruction to the proposed compression, not to that proof or general quantum chemistry.

## The attempted shared tail inequality

Keep the same positive charge-product metric as the accepted H6 proof. For each grouped transition g, certify a fixed sign sigma_g on its allowed source occupations and bound its metric ratio by rbar_g over the fixed-spin domain with at least two doublons. Define

\[
F(n)=V(n)-\sum_g \bar r_g I_g(n)\sigma_g A_g(n).
\]

On that domain, sigma_g*A_g=|A_g| wherever the event I_g is active. Dropping the target-valence exclusion only adds nonnegative penalties. Hence the actual weighted Q row R(n) obeys R(n)≥F(n). All transitions, including transitions to the one-doublon sector, remain included. Combining a successful tail inequality with the existing one-doublon proof would therefore use one consistent global metric; there is no invalid stitching of unrelated block bounds.

For the admitted one-/two-body Hamiltonian, the compiler reconstructs a degree≤4 polynomial and checks that degree limit explicitly. A single-excitation event times an affine spectator amplitude has degree≤3; a pure double-excitation event has degree≤4.

## Exact ratio maximization with finite memory

The initial product of independent local-factor maxima gives an LP proposal around−30.5283 Ha. It can combine local maxima from incompatible charges. The replacement DP eliminates that particular loss.

For each transition, process sites in order. Its state stores the assigned alpha count, beta count, capped doublon count, and the previous R site charges. Each positive metric factor is multiplied exactly once when its last site is assigned. Paths with the same state are merged by taking the larger rational product. Positive multiplication preserves that ordering. Population feasibility is enforced throughout and exactly at the end.

For the actual range2 metric, each ratio DP has memory2. Across1140 transitions the calculation visits **103020 DP states**, makes **134188 state transitions**, and computes **220192 scalar metric powers**. These are operation counts across separate transition problems, not counts of distinct physical configurations eliminated. No full occupation endpoints or determinant actions are used by the DP. The implementation refuses interaction ranges beyond4.

The largest ratio bound falls from about227.389 for independent factor bounds to the exact source-domain maximum about126.668. However, different transitions still maximize their ratios independently. Their maximizing charge arrangements need not agree.

## A checked numerical limit, not merely an unsuccessful LP

The new envelope has553 nonzero monomials. The degree4 positivity LP has794 coefficient equations and10857 columns. Its numerical proposal is approximately−10.204432073014774 Ha.

Exact replay checks

\[
F-b=P_4+(D-2)P_2+(N_\alpha-3)L_\alpha
                    +(N_\beta-3)L_\beta+\varepsilon,
\]

where P4 is a nonnegative sum of degree≤4 occupation indicators, P2 is a nonnegative sum of degree≤2 indicators, and the number multipliers have degree≤3. The absolute coefficient sum of epsilon gives a rigorous residual allowance. The exported witness uses240 positive indicators,43 charge localizers, and509 number-multiplier terms.

A single admissible occupation supplies the matching polynomial upper bound. Exhaustive discovery scanned200 D≥2 occupations, but final upper-bound replay checks only the selected witness. This is an upper bound on **min F**, not an upper bound on a Hamiltonian eigenvalue or the true weighted row minimum.

The bracket proves that the poor threshold is inherent in F. Increasing the degree of a positivity certificate for unchanged F cannot lift its lower bound above its actual minimum.

## The exact loss at the witness

The witness has occupation integer819 and charges

\[
q=(1,-1,1,-1,1,-1).
\]

There are no singly occupied sites, so no residual spin ambiguity is involved. Direct compilation of the same metric-weighted Q row gives

| Quantity at the same occupation | Value (Ha) |
|---|---:|
| True grouped weighted Q row | −6.185681551371275 |
| Independently maximized-ratio envelope | −10.204432073014694 |
| Lost lower-bound strength | 4.018750521643418 |

The true row passes−6.264 Ha. The envelope fails because it combines conservative penalties that are not simultaneously realized at this charge pattern. Target-valence penalties are also retained conservatively; the diagnostic does not attribute every part of the loss solely to ratio maximization.

## Consequence for the next attack

The finite-memory DP solves each individual ratio maximization exactly. The remaining issue is joint dependence across transitions. A next candidate must retain that dependence in a shared charge expression, or choose a metric with simpler tail behavior and reverify every source row. Repeating the LP at higher precision or adding positivity constraints to the unchanged envelope cannot solve this obstruction.

No new energy certificate is accepted in this experiment. The preceding charge-spin certificates and direct perturbed energy interval remain the accepted finite results. Charge enumeration, full-sector witness/response construction, larger-system transfer and general marginal representability remain open.


## Two simpler metric families are excluded exactly

A second attack tried simplifying the metric itself. Numerical feasibility searches used the explicit380-state Q matrix and checked the normalized physical row bounds, rather than trusting an LP status alone. The best located bounds were approximately−6.52227 for a doublon-count-only metric and−6.45715 for a constant-tail metric with independent one-doublon weights. Those numerical values are diagnostics, not exact family optima.

Exact physical-row certificates now establish the more useful statement: **neither family can reach−6.264 Ha through weighted Q-row dominance**, even with unrestricted strictly positive class weights.

| Metric family | Free positive class weights | Exact witness rows |
|---|---:|---:|
| Depend only on doublon count1,2,3 | 3 | 3 |
| Independent weights for each ordered doublon/holon pair at D=1; one common weight for all D≥2 | 31 | 31 |

For a row s and class map c, define

\[
A_{sj}=(H_{ss}-\gamma)\mathbf1_{c(s)=j}
       -\sum_{t\in Q,t\ne s}|H_{ts}|\mathbf1_{c(t)=j}.
\]

A successful metric must satisfy Aw≥0 with every w_j>0. Each certificate supplies positive rational source-row weights y_s for which **every component of yᵀA is strictly negative**. Then yᵀAw<0 for every positive w, contradicting the required inequalities. No numerical normalization, lower weight cutoff, or upper weight cutoff is needed by the exact proof.

The verifier rebuilds the selected rows from the actual Hamiltonian, validates their fixed-spin ionic occupations, and removes target-valence transitions exactly. Three or31 determinant actions suffice for replay. Discovery used all380 ionic states and remains explicitly finite. The least negative coefficients are approximately−0.0898886 and−0.00592866, respectively, so rounding is not near a sign boundary.

These are restrictions on two families of weighted-row certificates. They do not exclude a nonconstant tail metric, a different operator positivity decomposition, or a scalable representation. The accepted range2 charge-product metric already lies outside both excluded families. Within this proof route, the tail must retain some charge-arrangement dependence.


## What is tractable, and which dependence to restore

At fixed interaction range R, each ratio DP has at most (p+1)^2(d0+1)3^R states per site, where p is the target spin population and d0 the capped doublon threshold. It explores at most four local choices per state. Including m layers gives a polynomial state-count bound in m for fixed R; exact rational bit costs and the number of separate transitions must also be counted. This is a bound for the individual ratio subproblem, not for minimizing the joint row polynomial.

The fixed31-row obstruction also supplies a local metric-feature diagnostic. At the uniform metric, the largest normalized derivatives among the accepted range2 charge features favor increasing q1*q3+q2*q4, decreasing q2*q3, and decreasing q1²+q4² (zero-based spatial sites). Their normalized derivative magnitudes are approximately0.30592,0.30413,0.23595. An independent aggregation by target occupation reproduces every exact derivative. These are spatial charge features, not alpha/beta spin correlations. The calculation reuses only31 physical source actions and is preserved in `dual_feature_directions.json`. Improving one separating functional locally does not establish full metric feasibility.


Verification: **339 marginal regression tests passed in262.594 seconds**. Three new standard-library replays match every saved receipt field: the exact tail-envelope bracket and both physical metric-family obstructions. The prior successful energy certificates are unchanged. See `results/marginal_h6/charge_tail/progress.json`.
