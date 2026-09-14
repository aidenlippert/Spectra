# Compact charge metric and an exact obstruction to the simplest implicit bound

Follow-up: [locality and coherent occupation proofs](marginal_coherent_locality.md) reduce the metric to16 factors at pair range two, retain both energy proofs, and replace complement-matrix construction with a checked coherent tree. Only16 of380 ionic configurations are covered in shared branches, and all380 still appear as conditioned transition sources; implicit scaling remains unresolved.

The H6 complement proof now has a **24-factor occupation rule**, replacing 380 independent metric entries in its serialized recipe. It certifies the nearby localized Hamiltonian's ionic complement at **−6.24 Ha**, with an achieved exact weighted-DD floor of approximately **−6.23559189866444 Ha**. The full original-H energy proof uses 16 fresh response directions and retains its **7.416980322924775e−12 Ha** interval width.

This is a reduction in metric description and a tested route toward an operator representation. Complete complement coverage still rebuilds all 380 ionic rows. It is not an implicit proof of those rows or a general representability theorem.

The same 24-factor recipe also transfers under the localized-orbital hopping perturbation t=1/50. The independently recomputed norm shift is 0.04 Ha, so the reference threshold becomes −6.28 Ha. A new 24-step physical upper refinement and 18 fresh response directions recover the perturbed interval width **3.9293895949815965e−12 Ha**, with no dense complement factor. This is the same physical perturbation as the preceding valence-reference transfer, not an additional molecular system.

## The metric is a product of charge factors

At spatial site i define q_i=n_(i,alpha)+n_(i,beta)−1. Its eigenvalues are −1,0,+1 for a hole, single occupation, and a doublon. Use reflection-orbit sums of the features

\[
q_i,\quad q_i^2,\quad q_iq_j,\quad q_i^2q_j^2.
\]

For six sites these give 24 features: three linear, three square, nine charge-pair, and nine square-pair orbit sums. Reflection ties metric parameters only. The verifier makes no reflection assumption about H.

For positive rational factors r_k, the unrounded metric is

\[
u(s)=\prod_{k=1}^{24}r_k^{f_k(q(s))}.
\]

The accepted recipe uses exact nearest-integer rounding, w(s)=round(10^8 u(s)). Its verifier reconstructs these integers from the factors, retains the existing positive-integer bound of 10^12, and checks the actual rounded metric against every rebuilt physical Q row. The recipe stores the state labels but no expanded list of weights. Negative feature exponents are allowed because every factor is strictly positive.

Numerical construction starts with every factor equal to one. It optimizes the weighted row inequalities using explicit Q rows and an analytic gradient, without inheriting the earlier comparison eigenvector or 73-class metric. Rounding is accepted only after exact replay.

The present proposer still forms the dense 380-by-380 complement matrix, and exact replay rebuilds the complete blocks. Eliminating the dense complement factorization and the stored metric vector does not eliminate this matrix construction.

| Recipe | Factors | Exact achieved DD floor, Ha | Exported threshold, Ha |
|---|---:|---:|---:|
| Linear, square, and charge-pair features | 15 | −6.275430349384611 | −6.28 |
| Add square-pair features | 24 | −6.23559189866444 | −6.24 |

These are achieved bounds, not exact optimum claims. The 15-factor version suffices for the earlier unperturbed excitation threshold, but its full-norm shifted bound would miss the tested perturbed excitation threshold. Adding the square-pair features supplies the needed margin.

The factor count for this specified family grows quadratically with the site count. That observation applies to the metric description only: factor discovery, row verification, reference dimension, response dimension, and rational precision have not been shown to scale efficiently.

## Rounding is not the source of the successful metric

There is also a direct bound for the exact unrounded product. Set a=1/(2 w_min). Exact rounding gives (1−a)w_i ≤10^8 u_i≤(1+a)w_i. Thus each unrounded weighted row radius is at most (1+a)/(1−a) times the rounded radius. If all rounded rows have lower bound gamma, the unrounded product has lower bound

\[
\gamma-\frac{2a}{1-a}R,
\]

where R is any upper bound on the rounded weighted row radius. A conservative bound is the off-diagonal CAR coefficient l1 sum times w_max/w_min. Here w_min=2,881,040 and w_max=94,733,721, giving a rounding transfer penalty of approximately **5.186203256566569e−5 Ha**. Consequently the exact product also supports **−6.240051862032566 Ha**. This consequence still relies on the finite complete row proof and its verified weight range; it is not separate implicit coverage.

## The first fully implicit occupation route has a certified ceiling

`marginal_charge_polynomial.py` constructs and replays an occupation-polynomial lower bound without calling a determinant action. For every off-diagonal normal-ordered monomial, its source-occupation indicator and absolute coefficient give a row penalty. Hermiticity pairs the source and target penalties. The elementary two-state positive-square inequality therefore gives an operator inequality H≥F(n), where F is a diagonal occupation polynomial of degree at most four.

On the half-filled equal-spin sector, write D=Σ_i n_(i,alpha)n_(i,beta). On Q, D≥1. The exact verifier checks a proposed decomposition

\[
F-b=\sum_a c_a I_a+(D-1)\sum_b d_b I_b
+(N_\alpha-m/2)X_\alpha+(N_\beta-m/2)X_\beta+r,
\]

where c_a,d_b≥0, each I is a local occupied/empty indicator, and X are unrestricted occupation-polynomial multipliers. The certified lower bound is b−||r||_(coefficient,1). The implementation reconstructs the envelope and every polynomial coefficient from H; it never trusts an exported scalar bound.

The degree-four proposal has 238 positive indicators, 44 charge-localizing indicators, and 509 number-multiplier terms. Independent standard-library replay certifies

\[
QHQ\succeq -6.8803217648028\,Q.
\]

The explicit determinant label **3174** has six particles, three alpha electrons, and one doublon. Its exact value of F is −51602413232767/7500000000000, approximately **−6.880321764368933 Ha**. This one admissible occupation is a ceiling for every sound lower certificate of F. Together with the lower certificate it brackets min_Q F to **4.338666666666667e−10 Ha**.

Therefore refining the positivity decomposition of this particular envelope cannot reach −6.26489910104 Ha. The loss is already in taking absolute values term by term, before the occupation optimization. This ceiling does not constrain the physical Hamiltonian, grouped operator inequalities, the successful weighted-DD cone, or a different reference.

## Implementation and boundaries

`marginal_charge_product.py` proposes and independently replays the factor recipe. The existing `spin_gap` dispatch expands this family and rebuilds all physical blocks using the caller's actual H and P, both directly and inside a transferred reference. Conflicting proof families and incomplete or unphysical coverage retain their existing refusals. No dense complement factor is introduced.

The unperturbed upper witness is inherited from the already certified valence construction; its original-H energy and variance are rechecked. The response basis is rediscovered from zero. The perturbed upper is refined again using its actual Hamiltonian. Thus this experiment isolates the new metric description without claiming a new upper-state discovery.

The polynomial verifier supports up to 64 modes under its stated half-filling and spin-conservation assumptions; the numerical LP proposer is bounded to 16 modes. A 64-mode diagonal Hubbard certificate is checked exactly with determinant actions forbidden by a test. This demonstrates that the verifier itself can avoid sector enumeration on that simple family; it says nothing about achieving useful interacting molecular bounds at that size.

The full marginal suite passes **311 tests in 319.354 seconds**. Six fresh `python -S` replays match every corresponding saved replay field: the degree-four occupation bound, the 15- and 24-factor gap recipes, the 24-factor energy interval, the perturbed energy interval, and its physical upper witness. Focused tests cover operator-envelope positivity with interfering terms, residual protection against a forged bound, local positivity refusals, exact charge-factor rounding, direct and transferred H/P binding, and compact-recipe construction.

The next mathematical target is specific: retain the charge-product weighting and coherent transition amplitudes while certifying its row inequality from occupation rules. Merely increasing the degree of the failed termwise envelope is now ruled out as a solution for this H6 threshold.
