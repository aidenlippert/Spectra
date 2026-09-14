# V5: two costed transfers from acquired modes

## Claim

In the finite instrument family below, one binary observation identifies an unknown dense mode with bounded error. Retaining that direction compiles a measurement for a new coupled system. A second acquired mode enables a further transfer. Each transfer saves expected scalar reads after paying its acquisition cost, compared with any investigator starting each target without information about its immediate predecessor and achieving target error at most 1/20.

This is a theorem about **two successive costed transfers of acquired physical information**. Section 8 additionally proves strictly complementary information value on one common downstream task, including acquisition errors. The model grammar and nulling algorithm are supplied. Structured record replay and a stateful conventional Bayesian method can match the learner. This is not a theorem of algorithmic superiority, unbounded superadditive scaling, unknown-ontology discovery, many-body universality, or unrestricted scientific intelligence. The mathematical probability laws are ideal laws within a declared tolerance; numerical simulations do not certify those laws or their physical applicability.

## 1. Model and operations

Let Q₀,…,Q₃ be the dense orthonormal columns of Hadamard₄/2, c=3/5 and s=4/5. For hidden signs aᵢ∈{−1,+1}, set

\[
r_0=Q_0,\quad u_i=cr_{i-1}+a_i sQ_i,\quad
r_i=-a_i sr_{i-1}+cQ_i.
\]

Induction proves u₁,…,uⱼ,rⱼ,Qⱼ₊₁,…,Q₃ orthonormal. The supplied grammar is common to every method; the signs and resulting laboratory directions are unknown. With public α∈[1,2], define

\[
C_j=I+\sum_{i=1}^j\lambda_i u_i u_i^T,\qquad
(\lambda_1,\lambda_2,\lambda_3)=\alpha(10^{14},10^{10},10^6).
\]

A stage-j component offers only prefix sources 1,…,j. One source read prepares a fresh state and returns 1{|Y+ω|>20}, where Y is a scalar projection along an applied unit vector within norm η=10⁻⁹ of the commanded unit vector; |ω|≤10⁻³. Arbitrary removal of earlier modes, simultaneous/vector observations and uncharged source inspection are excluded from the instrument class. Actuation and readout cannot depend on hidden signs except through acquired records. Fresh baths are conditionally independent.

The admitted harmonic OU source is dq=−Cⱼ⁻¹q dt+√2 dW, with ||q₀||≤1 and stationary law N(0,Cⱼ). Each prepared law must be within total variation εprep=10⁻⁹ of stationarity. [Primary physical context](linear_source_context.md) gives the linear-SDE provenance. Contracts carry component lineage, α and uncertainty bounds. Their physical truth is assumed, not automatically discovered or certified.

## 2. Finite preparation and resource scope

Let V=1+2·10¹⁴. Every read pays relaxation time T=30V=6000000000000030. Its exact finite-time mean and covariance are m=exp(−TC⁻¹)q₀ and Σ=C[I−exp(−2TC⁻¹)]. Whiten by C⁻¹/². The covariance becomes I−R, with eigenvalues r≤exp(−60), and squared mean norm at most exp(−60). Since −r−log(1−r)≤r²/[2(1−r)],

\[
\mathrm{KL}(N(m,\Sigma)\Vert N(0,C))
\le\tfrac12e^{-60}+\frac{e^{-120}}{1-e^{-60}}<2^{-60}.
\]

The last inequality uses e>2. Pinsker gives TV<2⁻³⁰<10⁻⁹. This proves the ideal SDE preparation bound. Coefficient errors, bath deviations and numerical sampling are not certified by that calculation; they need a separate allocation in an actual instrument's total error budget.

Per read, the model charges one reset to the bounded initialization class, four coordinate preparations, T time units, one unit-vector command (at most 64 bits per rational coordinate), and one binary record. Reset and source construction are available primitives. Their hardware costs are **not derived here**. If their common finite costs are Rreset, Rsource and Rread, every read comparison is multiplied by Rreset+Rsource+T+Rread, with digital processing added separately. This is primitive accounting within an instrument model, not an apparatus or numerical energy-cost theorem.

Gaussian amplitude is unbounded. Expected quadratic energy at readout is at most 5/2 kBT, and E||q||²≤5+2(10¹⁴+10¹⁰+10⁶). For n reads, Markov and a union bound give probability at most n[5+2(10¹⁴+10¹⁰+10⁶)]/R² of any **readout endpoint** exceeding radius R. This is not a pathwise bound during relaxation. A bounded apparatus needs an additional tail/model budget. The saturating binary sensor does not remove this assumption. The enormous condition number and time and tight actuation tolerance preclude a claim of practical efficiency.

The digital learner has no hidden search oracle. At depth j≤3, all reconstructed mode/probe coordinates have denominator dividing 2·5ʲ and numerator magnitude at most that denominator, so their stored rational integers need at most eight bits. A call validates and replays at most two earlier records, constructs three four-coordinate vectors, and stores the current bit and contract. If ℓ counts the public contract's complete description, including identifiers and α, the fixed-depth learner uses a constant number of rational operations and conversions on at most O(ℓ+1)-bit contract integers, plus bounded-size mode arithmetic. Elementary schoolbook arithmetic gives the conservative bound O(K(ℓmax+1)³) over K acquisitions, with O(ℓmax+1) live space per prefix. This describes accepted finite inputs; host input-size limits may reject larger encodings. This is a bound for this finite grammar, not a scaling theorem for arbitrary matter. Source execution is charged separately. [Computational accounting](../../results/v5/computational_accounting.json) profiles fixture generation, both studies, every baseline and serialization, with excluded costs stated.

## 3. Warm identification and robustness

Conditional on a correct retained prefix, command xⱼ=−srⱼ₋₁+cQⱼ. For aⱼ=+1 this is rⱼ, orthogonal to all active modes. For aⱼ=−1 it has projection −24/25 on uⱼ and zero on all earlier modes. The measured fluctuation therefore identifies the new sign without receiving a simulator label.

Actuation error implies null variance ≤1+η²Σλᵢ<1001/1000; the alternative projection magnitude is at least 24/25−η>19/20. A null false positive requires |Y|>19.999, and Chebyshev bounds its probability by (1001/1000)/(19.999)²<1/399. An alternative miss requires |Y|≤20.001. Its Gaussian density is bounded by 1/(√(2π)σ), giving miss probability <20.001/[(19/20)√λⱼ]. Thus valid stage error bounds are

\[
e_1=e_2=1/399,\qquad e_3=11/500.
\]

TV errors add εprep per read. Conditional union accounting includes all earlier mistakes:

\[
P(\text{any wrong sign on a three-mode lineage})
\le5389/199500+3\cdot10^{-9}<1/20.
\]

This is a per-lineage guarantee, not simultaneous coverage of every saved target. The controller checks stored evidence, lineage, precision and confidence before executing. It cannot infer every undeclared model mismatch from those checks.

## 4. Cold one-read lower bound

Grant all older modes exactly. In an orthonormal coordinate system write the unknown predecessor and residual as uₐ=ce₁+as e₂, tₐ=−as e₁+ce₂, and target mode vₐᵦ=ctₐ+bs e₃, with independent fair a,b. The covariance is B+Λuₐuₐᵀ+λvₐᵦvₐᵦᵀ, B≥I independent of a,b, Λ/λ≥10⁴.

For every scalar probe x,

\[
\sum_a(x\cdot u_a)^2\ge(9/16)\sum_a(x\cdot t_a)^2.
\]

At least one a satisfies the corresponding pointwise inequality. For that a, the two target-sign variances have midpoint D and half-difference E with

\[
D\ge\Lambda(9/16)(x\cdot t_a)^2+\lambda s^2(x\cdot e_3)^2,
\quad |E|=2\lambda cs|(x\cdot t_a)(x\cdot e_3)|.
\]

The arithmetic-geometric mean inequality yields ρ=|E|/D≤(4/5)√(λ/Λ)≤1/125. For centered Gaussians of variance ratio 1+t with t=2ρ/(1−ρ), KL≤t²/4 and Pinsker gives TV≤ρ/(1−ρ). The other predecessor sign has TV≤1. Convexity gives mixture TV≤[1+ρ/(1−ρ)]/2, hence target Bayes risk

\[
r\ge\tfrac14-\frac{\rho}{4(1-\rho)}\ge123/496.
\]

Finite preparation reduces this bound by at most εprep, so r₀=123/496−εprep. The lower bound allows arbitrary scalar probes, full analog outputs and arbitrary postprocessing, which are stronger than the implemented binary sensor. Earlier-prefix reads have no fresh-target information. Randomized first probes are permitted but must not use hidden-state-correlated prior records. The bound does not apply to an investigator already possessing the predecessor information: that is the resource being acquired.

## 5. Adaptive expected read bound

For any variable-horizon investigator with target error R≤δ=1/20, let x=P(N≥1), y=P(N≥2). The decision to start and first query must be made without hidden-state-correlated records. Zero-read error is 1/2, so x≥1−2δ. Immediately after the first read, risk is at least r₀. Additional reads can reduce posterior error by at most 1/2 on the event of continuation. Therefore

\[
R\ge\tfrac12(1-x)+r_0x-\tfrac12y,\qquad
E[N]\ge x+y\ge(1-2\delta)(1+2r_0)=L.
\]

Numerically L≈1.3463709659>1.25. This permits randomized selective stopping and unlimited later reads. E[N]≥x+y is an inequality when three or more reads are possible.

## 6. Two acquisition-inclusive improvements

Acquire u₁ once and solve four new stage-2 components with independent signs and changed α in one read each: 1+4=5. Four targets starting without predecessor information require expected cost ≥4L>5.38548. Net advantage exceeds 0.38548 reads after acquisition.

Retain each newly observed u₂. Four fresh stage-3 components now take one read each; a frozen-after-u₁ investigator lacks u₂. Charging the entire one-read u₂ acquisition again in this second comparison still yields 1+4<4L. This proves a second positive costed increment. These are separate counterfactual comparisons; their savings must not be summed as independent costs against a single baseline.

The executed four-by-four tree contains 20 targets and costs 21 cumulative reads. Implemented frozen-after-first costs 37; implemented scratch and exact-key retrieval cost 56. These are algorithm costs, **not minimum-cost theorems** of two/three reads. The universal cold lower bound is L. Stateful Bayes and structured replay can acquire the same information; their observed cost is also 21.

Later component identities and couplings differ; target signs are fresh. Records determine a direction that changes the next physical observation. Unlimited computation on a single uninformed scalar observation cannot replace that acquired direction by the lower bound. Nevertheless the recursive grammar and conversion of signs into modes are supplied. The theorem establishes bounded experimental transfer of acquired information, not invention of that grammar or algorithm.

## 7. Execution and failure conditions

The fixed 64-world protocol contains 256 stage-2 and 1024 stage-3 targets. Identifiers are generated independently of signs. Learner functions receive contracts and binary callbacks; evaluator truth is recorded only after decisions. No teacher labels or parameter fitting occur. Methods have equal source access and paired normal variates. Shared-prefix errors induce task dependence; pooled errors are descriptive counts, not iid confidence intervals.

Verification replays source calls and learned evidence, recomputes ledgers, and checks exact inequalities. Negative controls cover declared precision/lineage violations, corrupt evidence and intentionally undeclared control mismatch. Success outside the contract is not promised or automatically diagnosed. The source family, prefix availability, reset primitive and range assumptions remain substantive limitations. No hardware experiment has been performed.

## 8. A common-baseline complementarity theorem

Two separate savings comparisons alone do not prove complementary knowledge value. We therefore define a second, stricter quantity on the **same final-stage target**. Let R∅, R_A, R_B and R_AB denote the minimum average error of a one-scalar investigator given neither prefix sign, only a, only b, or both. The target c is freshly sampled independently after calibration. Only the declared sign bits are retained; B-only access does not include a probe or full mode vector that encodes a. The investigator knows the family and can optimize any scalar probe and postprocessing. For gain G(S)=R∅−R_S, define

\[
J=G(AB)-G(A)-G(B)=R_A+R_B-R_{AB}-R_\varnothing.
\]

Positive J means the second piece of information is more valuable when the first is available. It is a common-task statement, not the sum of the two resource savings in section 6.

**A-only lower bound.** Section 4 applies to unknown b with a known, giving R_A≥123/496−εprep.

**B-only lower bound.** For fixed b, sum over the two unknown a values. In span(Q₀,Q₁,Q₂), orthonormality gives

\[
\sum_a(u_1u_1^T+u_2u_2^T)=2I-\sum_a r_2r_2^T.
\]

The residual sum has eigenvalues 2s⁴=512/625, 2(1−s⁴)=738/625 and 0. Its largest eigenvalue is below 32/25, so

\[
\sum_a[(x\cdot u_1)^2+(x\cdot u_2)^2]
\ge(9/16)\sum_a(x\cdot r_2)^2.
\]

Since λ₁≥λ₂, at least one a has sufficient nuisance variance relative to the target residual. The same D,E and Gaussian-TV argument uses λ₂/λ₃≥10⁴ and proves R_B≥123/496−εprep. This is a different unknown-sign conditioning from the A-only proof; the residual Gram bound is essential.

**Both upper bound.** With true signs, the warm target probe gives R_AB≤11/500+εprep.

**Neither upper bound.** Use the fixed unit probe null_probe((+1,+1)). Its correct-prefix branch has weight 1/4 and error at most 11/500+εprep. For each of the other three prefixes, write t=x·r₂ and z=x·Q₃. Uniformly for a unit actuator perturbation of size η,

\[
D\ge\sum_{i=1}^2\lambda_i\max(|x\cdot u_i|-\eta,0)^2,
\quad
|E|\le2\lambda_3cs(|t|+\eta)(|z|+\eta).
\]

Exact rational arithmetic verifies D≥125|E| for all three mismatches; the common factor α cancels. [The executable certificate](../../experiments/v5_complementarity.py) records every numerator and denominator. Thus their target TV is at most 1/124. Any fixed binary classifier has equal-prior error at most 1/2+TV/2, and preparation adds εprep. Consequently

\[
R_\varnothing\le\tfrac14\tfrac{11}{500}
+\tfrac34(\tfrac12+\tfrac1{248})+\epsilon_{prep}
=47557/124000+\epsilon_{prep}.
\]

Combining lower and upper bounds gives, for true sign information,

\[
J\ge2243/24800-4\epsilon_{prep}.
\]

**Acquisition errors.** The actual first stored bit differs from a with probability ≤e=1/399+εprep. The second is obtained using the first, so its error is ≤2e, as is the joint-prefix error. Couple actual and true side-information experiments using the same hidden prefix, calibration randomness and future target. Their transcripts differ only when a retained bit is wrong; under an identical policy, all further probe/output laws agree on the matching event. This coupling is uniform over policies, so optimal risks differ by at most e for A, 2e for B and 2e for AB. Taking the unfavorable directions yields

\[
\boxed{J_{\mathrm{acquired}}
\ge\frac{2243}{24800}-\frac5{399}-9\cdot10^{-9}
=\frac{963696138679}{12369000000000}>0.077.}
\]

This is a conservative positive lower bound on complementary value of **noisily acquired** information. Preparation contributions occur in both target risk and acquisition coupling because they concern different reads; none is silently omitted. A and B acquisitions precede independent fresh targets and have separate randomness. Actual full mode records are retained for the main transfer experiment; the controlled-retention study explicitly strips them down to the allowed bit subset.

The simulation evaluates four explicit policies, not the unknown global optimum over all scalar probes. Its observed contrast is therefore a descriptive cross-check; the optimal-risk guarantee comes from the mathematical bounds. Calibration uses two reads per source lineage. Four counterfactual target reads are charged in the study; the result does not turn four counterfactual experiments into one free physical read. Together with section 6, the two results establish a bounded operational sense of compounding: acquired pieces have strictly complementary downstream value and support two successive acquisition-inclusive cost reductions. General discovery of the representation grammar remains unproved.
