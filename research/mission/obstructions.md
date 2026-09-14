# Finite-budget obstructions to universal physical-invention certification

This note makes the negative claim precise. A finite experimental/computational
budget cannot certify a uniform statement about arbitrary physical systems when
the admissible model class has no quantitative regularity, complexity, or noise
assumptions. The result is an information limitation, not a claim that useful
materials science is impossible.

## 1. Adaptive compact-bump indistinguishability

Let the intervention/design domain be \(X=[0,1]^d\), and suppose an algorithm
makes at most \(q<\infty\) adaptive noiseless queries of an unknown response
\(f:X\to\mathbb R\). Exact sampled values can be certified; the impossibility
concerns uniform properties whose completeness requires conclusions about
unsampled points. No such nontrivial guarantee is valid over \(\mathcal F=C(X)\).

For a deterministic algorithm, run it against (f_0\equiv0), obtaining query
points (x_1,ldots,x_q). Since this is a finite set, choose a nonempty open
ball (B\subset X\setminus\{x_1,ldots,x_q\}). A continuous bump

\[
 h(x)=a\,\max\{1-\|x-z\|/r,0\},\qquad \overline B(z,r)\subset B,
\]

has (h(x_i)=0) for every query and (|h|_\infty=a). Consequently the
entire transcript under (f_0) and (f_1=h) is identical: induction on the
query index gives the same previous answers, hence the same next query, and all
answers are zero. The algorithm therefore returns the same certificate for two
functions separated by (a) in sup norm. It cannot certify, for example,

\[
\|f\|_\infty\le \varepsilon
\]

for every continuous (f) when (a>\varepsilon), nor distinguish the
existential claim (exists x:f(x)\ge a/2) from its negation.

The argument survives arbitrary adaptivity: adaptivity selects points from the
observed transcript, and the two transcripts are equal. It also survives any
finite amount of internal computation. If the algorithm is randomized, apply
the same construction conditionally to each random seed to show that no
*pointwise* guarantee can hold for all seeds. For a distributional statement,
use the finite packing below, which gives a quantitative failure probability.

### Quantitative randomized form

Pack \(M\) pairwise disjoint balls \(B_1,\ldots,B_M\) of radius \(r\) in \(X\),
and choose bumps \(h_j\) of height \(a\), supported in \(B_j\). Use a prior
that assigns probability \(1/2\) to \(f_0\) and \(1/(2M)\) to each \(h_j\).
Couple every run to a baseline run under \(f_0\), using the same algorithmic
random seed. Until a baseline query lands in \(B_J\), the transcripts and
queries are identical. At most \(q\) balls can be hit on any baseline path;
hence

\[
\Pr(\text{hit }B_J)\le q/M.
\]

The coupling fails with probability at most \(q/M\), so the total variation
distance between the transcript law under \(f_0\) and the mixture
\(Q=M^{-1}\sum_jP_{h_j}\) is at most \(q/M\). Binary testing with prior
probabilities \(1/2,1/2\) therefore has Bayes error
\[
R^*\ge \frac{1-\mathrm{TV}(P_{f_0},Q)}2\ge\frac{1-q/M}{2}.
\]
This Bayes risk lower-bounds the worst-case (minimax) risk over
\(\{f_0,h_1,\ldots,h_M\}\). Taking \(M>2q\) gives a constant failure
probability. This is a coupling/TV statement; conditioning on “no hit” can
still reveal which balls were ruled out, so it is not literally no information
about \(J\).

The same proof applies to physical design if an unqueried intervention can
change the relevant outcome while agreeing on all queried interventions. Thus a
certificate based only on a finite set of simulations or experiments cannot
establish uniform behavior over an unrestricted intervention space.

## 2. Noise imposes a separate budget lower bound

The bump obstruction is present even with exact observations. Noise creates a
second, independent limit. Consider two hypotheses whose queried scalar outcome
is Gaussian:

\[
H_0:Y_i\sim N(0,\sigma^2),\qquad
H_1:Y_i\sim N(\Delta,\sigma^2),
\]

independently across (i=1,ldots,n). Their KL divergence is

\[
 D_{\rm KL}(P_0^n\|P_1^n)=n\Delta^2/(2\sigma^2).
\]

By Pinsker, (\mathrm{TV}(P_0^n,P_1^n)\le
\sqrt{n\Delta^2/(4\sigma^2)}). Le Cam's two-point inequality therefore gives,
for every test (T),

\[
\max_k P_k(T\ne k)\ge \tfrac12(1-\mathrm{TV})
\ge \tfrac12\left(1-\frac{\sqrt n\,|\Delta|}{2\sigma}\right).
\]

In particular, achieving error at most (\delta<1/2) requires

\[
 n\ \ge\ 4\sigma^2(1-2\delta)^2/\Delta^2.
\]

Constants are inessential; the (\sigma^2/\Delta^2) dependence is the point.
For adaptive experiments, the chain rule for KL gives

\[
 D_{\rm KL}(P_0^{T}\|P_1^{T})
 =\mathbb E_0\!\left[\sum_{i=1}^{T}
 D_{\rm KL}(P_0(Y_i\mid H_i)\|P_1(Y_i\mid H_i))\right],
\]

so the same lower bound applies whenever every allowed measurement separates
the hypotheses by at most (\Delta) at noise level (\sigma). More generally,
the required budget is controlled by the accumulated distinguishability (KL,
Hellinger, or total variation), not by the number of model parameters alone.

For invention, this rules out a uniform claim that two candidate processes or
materials differ by a tolerance smaller than the measurement floor without
paying the corresponding replication cost. It also means that a favorable
optimizer output cannot be certified merely because its predicted margin is
small relative to unquantified experimental noise.

## 3. Infinite model classes and no finite uniform sample complexity

The class (C([0,1]^d)) contains the bump family above at arbitrarily small
scales. For every finite budget (q), choose (M>2q) disjoint supports (or,
for a fixed target region, choose a support avoiding the realized query set).
Thus the minimax error for sup-norm recovery over this class is bounded away
from zero for every finite (q). This is stronger than saying the class is
“large”: it has infinite metric entropy at every sufficiently fine scale, and
there is no uniform query bound independent of a modulus of continuity or a
complexity budget.

The same obstruction can be encoded in an infinite physical model class: models
agree on every queried intervention and differ only on one of infinitely many
unqueried, physically admissible regimes. A finite transcript identifies only an
equivalence class of models. A universal certificate is possible only if the
property is constant on that equivalence class, which is precisely an
assumption about the model class, the intervention coverage, or both.

## 4. What assumptions make a constructive version possible?

There is no single magic assumption. A useful finite-budget theorem needs the
following ingredients, each tied to a specific failure above.

1. **A declared scope and target tolerance.** Specify a compact intervention
   domain (X), operating conditions, admissible preparation policies, loss or
   constraint functionals, tolerance (\varepsilon), and failure probability
   (\delta). “All matter” without a domain and tolerance is not a statistical
   estimand.

2. **Quantitative restriction on the response/model class.** Examples include
   a known Lipschitz/Hölder modulus, bounded variation, a finite-dimensional
   parametric family with identifiable parameters, or a finite metric entropy
   bound (\log N(\mathcal F,\|\cdot\|,\alpha)\). For a Lipschitz class on a
   (d)-dimensional domain, a covering-grid argument gives finite uniform
   certification: sample on a mesh of spacing at most
   (\varepsilon/(2L)), estimate each value to (\varepsilon/2), and use the
   modulus to extend the bound between mesh points. The number of sites scales
   as (O((L/\varepsilon)^d)), before noise replication. The curse of
   dimensionality is explicit rather than hidden.

3. **A quantified observation/noise model.** State independence or dependence,
   tail bounds or a likelihood, and a calibration bound for each measurement.
   Confidence intervals must be simultaneous over the adaptively selected
   interventions (for example via a valid confidence sequence, a union bound
   over a finite cover, or an explicitly justified adaptive-design theorem).

4. **Coverage and identifiability.** Every property to be certified must be
   distinguishable by allowed interventions and measurements. If two admissible
   physical models induce the same law for all permitted observations but differ
   on the claim, no algorithm can certify that claim; this is an exact
   necessary condition, not merely a practical warning.

5. **A realizable execution and validity bridge.** The certificate must attach
   to the physical preparation and operating policy, including process history,
   not only to an unchecked simulator. One needs either a validated forward
   model with a bounded discrepancy, or a physical test whose measurement model
   directly bounds the required behavior. If
   ( |J(d)-\widehat J(d)|\le e(d)\) uniformly on the searched set and the
   optimizer is (\eta)-optimal for \(\widehat J), then the selected design is
   (2\sup e+\eta)-optimal; without this bridge the optimization guarantee is
   purely internal to the model.

Under these assumptions, a constructive workflow is defensible: maintain a
confidence region over the restricted model class, choose interventions that
either shrink that region or produce a certified feasible policy, and stop only
when the simultaneous error budget is below the requested tolerance. The
resulting theorem is scoped to the declared domain, class, noise, and budget.
It does not certify arbitrary matter or extrapolation outside them.

## Primary references

- L. Le Cam, “Convergence of estimates under dimensionality restrictions,” *Annals of Statistics* 1 (1973), 38–53. DOI: [10.1214/aos/1176342360](https://doi.org/10.1214/aos/1176342360) (two-point/minimax method).
- A. B. Tsybakov, *Introduction to Nonparametric Estimation*, Springer (2009), especially Chs. 2–3. DOI: [10.1007/978-0-387-79052-7](https://doi.org/10.1007/978-0-387-79052-7) (packing, testing, and entropy lower bounds).
- T. M. Cover and J. A. Thomas, *Elements of Information Theory*, 2nd ed., Wiley (2006), Ch. 11. DOI: [10.1002/047174882X](https://doi.org/10.1002/047174882X) (KL chain rule and Pinsker inequality).
- D. Haussler, “Decision theoretic generalizations of the PAC model for neural net and other learning applications,” *Information and Computation* 100 (1992), 78–150. DOI: [10.1016/0890-5401(92)90010-D](https://doi.org/10.1016/0890-5401(92)90010-D) (finite complexity/covering requirements for uniform learning).
- J. L. Doob, “Application of the theory of martingales,” in *Le calcul des probabilités et ses applications* (1949), 23–27. (Foundational confidence-sequence/martingale perspective; modern statements should specify the chosen concentration theorem.)

The bump proof above is self-contained. The cited works support the standard
testing, information, and complexity principles; they do not imply that the
physical invention objective itself is impossible.
