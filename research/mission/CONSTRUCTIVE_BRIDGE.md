# From physical evidence to reachable, certified design

This is a constructive reduction with explicit limitations, assembled from
standard concentration, covering, control and optimization arguments. It is not
a novelty claim, a universal solver, or physical validation. Its purpose is to
expose the shared mathematical obligations in prediction, preparation and
invention, and make incorrect compositions falsifiable.

## 1. The object being optimized is an execution policy

Let D be a compact, effectively coverable space. A decoder maps d in D to a
finite, physically admissible preparation-and-operation policy pi_d, including
feedstocks, transformations, intermediate measurements where represented, and
operation through a stated horizon. Establishing this decoder's physical
executability is an obligation, not a free synthesis oracle. A finite union of
bounded control-word spaces is one admissible construction. It covers only those
words, not arbitrary chemistry. Feedback policies require an explicit causal
parameterization and sufficient history state.

Let J(d) be an actual physical expected resource/loss and g_j(d)<=0 the physical
requirements, j=1,...,m. A pathwise requirement may be represented by
g_j(d)=Pr_pi_d(path violates requirement j)-delta_j. An expected-temperature
constraint is not a pathwise temperature guarantee. Budget limits and operating
conditions are included explicitly. Write F={d:g_j(d)<=0 for all j}, and
V*=inf_{d in F} J(d). Initially F may be empty.

There are two different probability levels. Confidence 1-alpha concerns whether
inferred bounds hold. Execution reliability 1-delta concerns a fresh run of the
returned policy. Conditional on valid bounds, certification of a failure
probability at most delta gives marginal failure at most alpha+delta (more
precisely alpha+(1-alpha)delta). Neither parameter absorbs unknown model error.

## 2. The deterministic certificate, valid after adaptive selection

Maintain a finite covering by cells C of D. For every cell and quantity f_j
(f_0=J), let [L_j(C),U_j(C)] enclose *every* f_j(d), d in C, simultaneously.
Also maintain point bounds for actually executable candidate policies. Define

\[
 A=\{C:\ L_j(C)\le0\ \text{for every constraint }j\},\qquad
 \underline V=\min_{C\in A} L_0(C).
\]

A candidate a is certified feasible if U_j(a)<=0 for every constraint. Among
these candidates choose one with smallest upper objective, denoted U_best.

**Certificate theorem.** On the simultaneous enclosure event:

* if A is empty, F is empty within the decoded, covered policy class;
* if a exists, it is feasible, V* lies in [underline V,U_best], and
  J(a)-V* <= U_best-underline V;
* otherwise the correct answer is unresolved, with the surviving cells and
  uncertainty widths retained.

**Proof.** Every feasible d belongs to a cover cell C. Since L_j(C)<=g_j(d)<=0,
that cell is in A. Thus underline V<=J(d) for every feasible d; take the infimum.
Point upper bounds make a feasible and give J(a)<=U_best. Both claims follow
without assuming a minimizer exists. If A is empty there can be no feasible d.
The reasoning is deterministic on the enclosure event, so arbitrary adaptive
choices and stopping do not invalidate it if that event is simultaneous. QED.

An outer trace enclosure **contained in** the specification can certify a fixed
robust policy. An inner subset of trajectories contained in the specification
cannot: omitted trajectories may fail. Similarly, failure to find a feasible
candidate is not proof that every covered cell is infeasible.

## 3. A finite experimental construction of the enclosures

Here is an explicit sufficient construction, deliberately without an efficient
high-dimensional claim. Fix an r-net of N physically executable policies
d_1,...,d_N, with associated covering cells of radius at most r. Let each true
response f_j have known uniform modulus omega_j on D:

\[
 |f_j(d)-f_j(d')|\le\omega_j(\operatorname{dist}(d,d')).
\]

At a selected site i, repeated observations have mean f_j(d_i), conditional
mean-zero noise with mgf bounded by exp(lambda^2 sigma_j^2/2), and a separately
bounded systematic measurement bias b_j. Selection is predictable from past
data; resets preserve the declared physical preparation distribution. If an
experiment cannot be carried out within its own constraints, it is not an
available sample. No safe-access conclusion follows from statistical precision.

For K=m+1 quantities and n>=1 observations at each site use

\[
 r_j(n)=\sigma_j\sqrt{\frac2n
 \log\frac{\pi^2 NK n^2}{3\alpha}},\qquad
 w_j(n)=r_j(n)+b_j.
\]

For a fixed i,j,n the conditional sub-Gaussian martingale bound has two-sided
failure at most 6alpha/(pi^2 NK n^2). Summing over i,j and all n gives at most
alpha. Predictable allocation permits the same exponential-supermartingale
argument on each site's sampled subsequence (including a stopped/infinite
subsequence by optional-stopping bounds). Therefore all the intervals
[sample_mean-w_j,sample_mean+w_j] hold simultaneously at all attained counts.
This conservative union-bound construction is weaker than modern confidence
sequences but completely explicit. See [Howard et al.](https://arxiv.org/abs/1810.08240).

Expand the site's interval by omega_j(r) to enclose its entire cell. Bounds from
other sites may be intersected using the maximum site-to-cell distance.
Contradictory intervals flag failure of at least one assumed bound, measurement
condition or implementation; they do not by themselves identify a new law.

The countable-net extension predeclares all possible sites and assigns summable
weights to site/quantity/sample-count triples. A data-dependent representation
cannot create uncharged confidence claims from old pointwise intervals.

This construction estimates the *physical response directly* if the observations
are physical and the observation assumptions hold. A simulator can replace or
augment those observations only with an independently justified discrepancy
bound. A numerical residual certificate contributes numerical error; a fitted
training residual does not supply that discrepancy bound.

## 4. What finite termination really costs

Small error alone does not guarantee a budgeted constrained optimum. The geometry
of feasibility enters twice. Suppose site intervals have halfwidth at most w_j
and contain their true values. Let

\[
 s_j=2w_j+\omega_j(r),\quad
 F_s=\{d:g_j(d)\le s_j\},\quad
 \kappa(s)=V^*-\inf_{d\in F_s}J(d)\ge0.
\]

Assume a reachable near-optimal policy d_eta exists with
J(d_eta)<=V*+eta and g_j(d_eta)<=-rho_j, where rho_j>=s_j. Assume further a
known upper bound on kappa(s), or accept that no quantitative stopping rate has
yet been obtained.

The endpoint convention matters: a site's interval is [mu_ij-w_j,mu_ij+w_j];
the candidate uses the **point** upper endpoint mu_ij+w_j. Only the cell lower
endpoint is expanded to mu_ij-w_j-omega_j(r). Containment implies that a point
upper endpoint can exceed its true value by 2w_j, not merely w_j. Requiring an
entire candidate cell to be feasible would be a different, stronger test.

**Finite-cover bound.** The point nearest d_eta is certifiably feasible, and the
certificate gap satisfies

\[
 U_{\rm best}-\underline V
 \le \eta+4w_0+2\omega_0(r)+\kappa(s).
\]

**Proof.** For a site i within r of d_eta, its constraint upper endpoint is at
most g_j(d_eta)+omega_j(r)+2w_j<=0. Its objective upper endpoint is at most
V*+eta+omega_0(r)+2w_0. For any surviving cell centered at site k, its lower
constraint endpoint is site_lower_j-omega_j(r)<=0. Thus
g_j(d_k)<=site_lower_j+2w_j<=s_j, so d_k lies in F_s. The objective lower cell
endpoint is at least J(d_k)-2w_0-omega_0(r), hence at least
V*-kappa(s)-2w_0-omega_0(r). Taking minima and subtracting proves the claim.
Intersections with additional valid bounds can only improve the gap. QED.

This gives a finite constructive scheme: choose r and target sample widths so
the displayed bound is <=epsilon, sample the finite net to those widths,
construct the cells, and run the deterministic certificate. Adaptive allocation
may save experiments, but cannot claim a better worst-case rate without proof.
The required sample count n_j is computable by increasing n until
sigma_j sqrt((2/n)log(pi^2 NK n^2/(3alpha))) <= w_j-b_j. If the desired width
is below the bias floor, this construction cannot achieve it.

For D a p-dimensional cube, a uniform infinity-norm grid uses at most
(ceil(1/(2r))+1)^p sites. Sampling all K outputs jointly takes at most
N max_j n_j executions; separately measured outputs may take N sum_j n_j.
Charge each policy's preparation, reset, sensing, computation and disposal cost,
plus construction/checking and any unsuccessful search. This is exponential in
the policy dimension; it is a baseline construction, not mastery of matter.

For an empty feasible set with a uniform violation margin
min_d max_j g_j(d) >=rho>0, cell bounds eventually exclude every cell when
2w_j+omega_j(r)<rho for all j. Exact boundary cases can remain unresolved.

Compactness and continuity make kappa(s) tend to zero as s tends to zero when F
is nonempty, by taking a convergent subsequence of relaxed minimizers. They do
not supply a useful computable rate from a Lipschitz constant alone. One useful
additional condition is a feasibility error bound

\[
 \operatorname{dist}(d,F)\le H\max_j[g_j(d)]_+.
\]

For Lipschitz J, this gives kappa(s)<=L_J H max_j s_j: choose a nearby feasible
point and bound its objective difference. Such bounds are classical for linear
inequalities; see [Guler, Hoffman and Rothblum](https://doi.org/10.1137/S0895479892237744).
They must be proved for the nonlinear physical policy class, not assumed because
the optimizer converges on examples.

Two exact diagnostic counterexamples explain these obligations. First, with
g(x)=x^2 on [-1,1], only x=0 is feasible. Nonzero two-sided measurement error
need never certify this boundary point, even though the function is smooth.
Second, with J(x)=-x and g(x)=a*x on [0,1], a>0, the optimum feasible value is
0 but allowing violation s>=a makes the relaxed value -1. Arbitrarily small a
defeats a uniform value-sensitivity rate even though all these functions are
1-Lipschitz. These are different failures: lack of interior margin and poor
constraint conditioning.

## 5. Where physical structure enters, rather than being named away

For a finite-horizon discrete controlled system with a *sufficient state*, suppose
the same physical dynamics obey

\[
 \|f(x,u)-f(x',u')\|\le L_x\|x-x'\|+L_u\|u-u'\|.
\]

Two length-H open-loop control words within r in the maximum control metric
then satisfy e_H<=L_x^H e_0+L_u r sum_{k=0}^{H-1}L_x^k, by induction. Lipschitz
trace objectives give policy-response moduli from these state bounds. These are
bounds on the declared dynamics; physical model discrepancy adds a separate term
at each step. Unknown history and omitted couplings invalidate the premise.

A threshold failure probability needs more: an indicator is not Lipschitz.
Suppose executions of policies d,d' can be coupled so their trace distance is at
most e except on an event of probability q. Let failure be h(trace)>0, with h
L_h-Lipschitz. Whenever the failure indicators differ outside the exceptional
event, |h(trace_d)|<=L_h e. Consequently

\[
 |p_{\rm fail}(d)-p_{\rm fail}(d')|
 \le q+\Pr\{|h(\text{trace}_d)|\le L_h e\}.
\]

This follows by bounding the difference of expectations by the probability of
indicator disagreement. A uniform anti-concentration bound on the right converts
trajectory control into a risk modulus. Without it, even stable, smooth dynamics
can have a discontinuous deterministic failure probability at a threshold. This
additional premise must be established for reliability claims; it does not follow
from a numerical integration certificate or a known state Lipschitz constant.

For coupled reduced components a legitimate interface theorem needs component
state sensitivity, interface-input sensitivity, *and* coupling gains. If those
give a nonnegative matrix A and local error vector b with
e_{t+1}<=A e_t+b, then e_H<=A^H e_0+sum_{k=0}^{H-1}A^k b. A Lipschitz constant
for the coupling alone is insufficient. Stability/small gain can prevent growth,
but finite-horizon composition does not require asymptotic stability.

For example a mass-action reaction network on a proven bounded concentration
domain has polynomial vector fields whose derivative bounds can be computed.
A thermal transport network on a bounded state domain has explicit coupling
gains. These are distinct conditional members of the same theorem, not a claim
that their true unknown rates, valid domains, retained variables, fabrication
routes or measurement biases have already been discovered. Phase transitions,
rare events, hidden microstructure and chemistry-dependent preparation can defeat
the assumed cheap representation. That is the physical research obligation.

## 6. What the existing literature does and does not remove

[Huang, Chen and Preskill](https://arxiv.org/abs/2210.14894) establish efficient
prediction of local outputs with small average error for specified input
distributions. The uniform intervention moduli/enclosures needed here are
additional obligations. [SafeOpt](https://proceedings.mlr.press/v37/sui15.html)
already gives safe adaptive optimization within its assumptions and reachable
safe set. The concentration/optimization connection is therefore established
mathematics, not Spectra's invention. Physical executable-policy decoding,
history closure, cross-domain representation discovery and affordable physical
applicability bounds remain outside this reduction.

With no quantitative regularity or accessible distinguishing experiments, a
finite-query algorithm cannot rule out a response bump supported away from every
queried policy. With M disjoint possible bumps and q queries, transcript total
variation from the zero response can be at most q/M, yielding minimax testing
error at least (1-q/M)/2. [The self-contained obstruction](obstructions.md)
explains why additional experiments alone do not solve unrestricted extrapolation.

## 7. Disposition and next scientific obligation

The result here is a general **conditional constructive connection**, and a
sharper identification of the feasibility/conditioning and physical-coverage
obligations. The executable policy-cover checker validates its arithmetic and
refusal logic on mathematical examples. Neither is evidence of general physical
synthesis, discovery of new laws, or autonomous method acquisition.

The central unresolved target is to obtain affordable, intervention-valid
enclosures by discovering a sufficient representation and its physical
applicability on independently specified chemical, materials and engineering
families. A proposed extension must replace a costly premise above with a
constructive, evidenced result—for example a stable transfer relation that
reduces the covering/experiment cost while retaining preparation and coupling
validity. Simply invoking a confidence oracle, reducing a fixed Taylor cost, or
renaming this grid construction as an intelligence would leave that target open.
