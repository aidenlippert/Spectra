# Representation revision for omitted GF(2) interactions

The current physical benchmark exposes a finite linear calibration law. The
smallest rigorous extension is to let the hidden source bit obey

\[
f(x)=m\cdot x+\sum_{i<j}q_{ij}x_ix_j\pmod 2,
\]

with (x\in\{0,1\}^d), while initially fitting only the linear representation.
This is a representation-revision test, not an unrestricted law-discovery
claim.

## Constructive identification

Query the zero vector, every unit vector (e_i), and every pair vector
(e_i+e_j). In the noiseless case,

\[
m_i=f(e_i)+f(0),
\qquad
q_{ij}=f(e_i+e_j)+f(e_i)+f(e_j)+f(0).
\]

Thus (d+\binom d2+1) distinct contexts identify all coefficients in the
declared degree-two class. The learner need not receive (m) or (q); it
estimates each bit from repeated observations and constructs the revised
quadratic representation only after a linear residual witness appears.

With independent binary symmetric noise η<1/2, estimate each queried context
by (n) shots and majority vote. Hoeffding gives per-context error at most

\[
\exp[-2n(1/2-\eta)^2].
\]

Taking

\[
n\geq\frac{\log((d+\binom d2+1)/\delta)}{2(1/2-\eta)^2}
\]

gives simultaneous coefficient recovery with probability at least (1-δ),
subject to independent shots and the stated BSC model. If the source returns a
real-valued statistic instead, replace this bound with its declared bounded-noise
concentration inequality.

## Falsifying the old representation

Fit the best linear law on (0,e_i). For any pair (i,j), the parity residual

\[
r_{ij}=f(e_i+e_j)+f(e_i)+f(e_j)+f(0)
\]

is exactly (q_{ij}). A confidence interval excluding zero is a constructive
witness that the linear representation is inadequate. The policy should then
acquire the corresponding interaction feature and retest held-out contexts.
Do not add a feature merely because it improves training error under noise.

The minimal independent test is: discover one nonzero pair interaction in stage
one, discover a second interaction in stage two using the revised representation,
then evaluate both on held-out contexts of Hamming weight at least three and on
new coefficient settings. The learner sees only queried outcomes and residuals;
the hidden coefficients are used for evaluation only.

## Transfer and the limit of the claim

For a degree-two law, the revised coefficients transfer to every held-out input,
including Hamming weight ≥3, because evaluation is the same polynomial over
GF(2). A fresh reset learner must relearn all coefficients. A warm learner can
reuse the previously certified interaction features and estimate only new
coefficients or test new residuals.

This can yield a real resource reduction, but it is not automatically
superadditive compounding. Report the number of new context queries and the
number of retained certified features. Compare reset, retrieval-only, warm
revision, and a full quadratic baseline under the same total shots. A result
where warm revision merely saves repeated support search is reuse, not evidence
that the system has learned a deeper scientific operation.

No finite degree-two experiment can rule out cubic structure. Every cubic law
can agree with a chosen set of weight-≤2 calibration queries. To test the
declared degree bound, include held-out weight-three contexts. Failure there
must produce abstention or an explicit degree-expansion request. Agreement on
the selected probes certifies only the finite declared class and noise model.

The bridge to the physical module is direct: its masked calibration callback
supplies (f(x)) observations, while representation revision adds pair probes
and retains the same finite-shot confidence accounting. The current model's
linear law is therefore a baseline representation, not a hidden oracle for the
quadratic extension.
