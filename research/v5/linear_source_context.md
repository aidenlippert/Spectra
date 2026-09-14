# Linear Gaussian thermal source: primary context and scope

The coupled source in `experiments/v5_core.py` is an exactly solvable overdamped Ornstein–Uhlenbeck process. No physical apparatus was built.

## Primary sources

1. **Uhlenbeck and Ornstein, “On the Theory of the Brownian Motion,” Physical Review 36, 823–841 (1930).** [APS DOI](https://doi.org/10.1103/PhysRev.36.823). The original paper derives Gaussian transition laws and exact time-dependent means and variances for linear Brownian motion, including the harmonically bound/overdamped regime. This is the historical primary source for the OU/Langevin calculation used here.

2. **Särkkä and Solin, Applied Stochastic Differential Equations (2019), Chapter 6.** The transition mean and covariance of a general linear SDE are derived in equations 6.1–6.7 of the [authors' book PDF](https://users.aalto.fi/~asolin/sde-book/sde-book.pdf); see also the [publisher's chapter](https://doi.org/10.1017/9781108186735.007). Specializing these formulas to symmetric drift C⁻¹ and diffusion √2 I gives the expressions below. An earlier draft's unverified Scientific Reports citation was removed during source audit and is not evidence.

3. **Widrow et al., “Adaptive Noise Cancelling: Principles and Applications,” Proceedings of the IEEE 63 (1975), 1692–1716.** [Bibliographic primary record](https://doi.org/10.1109/PROC.1975.10036). This establishes adaptive cancellation as a signal-processing method using a reference correlated with interference; it is relevant context for optional experimental noise-nulling, not evidence that thermal sampling or Hamiltonian identification is free. Any physical implementation still needs a reference channel, calibration, actuator bandwidth, and a stability/error budget.

## Exact calculation and resource meaning

For

\[
dq_t=-C^{-1}q_t\,dt+\sqrt 2\,dW_t,
\qquad C=I+\sum_i\lambda_i u_i u_i^T,
\]

with positive-definite `C`, the stationary Gaussian has covariance `C`, since the Lyapunov equation is
`C^{-1}C+C C^{-1}=2I`. Conditional on `q_0`,

\[
\mathbb E[q_t\mid q_0]=e^{-C^{-1}t}q_0,
\qquad
\operatorname{Cov}(q_t\mid q_0)=C(I-e^{-2C^{-1}t}).
\]

Taking `T=30 lambda_max(C)` makes the slowest mean mode at most `e^{-30}` in the spectral norm, under the stated initialization and exact coefficients. It does not make the output an exact independent draw from the equilibrium Gaussian at finite time: the mean retains an initialization-dependent term and the covariance differs from `C` by an exponentially small but nonzero matrix. Preparation time, coefficient precision, numerical integration/sampling, and readout costs must remain in the ledger.

The harmonic energy `H(q)=1/2 q^T C^{-1}q` has density proportional to `exp(-H)` at `k_B T=1`, so the target Gaussian is the Gibbs law for this classical quadratic model. The calculation connects the statistical family to coupled harmonic matter; it does not derive an apparatus, arbitrary prefix-source availability, a reset mechanism, or an unknown-source identification procedure. The finite-time bound and resource assumptions are in [PROOF.md](PROOF.md). Sampling an analytic scalar marginal with NumPy is not a rigorous TV guarantee for its floating-point random generator.
