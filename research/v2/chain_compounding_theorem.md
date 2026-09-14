# Exact compounding for an `r`-factor response chain

This is a restricted statistical theorem extending the two-factor construction. It demonstrates positive marginal value from sequentially acquired response laws; it is not a claim of unrestricted scientific discovery or algorithmic superiority over an equally informed Bayesian designer.

## Finite physical model

The target is an `r`-qubit state

`rho_{s,a_1,...,a_r} = [I + s v (⊗_{j=1}^r P_{a_j})]/2^r`,

where `s∈{−1,+1}`, each `a_j∈{0,1}`, `P_0=X`, `P_1=Z`, and known `0<v≤1`. Hidden bits are independent and uniform. The laboratory may choose any product-Pauli action `(i_1,...,i_r)∈{0,1}^r`; its plus outcome has probability

`(1+s v ∏_j 1[i_j=a_j])/2`.

Operational preparation is assumed to be available as a finite source: sample a uniformly random eigenvalue string, prepare each qubit in the requested Pauli eigenstate, and classically correlate product eigenvalues to produce the displayed parity term. Prefix-reference sources are separately preparable sources retaining the first `j`-body parity correlation; they are not obtained by tracing later coordinates, which would erase the correlation. Preparation and every measurement shot are counted. This explicit source assumption means the theorem studies transfer within the family, not derivation of the family from a microscopic theory.

## Risk with `k` learned laws

Suppose the first `k` response laws are known at the target knobs, while the remaining `r−k` axes are unknown. Under the uniform prior choose arbitrary actions on unknown coordinates, or choose them uniformly at random for a pointwise guarantee over every fixed hidden mask. The match probability is `2^{-(r-k)}`. A matched product measurement gives optimal sign error `(1-v)/2`; a mismatch gives `1/2`. Therefore the exact Bayes risk is

`R_k = 1/2 - v / 2^{r-k+1}`.

This is optimal over the complete `2^r` product-Pauli action menu: conditionally on a mismatch, the outcome law is identical for `s=±1`; conditionally on a match, the likelihood-ratio decision is optimal. The marginal value of learning law `k+1` is

`R_k-R_{k+1}=v/2^{r-k+1}`,

which doubles at every successive acquisition. The positive discrete complementarity is exact in this restricted model. A conventional Bayesian designer with the same observations can attain the same risks.

## Sequential acquisition protocol

Only the first direct calibration source is assumed. Each law is a public GF(2) response `a_j(x_j)=m_j·x_j mod 2`, with a `d`-dimensional mask `m_j`; `d` linearly independent contexts identify it noiselessly. Under BSC(`eta`) calibration noise, repeat each context and majority-decode.

To acquire law `j>1`, use the separately preparable known-positive prefix-`j` reference source whose first `j−1` axes are selected using already learned laws, and set the final candidate action to axis `0`. Conditional on prefix correctness, the plus probability is `(1+v)/2` when `a_j=0` and `1/2` when `a_j=1`. Apply the same random postprocessing: plus outputs bit `1` with probability `q=v/(2+v)`, minus always outputs `1`. Exactly,

`Pr(1|a_j=0)=[1-v+q(1+v)]/2=1/(2+v)`,

`Pr(1|a_j=1)=(1+q)/2=(1+v)/(2+v)`,

so each acquired law is a BSC with crossover `1/(2+v)` and bias `v/[2(2+v)]`. The first direct calibration law has its stated BSC noise. A sequential union bound over `r d` mask bits gives total acquisition failure at most `δ` using

`O(r d log(rd/δ)/v²)`

shots, since the compiled-label bias is `v/[2(2+v)]=Θ(v)` for `0<v≤1`. Memory is `O(rd)` bits plus the current reference state; each prefix protocol uses `O(r)` qubits. If any prior law is wrong, later conditional guarantees fail; allocating per-stage failure `δ/r` and conditioning on earlier successes yields the stated global bound. The earliest marginal gains are `v/2^{r-k+1}`, so proving every noisy marginal positive requires error below that exponentially small scale; fixed `δ` alone does not certify all margins as `r` grows.

## Scope

The increasing marginal is a property of the supplied factorization, uniform hidden axes, product-Pauli menu, and parity-correlated source. Correlated axes, arbitrary POVMs, adaptive multi-shot target policies, or a source-preparation cost that scales differently can change the result. The learner infers finite GF(2) labels from observations and reuses them on unseen knob combinations; it has not discovered that this factorization is the correct law outside the admitted class. The theorem’s scientific contribution is therefore a precise benchmark for compounding transfer, not evidence of universal physical synthesis.
