# Positive-amplitude/Jastrow route for stoquastic fixed-(N) chains

## Exact identity

Let (H) be real symmetric in an occupation basis, with (H_{st}\le 0) for (s\ne t), and let (psi_	heta(s)>0). Define

\[
 E_\theta(s)=\frac{(H\psi_\theta)(s)}{\psi_\theta(s)}
 =H_{ss}+\sum_{t\ne s}H_{st}\frac{\psi_\theta(t)}{\psi_\theta(s)}.
\]

For every off-diagonal edge (\{s,t\}),

\[
 H=\sum_{\{s,t\}}(-H_{st})\left(\sqrt{\frac{\psi_t}{\psi_s}}|s\rangle-\sqrt{\frac{\psi_s}{\psi_t}}|t\rangle\right)\left(\cdots\right)^\dagger
 +\operatorname{diag}(E_\theta).
\]

The edge terms are PSD, hence (H\succeq \operatorname{diag}(E_\theta)) and

\[
 E_0(H)\ge \min_{s\in\mathcal S_N}E_\theta(s).
\]

For a strictly positive finite-range Jastrow ansatz

\[
 \psi_\theta(s)=\exp\!\left(\sum_a\theta_a f_a(s_{a:a+r})\right),
\]

each (E_\theta(s)) is a sum of local terms involving exponentials of linear functions of (	heta). It is concave in (	heta), because each off-diagonal coefficient is nonpositive. Therefore

\[
\max_{\theta\in[-B,B]^p}\min_s E_\theta(s)
\]

is a convex optimization problem (maximize a concave objective over a convex box). This is a concrete lower-bound family, closely related to the stochastic-matrix-form/Rokhsar–Kivelson construction ([Castelnovo et al., 2005](https://doi.org/10.1016/j.aop.2005.01.006)) and to Perron–Frobenius stoquastic Hamiltonians ([Bravyi, DiVincenzo, Oliveira & Terhal](https://arxiv.org/abs/quant-ph/0606140)). It is not a general fermionic SOS certificate: positivity of amplitudes is a severe sign-structure assumption.

## Efficient separation in one dimension

Assume a 1D chain with finite-range hopping/diagonal terms, Jastrow range (r), and a fixed particle-number sector. The ratio (psi(t)/psi(s)) changes only windows intersecting the local hop, so (E_\theta(s)) is a finite-range classical factor sum. Its exact minimum can be found by a junction-tree/dynamic-programming pass whose state is the last (O(r+\text{interaction range})) occupation bits plus the accumulated particle count. For fixed range and local alphabet, the arithmetic state count is polynomial in (N) (typically (O(N^2)) with a total-(N) counter), rather than exponential in the Hilbert-space dimension.

With (|\theta_a|\le B), coefficient bounds (|H_{st}|\le J), and (q) bits of interval arithmetic, every exponential ratio has a computable enclosure. The separation oracle should return a certified lower enclosure for the minimum and a witness configuration; the optimization must then include an outer precision/error budget. Exact symbolic separation is generally unavailable because (exp(\theta)) is transcendental, so “exact” means interval-certified to a requested tolerance.

For a 1D finite-range chain, the variational upper energy

\[
 U_\theta=\frac{\langle\psi_\theta|H|\psi_\theta\rangle}{\langle\psi_\theta|\psi_\theta\rangle}
\]

can be contracted by a finite transfer matrix, with a particle-number fugacity or coefficient-extraction DP for the canonical sector. Interval transfer-matrix arithmetic gives a certified upper enclosure. Thus the route can produce a two-sided bracket when both contractions and the lower separation oracle are implemented rigorously.

## What is and is not proved

The structural condition is efficiently checkable: stoquastic signs, bounded interaction/Jastrow range, bounded local degree, and a 1D (or bounded-treewidth) factor graph. The lower-bound identity is exact for every positive (psi). The DP complexity claim follows from finite-range factorization, and the transfer-matrix upper bound is standard finite-state contraction.

What is *not* automatic is accuracy. There is no theorem that a bounded-range positive Jastrow family approximates the ground state of every stoquastic Hamiltonian to chemical accuracy. Frustration, criticality, near-degeneracy, long-range Coulomb terms, and sign-changing fermionic ground states can make the bracket wide. Perron–Frobenius gives positivity only in a connected stoquastic basis; it does not imply a short-range Jastrow logarithm.

The practical falsification ladder is H4/H6/H8/H10 in a genuinely stoquastic fixed-(N) 1D family: optimize the lower bound with DP separation, contract the upper bound by transfer matrices, and report width versus range (r), parameter bound (B), precision, and (N). A successful scaling claim requires width decay on held-out couplings/geometries; polynomial separation alone is not an accuracy theorem.

## Existing implementation and remaining extension

The repository already has `stoquastic_chain_certificate.py`: a restricted
periodic spin-chain implementation with one positive domain-wall parameter,
exact local-energy min/max DP, and known parent controls through 256 sites.
Those controls are prior work, not a new result of this continuation. They
are not number-conserving fermionic molecular Hamiltonians.

The concrete proposed extension is local hopping in a fixed-N fermionic
chain, a particle-count DP, and a larger finite-range log-amplitude family.
An efficient sign check here assumes an explicit bounded-range occupation
basis representation; discovering a stoquastic basis is a separate problem.
General bounded-treewidth statements must apply to the expanded local-energy
factor graph after amplitude ratios, not merely to the original interaction
graph. The existing spin-chain code does not establish these extensions.

The two primary-source abstracts and bibliographic identities above were
checked directly on arXiv during the root audit. The DP/optimization proposal
is a derivation here, not a claim that either cited paper supplies an
accurate general chemistry certificate theorem.
