# Weak-coupling cluster expansions: a plausible scaling condition

## Closest constructive theorem

Bravyi, DiVincenzo, and Loss study (H=H_0+\epsilon V), where (H_0) is a product of single-qubit terms and (V) is a bounded two-qubit interaction on a bounded-degree graph. Their theorem gives a polynomial-time algorithm for ground-state energy and correlation functions when (|\epsilon|) is below a threshold depending on the gap of (H_0), graph degree, and interaction norms. The algorithm uses a generalized Kirkwood–Thomas ansatz and perturbative expansions ([arXiv:0707.1894](https://arxiv.org/abs/0707.1894); published DOI [10.1007/s00220-008-0574-6](https://doi.org/10.1007/s00220-008-0574-6)).

This is the strongest directly relevant precedent for the proposed idea: a checkable small-coupling/gap regime can turn an apparently many-body ground-state calculation into polynomial work. It is an energy/state approximation algorithm, however. The paper does not establish a rational fermionic SOS dual certificate, a CAR-ideal quotient construction, or an independently replayable lower-bound bracket in the Spectra sense.

## What the desired (O(\log(N/\varepsilon))) cluster argument would require

Suppose a separate linked-cluster theorem proves that the *total* contribution of connected clusters larger than (k) is bounded by

\[
  |R_k|\le N C e^{-\mu k},
\]

for a bounded-degree graph and constants (C,\mu>0) uniform in (N). Then choosing

\[
  k\ge \mu^{-1}\log(NC/\varepsilon)
\]

gives an absolute total-energy error at most (\varepsilon), rather than the common mistake of assigning (\varepsilon/N) independently to every local term. If the number of connected clusters of size (j) containing a root is at most (B^j), enumerating all clusters through (k) costs at most

\[
  N\sum_{j\le k}B^j = N\,B^{O(\log(N/\varepsilon))},
\]

which is polynomial in (N/\varepsilon) for fixed (B,\mu). This is a useful theorem target, not a consequence of Bravyi–DiVincenzo–Loss as cited in the abstract.

To make it a certificate result, each retained cluster calculation also needs a rigorously bounded local energy/error contribution, and the discarded tail must be bounded in the same norm used by the global lower bound. A convergent perturbation series or an approximate ground state by itself does not give a dual SOS lower certificate.

## Model limitations

The Bravyi–DiVincenzo–Loss assumptions are product (H_0), bounded-degree graph, bounded two-body norms, and sufficiently weak coupling. Ab-initio fermionic Hamiltonians generally have long-range Coulomb matrix elements, basis-dependent graph degree, particle-number constraints, and possible small gaps or near-degeneracies. Mapping fermions to qubits can also change locality and coefficient norms. Any application to H4/H6/H8/H10 must therefore report the actual gap, degree/norm bound, coupling ratio, and perturbative radius rather than labeling the family “weakly interacting.”

The impurity-model result of Bravyi and Gosset is a separate useful precedent: for a quadratic fermionic bath plus an (O(1))-mode interacting impurity, they obtain quasi-polynomial time in system size and exponential dependence on requested bits ([arXiv:1609.00735](https://arxiv.org/abs/1609.00735)). It does not cover general molecular active spaces, but it suggests testing certificate transfer first on impurity-like decompositions.

## Falsifiable structural experiment

For a family with a deliberately controlled product reference and bounded-degree perturbation, compute cluster contributions through order (k), independently bound the remainder, and compare the resulting interval with exact diagonalization. Sweep (N), coupling strength, and gap. Success requires the predicted (N e^{-\mu k}) tail to remain valid and the number of retained clusters to follow the polynomial bound. Only after that should the same machinery be tested on fermionic Hamiltonians; the required bridge is a theorem translating cluster truncation into a certified SOS dual certificate and CAR-ideal verification.
