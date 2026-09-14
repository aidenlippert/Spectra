# Compact local dynamics: result receipt

A conditional mathematical certificate now handles an interacting XX+ZZ chain without constructing its full Pauli space. It does not demonstrate learned scientific strategy or broad compounding capability.

Read [the theorem and constructive numerical method](../../research/v3/local_pauli_theorem.md).

| Qubits | Hamiltonian terms | Retained Pauli observables | Uniform physical expectation error bound |
|---|---:|---:|---:|
| 3 | 7 | 30 | 0.000238912440 |
| 6 | 16 | 98 | 0.000238912440 |
| 12 | 34 | 109 | 0.000238912440 |
| 24 | 70 | 109 | 0.000238912440 |
| 48 | 142 | 109 | 0.000238912440 |

The certificate applies at every time and to every density state for the supplied one- and two-site Hamiltonian and known local depolarization rate. The tested parameters are XX=1, ZZ=1/5, on-site Z=1/2, and depolarization gamma=300. The strong damping assumption is essential.

The retained set contains every word reached by four commutators from the local seed. Exact rational bounds combine a short-time Dyson tail with long-time contraction in the Pauli coefficient l1 norm. This norm directly bounds observable error without an exponential Hilbert-space conversion.

Rational uniformization supplies numerical predictions with an additional declared error below 1e-10. Certificate construction and verification use local frontier operations, not a full 4^N matrix. Dense full-space checks are limited to three qubits.

A product-state counterexample makes the physical interaction matter: at time 1/1000, damping-only evolution predicts zero, while the full expectation is certified above 0.00232. This exceeds the requested 0.001 tolerance. It is a sanity check against one simple baseline, not an advantage over strong stateful algorithms.

The full project now passes **65 tests**. These include missing-damping rejection, false compactness refusals, certificate tampering, independent dense dynamics, exact numerical enclosures, and the damping-only counterexample.

The Hamiltonian and noise parameters are supplied. The procedure for selecting the neighborhood is programmed. No acquisition-cost advantage, learned representation-selection rule, real apparatus validation or quantum-memory protection is claimed.

- [Exact certificates, prediction intervals and counterexample states](locality_results.json)
- [Adversarial mathematical review](../../research/v3/locality_adversarial_review.md)
- [Hilbert–Schmidt norm obstruction](../../research/v3/physical_norm_obstruction.md)
- [Primary-source context](../../research/v3/locality_sources.md)
- [Active research status](../../research/ACTIVE_STATUS.md)
