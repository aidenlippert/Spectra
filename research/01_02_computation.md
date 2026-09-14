# Targets 1–2: hypercomputation and nonlinear quantum computation

This note separates established mathematical results from the stronger physical claims required by the two keystone targets. The central standard is operational: preparation, evolution, readout, resources, and error bounds must all be specified.

## 1. Robust physical hypercomputation

The target is a uniformly specified, physically preparable procedure that decides `HALT(M,x)` in finite observer time with bounded error, for every finite machine description. An ideal Malament–Hogarth (MH) spacetime contains a worldline on which a Turing machine experiences infinite proper time before a finite observer event. Etesi and Németi explicitly frame this as a consequence of choosing classical general relativity as the background theory, not as an engineered device: [arXiv:gr-qc/0104023](https://arxiv.org/abs/gr-qc/0104023).

The finite-approximation obstruction is precise. Suppose a proposed device has finite approximants `S_k`, and, for each input, there is an effective procedure that (i) computes the output probability of `S_k` to arbitrary requested numerical error, (ii) supplies a computable, input-dependent modulus `e(k) -> 0` for the difference between the ideal probability and `S_k`, and (iii) guarantees a fixed gap: yes instances have probability at least `2/3`, no instances at most `1/3`. Search for `k` with `e(k)<1/12`, then simulate `S_k` with error `<1/12`. The total error is `<1/6`, leaving a strict threshold gap at `1/2`; this gives an ordinary terminating algorithm for HALT, contradiction.

“Computable descriptions” by themselves are insufficient: the effective tail/modulus and termination of the probability computation are essential. A noncomputable modulus is logically one escape, but it is not a mechanism or evidence that physics supplies one. Other escape routes are noncomputable initial data, an unbounded/infinite auxiliary process, an exact singular limit with no effective approximation, or a readout with no finite reliability margin.

The output channel is therefore the decisive engineering theorem. An MH construction would need a finite preparation history; a signal that reaches the observer before the MH event; a stable encoding of “halted”; bounded backreaction and disturbance; and a survival/error analysis for the machine over unbounded proper time. A resource audit must include memory, matter, energy, acceleration, cooling, and signal attenuation. A proposed computer that fails with probability tending to one along the infinite trajectory is not robust.

Closed timelike curves illustrate why model labels matter. Aaronson and Watrous prove in [arXiv:0808.2669](https://arxiv.org/abs/0808.2669) that polynomial-size Deutsch-CTC computations have classical and quantum power exactly `PSPACE`, through a causal-consistency fixed-point condition. That is decidable complexity amplification, not a halting oracle. The unrestricted-width/length model is different: Aaronson et al. [arXiv:1609.05507](https://arxiv.org/abs/1609.05507), version 2 (2024), identifies the computability level with `Delta_2`, languages Turing-reducible to HALT, while also noting that version 1 contained an erroneous proof later corrected. Neither result supplies a physically realized CTC or validates its fixed-point resource.

Cheap falsifiers for an alleged hypercomputer are: a computable finite approximation with an effective modulus; ordinary simulation substituted for the output channel; an unbounded resource hidden in the trajectory; or a readout margin that shrinks faster than all certified experimental precision.

## 2. Operational nonlinear quantum complexity amplification

The useful target is not “some nonlinearity.” It is an explicit family of state transformations `Phi_n` with a complete multipartite semantics and polynomial total resources such that an exponentially small input distinction is amplified to a constant. For trace distance `D`, require, for example,

`D(rho_0,rho_1) = 2^{-poly(n)}` and `D(Phi_n^m(rho_0), Phi_n^m(rho_1)) >= 1/3`, with `m=poly(n)`.

Abrams and Lloyd show in [Phys. Rev. Lett. 81, 3992 (1998)](https://doi.org/10.1103/PhysRevLett.81.3992) that suitable hypothetical nonlinear gates give polynomial-time algorithms for NP-complete and `#P` problems. This is an established conditional theorem; it does not say that every nonlinear evolution has the required gate, conditioning, precision, or fault tolerance.

For a scalar signal error `d_t`, a Lipschitz amplification step with factor `lambda` and additive uncertainty `nu` obeys

`d_{t+1} <= lambda d_t + nu`;

therefore `d_m <= lambda^m d_0 + nu (lambda^m-1)/(lambda-1)` for `lambda != 1`, and `d_m <= d_0 + m nu` for `lambda=1`. If the useful signal begins at exponentially small `d_0`, polynomial-time amplification requires a noise/fault-tolerance theorem ensuring that the second term remains below the final decision margin. A bare unstable map usually converts exponentially fine preparation errors into the same scale as the desired signal.

Success accounting is equally mandatory. If one complete run succeeds with probability at least `2/3`, majority voting over `R` independent odd repetitions has failure probability at most `exp(-R/18)` by Hoeffding. If the nonlinear procedure uses postselection with success probability `q`, `R` attempts fail with probability `(1-q)^R` and require expected `1/q` attempts; polynomial complexity requires `q >= 1/poly(n)`.

The map must preserve positivity and normalization, define evolution of entangled composites, and specify how mixed-state decompositions are treated. Weinberg-type nonlinear quantum mechanics has known causal problems: Polchinski’s analysis [Phys. Rev. Lett. 66, 397 (1991)](https://doi.org/10.1103/PhysRevLett.66.397) finds EPR/branch communication consequences. A no-signaling nonlinear extension may be possible, but it must be shown that its composite dynamics still implements the Abrams–Lloyd amplification; no-signaling alone does not imply a complexity advantage.

Cheap falsifiers are an effect explainable by ordinary state-dependent control; exponentially precise preparation; exponentially rare postselection; nonpositive or decomposition-dependent maps; or the absence of a uniform `Phi_n` family and reduction. The bounded research target is thus a consistency, resource, and fault-tolerance theorem for one sharply specified nonlinear operation, followed by a discriminating experiment—not an inference from an anomalous signal to a complexity collapse.
