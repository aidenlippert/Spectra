# Targets 3 and 6: spacetime connectivity and thermal passivity

## 3. Constructive control of macroscopic spacetime connectivity

Gao, Jafferis, and Wall (GJW) show a precise but narrow positive result: coupling the two boundaries of an eternal BTZ black hole with a double-trace interaction produces negative averaged null energy, whose semiclassical backreaction makes the Einstein–Rosen bridge traversable. Their abstract explicitly says that the resulting wormhole cannot be used to violate causality ([GJW, arXiv:1608.05687](https://arxiv.org/abs/1608.05687)). This is a source–geometry construction in a specially prepared asymptotic setting, not a laboratory recipe for creating a macroscopic throat, an asymptotically flat shortcut, or faster-than-light transport.

The relevant engineering theorem must solve the source problem and the metric problem together. For a declared physically preparable state/control class \(\mathfrak S_B(g)\), one needs a finite-resource solution

\[
G_{\mu\nu}+\Lambda g_{\mu\nu}=8\pi Gc^{-4}
\bigl(T^{\rm classical}[u]+\langle T_{\mu\nu}\rangle^{\rm ren}_{\psi}\bigr),
\]

with a macroscopic traversable region, bounded tidal stress, specified transit fidelity and lifetime, and stability under perturbations. The state \(\psi\), controls \(u\), and initial data must be preparable; choosing a metric and assigning it the required stress tensor is not a construction.

Quantum energy inequalities (QEIs) are the sharp relevant restriction. Fewster’s review emphasizes that QFT violates classical pointwise energy conditions while satisfying averaged remnants ([Fewster, arXiv:1208.5399](https://arxiv.org/abs/1208.5399)). Fewster and Kontou derive state-dependent lower bounds for averaged effective energy density, including thermal-state dependence ([Fewster–Kontou, arXiv:1809.05047](https://arxiv.org/abs/1809.05047)). Thus negative energy is allowed in restricted spacetime averages, but is not an unlimited engineering resource. The exact sampling function, field, state, geometry, and preparation cost belong in the bound. Global-causality and energy assumptions in topological-censorship arguments must be stated explicitly; evading one assumption changes applicability rather than disproving the theorem.

The smallest useful feasibility lemma is: for a fixed QFT, control budget \(B\), and source class \(\mathfrak S_B\), characterize whether the attainable averaged null-energy profiles intersect the profiles required by a chosen throat after self-consistent backreaction. The first decisive test is a finite-resource source–geometry solution with a predicted transit signal and perturbation response distinguishable from lensing, a waveguide, or an ordinary causal channel. Stop if the proposal requires an infinite pulse, an unpreparable state, a prescribed metric independent of source dynamics, or causal behavior that exists only in the idealized boundary setup.

## 6. Resource-complete violation of thermal passivity

Passivity has two related but distinct formulations. For a system with Hamiltonian \(H\), ordinary passivity says that no cyclic external drive represented by an arbitrary unitary \(U\) can lower the system’s mean energy: \(\operatorname{Tr}(H\rho)\leq\operatorname{Tr}(HU\rho U^\dagger)\). Complete passivity requires this for every finite tensor power. Pusz and Woronowicz prove that, under their C*-dynamical assumptions, the completely passive states are ground states and KMS states ([Pusz–Woronowicz, DOI](https://doi.org/10.1007/BF01614224)). This result concerns cyclic unitaries on the tested system; it should not be conflated with the separate thermal-operations model, where system, bath, battery, and catalyst undergo global energy-conserving dynamics.

A one-reservoir heat engine must therefore specify its global process. Let \(S\) be the working system, \(B\) an equilibrium bath, \(W\) a work store, and \(C\) controllers/catalysts/records. For a closed global energy-conserving unitary, total energy and entropy are both constant, so the total nonequilibrium free energy is exactly conserved:

\[
\Delta F_\beta(SBWC)=0,
\qquad
\Delta F_\beta=\Delta\langle H\rangle-\beta^{-1}\Delta S,
\]

with the total Hamiltonian including interaction terms. At decoupled endpoints,

\[
F_\beta^{\rm tot}=\sum_i F_\beta(i)+\beta^{-1} I_{\rm tot},
\]

so \(\Delta F_\beta(W)=-\sum_{i\neq W}\Delta F_\beta(i)-\beta^{-1}\Delta I_{\rm tot}\) (with any endpoint interaction-energy change included in the corresponding total-energy term). Thus a claimed work gain is valid only after accounting for work-store entropy, bath and controller free energies, interaction energy, and all correlations. “The catalyst returned” must include its marginal, correlations, and any decoupling cost. Resource-theoretic treatments explicitly require returned devices in catalytic thermal operations ([Brandão et al., DOI](https://doi.org/10.1038/ncomms6407)); later work shows that approximate or multicopy catalysis can activate transformations, so catalyst dimension, error, copy count, and cumulative degradation must be bounded ([Lipka-Bartosik–Skrzypczyk, PRX](https://doi.org/10.1103/PhysRevX.11.011061)).

For an initially Gibbs bath at the reference \(\beta\), \(F_\beta(B)\) is minimal: \(\Delta F_\beta(B)=\beta^{-1}D(\rho_B'\|\tau_\beta)\geq0\). Cooling that initially Gibbs bath raises its nonequilibrium free energy relative to the initial temperature. Its loss of energy alone cannot account for cyclic entropy-free work: the entropy and free-energy balance still needs another changing resource. A bath initially out of equilibrium is a different case and may supply consumable free energy. The infinite-bath idealization is different: the bath can remain effectively Gibbs while exporting entropy, but then the global bath limit and heat/work bookkeeping must be stated. If controllers are restored and the work-store entropy is unchanged, the endpoint identity gives \(W\leq0\) when the bath and all other non-work resources have nonnegative free-energy change and correlations are not used as a hidden resource. Strong coupling does not create a loophole: changes in \(H_{SB}\), correlations, controller energy, and record reset must be included.

The smallest decisive lemma is a full cyclic map on \(S,B,W,C\) proving positive \(\Delta E_W\) while bounding every non-work free-energy and correlation term, including the work store’s entropy. The decisive experiment is repeated operation with calibrated bath temperature and independent calorimetry of controller power, interaction energy, catalyst drift, record reset, and bath state. Stop when gain is explained by hidden drive work, nonthermal bath consumption, catalyst degradation, initial coherence, residual correlations, or interaction-energy bookkeeping.

### Sources

- https://arxiv.org/abs/1608.05687
- https://arxiv.org/abs/1208.5399
- https://arxiv.org/abs/1809.05047
- https://doi.org/10.1007/BF01614224
- https://doi.org/10.1038/ncomms6407
- https://doi.org/10.1103/PhysRevX.11.011061
