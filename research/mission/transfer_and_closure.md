# Transfer of physical abstractions under interventions and coupling

## Scope

This note addresses one narrow part of the mission: when a learned mechanism or physical abstraction may be transferred from one system to another and then composed with other systems. It does not claim a universal representation or a universal chemistry theorem. The guarantee is conditional on an explicit interface contract, a declared intervention class, and a history state that is sufficient for the contract.

The requirement is stronger than observational fit. Let source and target updates
be F_S:X_S x U_S -> X_S and F_T:X_T x U_T -> X_T, with abstractions
alpha_S:X_S -> Z and alpha_T:X_T -> Z into a common metric space. For related
states (x,y) and translated intervention w(u), an appropriate one-step condition
compares alpha_S(F_S(x,u)) with alpha_T(F_T(y,w(u))). Both sides now have the
same type. The relation and error bound must hold for the interventions and
couplings for which transfer is claimed. Causal abstraction formalizes such
intervention requirements; see [Beckers and Halpern](https://doi.org/10.1609/aaai.v33i01.33012678).

## A conditional transfer-and-closure theorem

Model each physical component \(i\) as a metric transition system

\[
X_i\xrightarrow{(u_i,a_i)}X_i',\qquad y_i=O_i(x_i),
\]

where \(u_i\) is an interface input supplied by a coupling, \(a_i\) is an allowed intervention, and \(x_i\) contains the retained physical state **and any bounded history variable required by the contract**. Let \(\widehat S_i\) be a learned or reduced mechanism. Let \(R_i\subseteq X_i\times\widehat X_i\) relate concrete and abstract states.

Assume, for each component, an \(\varepsilon_i\)-approximate bisimulation relation uniformly over all \((u_i,a_i)\) in a certified set \(U_i\times A_i\): related states have output distance at most \(\varepsilon_i\), and every concrete transition under every allowed input/intervention is matched by an abstract transition (and conversely) whose successor remains related. This is the standard metric-transition-system notion of approximate bisimulation [Girard & Pappas, 2007](https://doi.org/10.1109/TAC.2007.902750); its compositional form under synchronization is given by Julius et al. (2008), including explicit accumulation of approximation parameters [PDF](https://sites.ecse.rpi.edu/~agung/Research/04177873.pdf).

Let the concrete network be closed by a coupling \(C\), and the abstract network by \(\widehat C\). Suppose:

1. **Interface translation:** translated interventions preserve the interface contract: applying \(a_i\) before abstraction and \(w_i(a_i)\) after abstraction yields related interface states.
2. **Coupling closure:** for every related network state and every admissible intervention, \(C\) produces inputs in the certified sets \(U_i\), and \(\widehat C\) produces the translated inputs. For component \(i\), input mismatch obeys
\[
d_U(C_i(y),\widehat C_i(\widehat y))\le \sum_j c_{ij}d_Y(y_j,\widehat y_j)+\rho_i.
\]
3. **History closure:** the interface state is Markov-sufficient for the declared horizon: two related interface states receiving the same future input/intervention traces remain comparable. If projection creates a memory kernel or noise term, that kernel/noise state is included in \(x_i\) or the contract is restricted to inputs for which the omitted term is bounded. Mori–Zwanzig reduction makes the issue explicit: exact reduced equations generally contain a memory function and a force from eliminated degrees of freedom [Hsu & Hsu, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2728514/).
4. **Output composition:** the network output map is \(L_O\)-Lipschitz in component outputs.

Assume additionally that the concrete/abstract component update maps are state/input Lipschitz on the certified domain:
\[
e_{i,t+1}\le a_i e_{i,t}+b_i d_U(u_{i,t},\widehat u_{i,t})+\varepsilon_i.
\]
Let e_i be the component state discrepancy and assume component output discrepancy is at most ell_i e_i. Define \(A_{ij}=a_i\mathbf 1_{i=j}+b_i c_{ij}\ell_j\) and \(q_i=b_i\rho_i+\varepsilon_i\). Then, for any horizon \(T\), any admissible intervention sequence, and any admissible initial related states, the coupled abstract execution has output error bounded by

\[
d(y_t,\widehat y_t)\le \sum_i o_i e_{i,t},
\qquad
e_{t+1}\le A e_t+q,
\]

and therefore

\[
e_T\le A^T e_0+\sum_{k=0}^{T-1}A^kq.
\]

If the spectral radius \(\rho(A)<1\), the steady-state bound is \(e_\infty\le(I-A)^{-1}q\); otherwise the finite-horizon expression is the guarantee. The proof is direct induction using the component state/input Lipschitz inequalities and the coupling gain bound. Approximate bisimulation alone is insufficient when the two executions receive mismatched interface inputs; the \(b_i c_{ij}\) terms are required.

The same statement gives a transfer rule. A mechanism learned in source context \(S\) may be installed in target context \(T\) only after establishing a relation \(R_{S,T}\) with the same uniform input/intervention contract and a transfer error \(\varepsilon_{S,T}\). In a network, replace the component's \(\varepsilon_i\) by its intrinsic approximation plus transfer error. Thus reuse is certified by a measured or proved interface relation, not by similarity of observational trajectories.

## Exact counterexample: observational equivalence fails under coupling

### Correction to an earlier draft

An earlier version stated \(E_{t+1}\le L_CE_t+\sum_i\varepsilon_i\) from coupling Lipschitzness alone. That implication was invalid: it omitted the component sensitivity to mismatched interface inputs. The authoritative theorem above replaces it with the state/input gain matrix \(A\). The earlier scalar recurrence should not be used.

Let each system expose one scalar output \(y\), have hidden state \(h\in\{-1,+1\}\), and accept a scalar coupling intervention \(c\). Define two deterministic mechanisms:

\[
S:\ y=c h,\qquad T:\ y=-c h.
\]

The observational regime is (c=0). For either system and either hidden state, the complete observed trajectory is (y_t=0) at every time. Hence the systems are exactly observationally equivalent in that regime, including all observational histories.

Now connect the component to a coupling that applies (c=1). With (h=+1), (S) outputs (+1) and (T) outputs (-1). No abstraction based only on the observational trace can select the correct transferred sign: the two mechanisms induce the same observational data and opposite intervention responses. A transfer claim based on observational equivalence is therefore false even in a one-step, noiseless, finite-state physical interface.

The counterexample also identifies the missing contract. The allowed intervention set must include (c=1), and the mechanism relation must commute with that intervention. If the coupling can expose or alter the hidden state, the hidden state (or an equivalent sufficient interface variable) must be retained. Restricting the contract to (c=0) makes the theorem vacuous for transfer: it certifies only the regime in which the systems were constructed to be indistinguishable.

## History version and closure failure

A one-state reduction can fail even when the coupling is unchanged. Let the physical state be \((x,m)\), with observed (x=0), hidden memory (m\in\{-1,+1\}), and input sequence (u_t\in\{0,1\}). Set (m_{t+1}=m_t(-1)^{u_t}) and (y_{t+1}=m_t u_t). Under the observational input (u_t=0), every history has (y_t=0), so a memoryless abstraction (x=0) appears exact. After the same observed history, the intervention (u_t=1) gives opposite outputs for (m=+1) and (m=-1). The omitted memory is therefore an intervention-relevant state. This is the discrete analogue of the memory term produced by projection of unresolved physical variables.

## What must be demonstrated for physical reuse

The theorem turns “a mechanism transfers” into four checkable obligations:

- specify the physical variables, history state, coupling ports, and intervention algebra;
- establish a uniform simulation/bisimulation or causal-commutation bound on that declared domain;
- prove that the coupling maps remain inside the certified domain and propagate errors with a known modulus;
- report the composed bound and the experiments or calculations that establish each premise.

Without these obligations, observational accuracy is compatible with arbitrarily bad behavior under the interventions chosen by an inverse-design search. With them, transfer is a conditional engineering theorem whose error budget composes across mechanisms and couplings. It is a route to reusable physical abstractions, not a claim that all matter admits one representation.

## Sources

- A. Girard and G. J. Pappas, “Approximation Metrics for Discrete and Continuous Systems,” *IEEE Transactions on Automatic Control* 52(5), 2007. [DOI](https://doi.org/10.1109/TAC.2007.902750).
- P. Julius et al., “Approximate Equivalence and Approximate Synchronization of Metric Transition Systems,” *Systems & Control Letters* 57, 2008. [Paper](https://sites.ecse.rpi.edu/~agung/Research/04177873.pdf).
- S. Beckers and J. Y. Halpern, “Abstracting Causal Models,” *AAAI* 33, 2019. [DOI/article](https://doi.org/10.1609/aaai.v33i01.33012678).
- A. Geiger et al., “Causal Abstraction for Mechanistic Interpretability,” *JMLR* 26, 2025. [PDF](https://jmlr.org/papers/volume26/23-0058/23-0058.pdf).
- D. Hsu and M. Hsu, “Zwanzig-Mori Projection Operators and EEG Dynamics,” 2009. [Open-access article](https://pmc.ncbi.nlm.nih.gov/articles/PMC2728514/).
