# A bounded dissipative closure lemma

Use a Hilbert–Schmidt-orthonormal Pauli basis, removing identity. The coefficient generator is `L(t)=K(t)-D(t)`, with real skew `K(t)` and diagonal damping `D(t)>=gamma I` on all retained and discarded coordinates.

For an orthogonal retained/discarded split write

\[
L=\begin{pmatrix}A&B\\C&E\end{pmatrix},
\qquad E=K_{QQ}-D_Q,\quad D_Q\succeq\gamma I.
\]

The homogeneous discarded equation gives

\[
\frac{d}{dt}\|y\|_2^2
=2y^T(K_{QQ}-D_Q)y\le-2\gamma\|y\|_2^2.
\]

Thus its propagator obeys `||Phi_Q(t,s)||_2<=exp(-gamma(t-s))`. This remains true for bounded measurable time-varying skew Hamiltonian controls and a uniform damping lower bound. The memory kernel satisfies

\[
\|B(t)\Phi_Q(t,s)C(s)\|_2
\le b c e^{-\gamma(t-s)},
\]

where b and c uniformly bound the coupling norms over the admitted controls. Zero initial discarded coordinates eliminate the separate inhomogeneous transient term.

If the retained propagator decays at rate alpha>0 and `||z(s)||_2<=M`, variation of constants yields

\[
\|z(t)-\widehat z(t)\|_2
\le bcM\int_0^t e^{-\alpha(t-s)}\frac{1-e^{-\gamma s}}\gamma\,ds
\le\frac{bcM}{\gamma\alpha}(1-e^{-\alpha t}).
\]

Under the stronger full-generator contraction and retained contraction at rate gamma, `M(s)<=||z(0)||exp(-gamma s)` gives the sharper envelope

\[
\|z(t)-\widehat z(t)\|_2
\le bc\|z(0)\|_2\frac{t^2}{2}e^{-\gamma t}.
\]

The finite block checker in [v3_memory.py](../../experiments/v3_memory.py) includes damping in **both** sectors. Damping only the discarded block does not imply the uniform error estimate used there.

These are mathematical coefficient-norm bounds. A physical observable claim can hide exponential dimension factors when converting from the Hilbert–Schmidt norm. The [physical-norm audit](physical_norm_obstruction.md) identifies that failure, and the [local Pauli l1 theorem](local_pauli_theorem.md) supplies a separate dimension-independent physical bound under stronger, explicit locality and damping conditions.

No self-correcting quantum memory or universal compactness follows. Finding small coupling blocks can still be expensive, and a physically useful source must independently justify its noise model and damping rates.
