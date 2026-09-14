# Constructive lemmas and their limits

These are elementary derivations that expose assumptions and verification obligations. They are not claims of new mathematical results.

## 1. Closure under interventions

Let $V$ be the span of distinct Pauli observables, orthonormal under $2^{-n}\operatorname{Tr}(P^\dagger Q)$. For a Hermitian Hamiltonian, project $i[H,O_i]$ into $V$ to obtain a real row generator $A$. Cyclicity of trace implies $A^T=-A$. For each discarded row, write $R_i=\sum_{P\notin V}r_{iP}P$ and set $b_i=\sum_P|r_{iP}|$. Because $\|P\|=1$, every density matrix obeys $|\operatorname{Tr}(\rho R_i)|\le b_i$.

For time-dependent controls, variation of constants gives

$$\|z(t)-\widehat z(t)\|_2\le\delta_0+\int_0^t\sqrt{\sum_i b_i(s)^2}\,ds.$$

If every control generator preserves $V$, the residual vanishes for every finite piecewise control sequence. The code checks these generator identities exactly. Learning a small $V$ in broad physical domains is a separate problem.

If the actual Hamiltonian differs by $\Delta H$ with $\|\Delta H\|\le\eta$, each unit-norm observable contributes at most $2\eta$ through $i[\Delta H,O_i]$. Thus the vector residual gains at most $2\sqrt{k}\eta$, multiplied by any declared coefficient-control amplification. This is the bridge from parameter intervals to dynamical uncertainty used in the integrated run.

## 2. Why a damped discarded sector can support bounded memory

Consider a finite linear block system

$$\dot x=Ax+By,\qquad \dot y=Cx+Dy.$$

Assume $\|e^{At}\|\le1$, $\|e^{Dt}\|\le M e^{-\gamma t}$ for $\gamma>0$, $\|B\|\le b$, $\|C\|\le c$, $\|x(t)\|\le X$ and $\|y(0)\|\le Y$. These are explicit hypotheses. Solving the second equation yields

$$y(t)=e^{Dt}y(0)+\int_0^t e^{D(t-s)}Cx(s)\,ds.$$

Hence

$$\|By(t)\|\le bMe^{-\gamma t}Y+\frac{bMcX}{\gamma}(1-e^{-\gamma t}).$$

If $\dot{\widehat x}=A\widehat x$ starts with error $\delta_0$, integrating this forcing gives

$$\|x(t)-\widehat x(t)\|\le\delta_0+
\frac{bMY}{\gamma}(1-e^{-\gamma t})+
\frac{bMcX}{\gamma}\left[t-\frac{1-e^{-\gamma t}}{\gamma}\right].$$

A checkable sufficient condition for $M=1$ is $(D+D^\dagger)/2\preceq-\gamma I$. The discarded sector of a purely Hamiltonian generator is skew-adjoint: rapid oscillation does not establish this strict damping. A useful research result must discover a physical state/observable class with bounded $X,Y$, preserve the assumptions under coupling, and compute these bounds more cheaply than the eliminated dynamics. The expression alone supplies none of that compression.

## 3. Rational many-body energy certificate

Partition a finite Pauli expansion into groups of mutually anticommuting nonidentity strings. For group $Q_g=\sum_{j\in g}c_jP_j$,

$$Q_g^2=\left(\sum_{j\in g}c_j^2\right)I.$$

Choose a rational $u_g\ge0$ with $u_g^2\ge\sum_jc_j^2$. Then $Q_g\succeq-u_gI$, so $H\succeq(c_I-\sum_gu_g)I$. No commutation between different groups is required. An explicit normalized trial state supplies an upper expectation. A model error of norm at most $\eta$ widens both endpoints by $\eta$.

The checker must verify the supplied term multiset, coefficients, anticommutation, square-root inequalities and state expectation. Merely trusting a numerical optimizer or checking coverage is insufficient. The present greedy partition is deliberately a basic constructive baseline; its interval may be loose and no scalable absolute-accuracy theorem is asserted.

## 4. Finite CSS local-generation test

For positive-sign CSS stabilizers, products in each X/Z sector correspond to binary row spans. Solve for combinations cancelling outside region $A$ to obtain global stabilizers supported in $A$. Repeat using only generators supported in buffer $B$. If the global subspace is not contained in the buffer-local subspace, Gaussian elimination returns a concrete missing constraint.

For generators $X_1X_2$ and $X_2X_3$, $A=B=\{1,3\}$ misses $X_1X_3$; including site 2 in $B$ restores that constraint. This is a manufactured diagnostic example with a disconnected region. It is not a counterexample satisfying the geometric hypotheses of a published TQO-2 theorem, and it is not an instance of the 2026 recursive memory construction. Extending the check to that construction requires its actual signed generators, valid region geometry, recursive scale and an induction proof.
