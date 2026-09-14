# Local norm obstruction and a conditional dissipative route

## The hidden exponential in Hilbert–Schmidt coordinates

For (n) qubits, write a traceless observable as

\[
 O=\sum_{P\ne I} c_P P,
 \qquad c_P=2^{-n}\operatorname{Tr}(PO).
\]

The Euclidean coefficient norm is tied to Hilbert–Schmidt norm:

\[
 \sum_P c_P^2=2^{-n}\operatorname{Tr}(O^2).
\]

For a state, the normalized Pauli expectation vector has squared norm

\[
 \sum_{P\ne I}\langle P\rangle_\rho^2
 =2^n\operatorname{Tr}(\rho^2)-1,
\]

which can be (2^n-1) for a pure state. Thus an (\ell_2) error bound on all Pauli coordinates can hide a factor (2^{n/2}) when converted to a general expectation or trace-distance claim. The Hilbert–Schmidt-to-operator-norm conversion has the same dimension dependence. A compact reduced-coordinate theorem based only on Euclidean error therefore does not establish a scalable local-observable guarantee.

Even an (\ell_2) bound can be misleading for a local observable after Heisenberg evolution: converting coefficient error to a worst-case expectation error uses

\[
 |\operatorname{Tr}(\rho\,\Delta O)|
 \leq \sum_P|\Delta c_P|,
\]

and (\ell_1\leq\sqrt{|\operatorname{supp}\Delta O|}\ell_2). If support becomes extensive, this restores an exponential factor. The relevant target is a locality-aware coefficient norm, not a global normalized-vector norm.

## Conditional local \(\ell_1\) route

Let (H=\sum_X h_XP_X) be a Pauli Hamiltonian with bounded interaction support. In the Heisenberg picture, a Pauli column (Q) satisfies

\[
 i[h_XP_X,Q]=0
 \quad\text{or}\quad
 i[h_XP_X,Q]=\pm 2h_XR_{X,Q},
\]

where (R_{X,Q}) is one Pauli string. Let (w(Q)) be Pauli weight and define

\[
 \|O\|_{1,w}=\sum_{Q\ne I}|c_Q|w(Q).
\]

Under independent local depolarization at rate \(\gamma\), a nonidentity Pauli of weight (w(Q)) receives damping (-\gamma w(Q)c_Q). If each qubit belongs to at most (D) Hamiltonian terms and \(|h_X|\leq J), the absolute Hamiltonian column growth obeys

\[
 2\sum_{X:X\cap\operatorname{supp}(Q)\ne\varnothing}|h_X|
 \leq 2DJ\,w(Q).
\]

For the unweighted coefficient norm \(\|O\|_1=\sum|c_Q|), the same bound gives a sufficient uniform contraction condition

\[
 \gamma>2DJ.
\]

More generally, for a weighted norm (\sum |c_Q|\alpha^{w(Q)}), one must include the support change in each commutator:

\[
 2\sum_X|h_X|\alpha^{w(R_{X,Q})-w(Q)}
 <\gamma w(Q)
\]

uniformly over allowed (X,Q). A local interaction of bounded size (k) changes weight by at most (k), so a conservative threshold can be stated using the worst allowed factor (\alpha^{k}), but that threshold may become unnecessarily strong. The unweighted (\ell_1) estimate avoids that artificial factor and is the first bound to test.

If the contraction holds, an initially local observable with coefficient (\ell_1) norm (C) has an expectation-error bound of the form

\[
 \sup_\rho|\operatorname{Tr}[\rho(O_t-\widehat O_t)]|
 \leq C e^{-\kappa t}\,\text{(initial error)}+\text{forcing/error terms},
 \qquad \kappa=\gamma-2DJ,
\]

with constants independent of total qubit number, provided the forcing norm is measured in the same local (\ell_1) norm. Since every Pauli expectation is at most one, the final conversion from coefficient (\ell_1) error to a local expectation error is dimension independent.

## Conditions that must remain explicit

This route is conditional on bounded interaction degree, bounded local coefficients, a specified Markovian local-depolarization generator, and an error source controlled in the same norm. It does not establish a generic physical law of dissipative contraction. A nonlocal Hamiltonian, long-range interactions with unbounded incident strength, correlated noise, or an error model specified only in Hilbert-Schmidt norm can reintroduce dimension dependence.

The threshold also cannot be advertised as a free scientific breakthrough. Taking \(\gamma\) arbitrarily large can erase the clock, signal, and controllable dynamics along with the error. Any useful theorem must retain a declared observable task, nonzero signal, preparation and measurement costs, and a comparison against an equally damped baseline. Otherwise the result is a stable trivial fixed point rather than a constructive representation of matter.

The appropriate v3 claim is therefore: under the stated local bounded-degree and damping assumptions, a local weighted-\ell_1 representation can have system-size-independent expectation-error bounds. It is not a global state-norm theorem, not a generic closure result, and not evidence for broad physical controllability.
