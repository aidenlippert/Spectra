# Local tensor proof for the recognized chain

The implemented family is

\[
H=a\sum_{i=1}^{M-1}[-(a_i^\dagger a_{i+1}+\mathrm{h.c.})+n_i n_{i+1}]+\mu N+c,
\qquad a>0.
\]

Its hidden supersymmetry is known: [Yang–Fendley](https://arxiv.org/abs/cond-mat/0404682), [Hagendorf–Liénardy](https://arxiv.org/html/1612.02951). Our verifier uses the spin-reversed local map

\[
q|0\rangle=|11\rangle,\quad q|1\rangle=0,\quad
Q_M=\sum_{j=1}^{M}(-1)^{j-1}q_j.
\]

It checks coassociativity, \((q\otimes I-I\otimes q)q=0\), and the local Hamiltonian identity. These imply

\[
Q_{M+1}Q_M=0,\qquad
H-[c+\mu N-a(M-N)]I
=a(Q_M^\dagger Q_M+Q_{M-1}Q_{M-1}^\dagger)\succeq0.
\]

For even half filling, let

\[
C=(|01\rangle+|10\rangle)^{\otimes M/2},\qquad
\omega=|01\rangle^{\otimes M/2}.
\]

The checked two-site identities give \(Q_MC=0\), \(Q_{M-1}^\dagger\omega=0\), and \(\langle\omega,C\rangle=1\). Hence C is outside \(\mathrm{im}\,Q_{M-1}\). Its orthogonal projection onto that image's complement is a nonzero simultaneous kernel vector. The particle grading keeps this vector in the required sector. Thus

\[
E_0=c+\mu M/2-aM/2.
\]

The proof checks constant-size tensors and the input coefficients. It does not construct the projected wavefunction. This reproduces a known solvable class; its auxiliary maps change chain length, so it is a separate certificate format from the main fixed-M CAR/SOS engine.
