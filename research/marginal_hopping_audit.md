# Bounded audit for the hopping deformation

Consider six spinless modes, (N=3), grouped into two triples, with

\[
H(t)=\sum_{T}\sum_{i<j\in T}n_i n_j-t\sum_{(u,v)\in E}(a_u^\dagger a_v+a_v^\dagger a_u),
\]

where (E) is the declared matching of inter-triple edges. This is a useful smallest perturbation of the diagonal witness, but it has a decisive completeness trap: in the (N=3), (M=6) sector, unrestricted cubic creation/annihilation words can span the entire 20-dimensional physical Hilbert space. A full degree-3 SOS search can therefore reproduce an arbitrary positive operator on the sector and make the result a tautology. It tests the solver, not sparse representation.

## Fair baselines

Run all of these at each (t), with identical exact or interval tolerances:

1. Full DQG ((D,Q,G)) relaxation.
2. DQG plus the original (2r=2) local cubic anticommutators.
3. A fixed sparse dictionary containing only local pair words, number words, and the selected inter-triple hopping words.
4. An adaptive dictionary with a strict budget (k) on cubic words and support size, selected by violation magnitude.
5. Full degree-3 SOS, labeled as a completeness upper baseline rather than a scalable result.
6. Exact diagonalization of the (\binom63=20) dimensional sector.

The meaningful curve is certified gap versus dictionary budget (k), support radius, and maximum ladder degree. Compare adaptive selection with random words and a deterministic largest-violation oracle. Report wall time only after separating SDP solve, candidate generation, symbolic reduction, and proof verification.

## What counts as sparse

Fix the orbital basis and edge set before training. Count a cubic word by the number of distinct orbital labels and by graph support on the triple-interaction/hopping graph. A dense orbital rotation may preserve an SOS identity while making every coefficient nonzero, so do not count a rotated dense coefficient tensor as sparse. Include coefficient bit length and number of nonzero monomials in the budget. Re-run after a declared symmetry-preserving orbital rotation as a covariance control.

Use three regimes: (t=0), where the local witness is known; small (t), where perturbative continuity is expected; and a scan through avoided crossings, where local occupancy cuts should fail sharply. A useful falsification is an extensive or abrupt budget increase at the first nonzero hopping scale, despite the exact problem remaining tiny.

## Exact certificate export

Represent every word in a canonical normal-ordered CAR basis with rational coefficients. Given a proposed Gram matrix (Q\succeq0), export a rational (Q_R), an exact normal-ordered residual

\[
R=H-bI-w^\dagger Q_Rw-Z,
\]

and a sector-specific bound. Verify Hermiticity coefficient by coefficient, verify (P_3 ZP_3=0) as an exact matrix identity (or retain an explicit number-ideal multiplier), and verify (Q_R\succeq0) by an exact (LDL^\dagger) factorization with nonnegative rational pivots or by an interval Cholesky certificate. For floating-point (Q), rational reconstruction must be followed by re-expansion; proximity to the numerical SDP solution is not evidence of validity.

If the remaining sector residual is nonzero, compute a rigorous operator-norm bound η on (P_3RP_3) and report (b_{safe}=b-η). The operator Frobenius norm is a valid, usually loose, spectral bound. A coefficient-array norm is not valid without the Gram matrix of the word basis. Include all SDP and rationalization tolerances in η.

## Smallest meaningful dictionary

For this six-mode test, use the original two local cubic annihilators and their adjoints only as a seed, then add at most (k=1,2,4,8) *mixed* cubic words with one inter-triple edge, together with the linear number/hopping words needed to match (H(t)). Impose particle-number charge homogeneity, Hermitian conjugate closure, and the fixed support budget. This keeps the test from silently becoming full-sector tomography while allowing the dictionary to discover dressed certificates.

The primary output is the table

\[
(t,k,d,\#\mathrm{words},\#\mathrm{nonzeros},\mathrm{bits},
E_{lower}^{safe},E_{upper},E_{upper}-E_{lower}^{safe}).
\]

Success means a small-budget certificate closes the gap over a nontrivial interval of (t), survives exact CAR re-expansion and rational PSD verification, and beats random-word and fixed-local baselines. A single six-mode success cannot support a scaling claim; the next control should replicate the same budget discipline on disjoint copies while keeping hopping degree bounded.

