# Frustration-free parent plus relative residual

Let (P\succeq0) be an exactly represented parent operator, with normalized
trial state (|\phi\rangle\) in its common kernel: (P|\phi\rangle=0). Set

\[
 E_\phi=\langle\phi|H|\phi\rangle,\qquad R=H-E_\phi I-P.
\]

If there are exact constants (0\le\eta<1) and \(\varepsilon\ge0) such that

\[
 R\succeq-\eta P-\varepsilon I,
\]

then (H\succeq(E_\phi-\varepsilon)I+(1-\eta)P\succeq
(E_\phi-\varepsilon)I). The variational principle gives

\[
 E_\phi-\varepsilon\le E_0(H)\le E_\phi.
\]

The theorem does not need a spectral-gap oracle. It does require: normalized
trial state; a common kernel (or a separately certified kernel expectation
bound if the kernel is multidimensional); exact operator ordering and CAR
replay; and a relative residual inequality valid on the full Hilbert space,
including the parent kernel. If (P=\sum_iP_i) is local and the inequality is
proved by local terms, the total \(\varepsilon\) generally scales as an
extensive sum of local errors. Physical locality alone does not make it
constant; a boundary-only theorem needs an independent boundary estimate.

## Molecular premise test

`research/wave2_20260913/theory/parent_test.py` uses the active-space ladder H4
fixture, fixed (N=4), dimension 70. It chooses the lowest-diagonal
occupation determinant as a simple HF parent and sets (P=I-|HF\rangle\langle
HF|), rather than fitting (P) to the residual. Exact rational LDL arithmetic
certifies the relative inequality at \(\eta=1,\varepsilon=4\), yielding

* (E_\phi=-3.525719766889),
Exact shift searches give widths (epsilon) 0.262389 (eta=0), 0.549556
(eta=1/2), and 0.141281 (eta=1). Even the best tested case is far above the
0.0016-Hartree target. The parent is a global diagnostic projector, not a
local sum of HF occupation penalties. Width is epsilon, not 2 epsilon.

This is a valid but useless certificate, and is a falsifiable rejection of the
claim that a cheap parent automatically gives chemical accuracy. The next
test must seek a physically motivated local parent and measure whether its
smallest exact \(\varepsilon\) falls near (10^{-3}) Hartree while preserving
the common-kernel condition. A dense 70-state matrix is only a diagnostic;
its cost is exponential in system size.

The checker is exact for the rational fixture and independently replays the
matrix. It does not use the saved FCI value or a quantum gap estimate.
