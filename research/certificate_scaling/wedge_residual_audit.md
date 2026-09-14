# Wedge residual bound audit

The independent controls construct occupation-basis matrices for balanced
normal-ordered residual bodies of orders 2, 4, and 6 and compare every fixed
particle-number minimum eigenvalue against `residual_bounds`.  The audit also
checks rejection of unbalanced/non-Hermitian input, a positive-projector
control, and the vanishing of a body-2 term on the one-particle sector.

The bound is a Gershgorin lower bound on each exterior-power coefficient
matrix, multiplied by `C(N,k)` for the fixed-number lift, then combined with
the coefficient L1 bound.  It is therefore a residual improvement for a fixed
N sector; it does not replace exact SOS replay or establish a smaller
certificate dictionary.
