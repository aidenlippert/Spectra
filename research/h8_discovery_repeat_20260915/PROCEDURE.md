# Repeat the successful H8 discovery before claiming transfer

The existing 0.767448135463 mHa interval is an achieved, preserved result. This
experiment asks whether a new optimization, under a procedure frozen before
execution, produces another accepted interval at most 1.6 mHa wide. Its trial
count is one. It is not a statistical reliability estimate.

The executable caller is [run.py](/Users/aidenlippert/Documents/Spectra/research/h8_discovery_repeat_20260915/run.py).
The exact numerical commands, input/source hashes, stopping rules, and environment
are fixed in [protocol.json](/Users/aidenlippert/Documents/Spectra/results/h8_discovery_repeat_20260915/protocol.json).
The supplied optimizer and accepting modules are unchanged.

## Meaning of a fresh run

The output directory and the Numba cache start empty. This pass recomputes the
actual degree-six MPS moments, their coefficient-functional pullback, the sparse
spin representation, and the sparse normal matrix. It then runs the supplied
130-second bound-maximization stage and 180-second refinement stage, with their
existing numerical stopping rules. Export, pruning, and exact acceptance follow.
Hard process limits also cover setup and export. An unsuccessful step is retained;
there is no replacement trajectory hidden behind a successful result.

The old compact checkpoint, physical coefficient maps/frame, rational MPS, and
nonsinglet proof are inherited. This is fresh optimization conditional on those
inputs. It does not remove their construction costs or demonstrate discovery
from a Hamiltonian alone. The previously successful spin-adapted checkpoint and
the stronger full-cubic factors are not optimizer inputs.

## What accepts the new energy

Let H_s be the exact SU(2) average of the rational input H and let delta_H be
the checked coefficient-norm bound on H-H_s. For the current H8 input,

    delta_H = 59/250000000000 Ha.

The rational exported factors are actual operators F_j. Their positive squares
are expanded with the original fermionic algebra. After subtracting b and the
declared number/spin ideal terms from the auxiliary Hamiltonian, the checker
computes the actual residual R. Exact spin averaging and the coefficient L1
allowance give the singlet bound

    L_0 = b - ||T(R)||_(coefficient,1).

This follows because the ideal terms vanish under singlet compression, averaged
squares remain positive, and every normal-ordered fermionic monomial has operator
norm at most one. It does not depend on the floating solver declaring convergence.

The separate M_S=1 certificate supplies L_1. With even electron number, each
integer-spin multiplet with S>=1 has an M_S=1 member. Therefore

    L = min(L_0, L_1) - delta_H.

The original checker explicitly rejects odd electron counts on this path.
The original rational MPS contraction supplies U anew. Acceptance requires

    L <= U,                  U-L <= 1/625 Ha = 1.6 mHa.

This is total electronic-energy width for the declared model. It is not a
per-electron tolerance, a physical-model error bar, or a certified family optimum.

## Cost boundary

The useful future measure is time until an accepted interval reaches the target:

    T_epsilon(H) = first elapsed time at which the complete checker accepts
                   an interval for H with width at most epsilon.

Report the starting inputs with that number. In this pass, process receipts cover
all seven new steps, and the enclosing execution receipt also charges protocol
freezing and repeated input checks. Historical dependencies and post-run report
generation are explicitly outside that elapsed time. A miss under the budget is
a recorded miss; it is not assigned a successful runtime or interpreted as a
proof that the mathematical family cannot succeed.

For a future fresh-problem run, the measured boundary must include integral and
model construction, state discovery, physical-map preparation, compact seed
construction (or its declared replacement), nonsinglet discovery, unsuccessful
searches, export, and replay. Shared setup must be reported separately if it is
amortized across related geometries. Historical stage times have overlapping
dependencies and different timing boundaries; adding selected successful-stage
numbers does not produce a complete cold runtime.

## Concrete transfer boundary found in the current code

The supplied scripts are specialized to H8. The raising-map constructor uses
physical group indices 42 through 57. The residual predictor uses 16 spin orbitals
and 8 electrons. Paths and the upper endpoint are fixed to this fixture.
An input-independent procedure must derive these from explicit model metadata.

A geometry change that breaks the current spatial parity must regenerate the
physical dictionaries and coefficient equations, retaining newly allowed
couplings. Replacing only the Hamiltonian vector inside H8's reduced equations
is not a valid transfer test. A different molecule also requires a new upper and
a new complete lower-sector argument. No current H8 nonsinglet bound transfers
automatically to a new Hamiltonian.

Frozen adaptation rules may change dimensions with the model. The point is to
freeze how charge/spin/parity groups are formed, which operator spaces are kept,
how the initial state and seed are produced, and the optimization/export rules
before inspecting held-out outcomes. Constant matrix dimensions are not required.

The next transfer design should contain a changed H8 geometry that removes a
used spatial symmetry, a different even-electron molecular model, and a modest
size or basis increase. Each needs the complete interval and all dependencies.
These cases have not been executed by this repeat. Accuracy-matched numerical
comparisons and guarantee-matched certificate comparisons are separate controls.

## Diagnostics needed when a later case misses

Keep the trial upper, optimization error, exported residual/spin allowances,
and possible limitation of the certificate family distinct. A verified dual
ceiling compatible with the actual residual rule is needed to establish a
limitation of that family; a nonconverged dual is insufficient.

To separate representation efficiency from mathematical freedom, compare
equivalent attainable certificate cones with equivalent initial physical points,
stopping rules, and hardware. Compare a restricted cone separately. This repeat
does not assign the improvement uniquely to sparsity, spin organization, or
the added multiplicities.

For noninteracting fragments at fixed individual electron counts, the exact
energies add. Compare a direct combined certificate with the sum of separately
certified fragment intervals. The diagnostic extra width is

    W_direct(A+B) - [W(A)+W(B)].

Fixing only the total electron count can permit charge redistribution, so it
does not justify the same additivity test. Experimental comparisons additionally
require specified states, geometry, physical-model errors, and the uncertainty
of the actual energy difference rather than each absolute energy alone.
