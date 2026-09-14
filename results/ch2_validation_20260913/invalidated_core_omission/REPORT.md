# CH2 validation pilot: implementation obstruction

## Scope fixed before calculation

The intended observable is the fixed-geometry electronic singlet–triplet
gap

\[
\Delta_{ST}=E_T-E_S.
\]

A positive value means the singlet is lower.  The intended first calculation
requires one explicit geometry, basis, frozen-core convention, and active
space, followed by separate total-spin `S=0` and `S=1` certificates.  A
calculation restricted only to `M_S=0` would not satisfy this requirement.

## Result

Using the approved molecule runtime (`.venv-molecule`, PySCF 2.14.0), I ran a
fixed-geometry STO-3G active-space control. The geometry is C at the origin,
H at `(0,0,1.117)` Å and H at `(0,1.047,-0.389)` Å. One RHF core spatial
orbital is frozen; four active spatial orbitals and four active electrons are
retained. The M_S=0 active determinant space has dimension 36.

The lowest separately identified roots are

```
S=0: E = -15.096934954026468 Ha, <S^2> = 0
S=1: E = -15.089349106876613 Ha, <S^2> = 2
Delta_ST = E_T - E_S = +0.007585847149855 Ha
                         = +4.762 kcal/mol
```

The sign is therefore singlet lower for this finite-basis active-space
Hamiltonian. The result is a numerical proposal/control, not a certified
chemical prediction: the current run does not yet export a rational Hamiltonian
matrix and independently replay its eigenvalue intervals. It consequently
does not satisfy the project's exact-certificate requirement or establish
agreement with the experimental CH2 gap.

## State-labeling control

`spin_sector_control.py` constructs the two-electron, two-spatial-orbital
`M_S=0` determinant block and diagonalizes its exact `S^2` operator. It
recovers one `S=0` eigenvector and three `S=1` eigenvectors, then checks an
exact synthetic Hamiltonian with `E_S=0`, `E_T=2`, hence
`Delta_ST=+2 Ha`. This control is deliberately not presented as a CH2 model;
its purpose is to prevent an invalid `M_S=0` state label when the chemistry
runtime becomes available.

## Accounting and next action

The earlier 4-electron and no-core proposals were invalidated and preserved.
STO-3G CH2 has seven spatial orbitals. The corrected calculation freezes one
RHF spatial core and treats the remaining six orbitals with six active
electrons (3 alpha, 3 beta), giving 400 M_S=0 determinants. Separate
spin-adapted PySCF solvers were used: `direct_spin0` for S=0 and an M_S=1
(4 alpha, 2 beta) calculation for S=1. It gives

```
S=0: -19.66114749725123 Ha (active electronic part)
S=1: -19.74038969818306 Ha (active electronic part)
Delta_ST = -0.079242200931827 Ha = -49.72 kcal/mol
```

The triplet is lower in this fixed STO-3G model. The exact rational export and
independent lower/upper replay were attempted as the next step but remain
incomplete; this remains a numerical state-resolved control,
not a certified interval or an experimental comparison.
Basis, geometry, frozen-core, zero-point, and reference uncertainties remain
separate from solver error.
