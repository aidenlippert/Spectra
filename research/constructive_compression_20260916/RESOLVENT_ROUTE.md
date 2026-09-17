# Block-local resolvent route

Let (A=Q(H-b)Q) be an eliminated block and (B=QHP). Choose a block-local
positive preconditioner (D) with an independently checkable inverse action,
and write (A=D+V). If

\[
 \rho=\|D^{-1/2}VD^{-1/2}\|<1,
\]

then the noncommutative Neumann series is norm convergent:

\[
 A^{-1}=D^{-1/2}\sum_{k\ge0}(-W)^kD^{-1/2},
 \qquad W=D^{-1/2}VD^{-1/2}.
\]

The truncation after (r) terms has the exact bound

\[
 \left\|A^{-1}-D^{-1/2}\sum_{k=0}^{r-1}(-W)^kD^{-1/2}\right\|
 \le {\rho^r\over(1-\rho)}\|D^{-1}\|.
\]

For a response, apply this operator only to (B). The omitted Schur response
is bounded by

\[
 \|B^*(A^{-1}-A_r^{-1})B\|
 \le {\rho^r\over1-\rho}\|D^{-1/2}B\|^2.
\]

This is a genuine compact construction only when (D^{-1}), (V), and (B)
are represented by local operator rules and their norm bounds are certified
without constructing the (Q)-space matrix. A finite-depth word expansion by
itself does not make (D^{-1}) local.

## Local norm criterion

If (D\succeq dI), a safe input-only estimate is

\[
 \rho\le {\|V\|\over d}
 \le {1\over d}\sum_C\|V_C\|.
\]

The sum is required even for disjoint supports; disjoint tensor factors do not
turn a sum's norm into the maximum local norm. Sharper bounds require a
separate commuting spectral or relative-form certificate. A block-local
relative SOS inequality (V^*D^{-1}V\preceq\rho^2D) is preferable and is
directly composable if its boundary terms are retained.

## Inconclusive elementary ladder estimate

For the half-filled (2\times4) repulsive Hubbard ladder with (U/t=8), a zero-shift
first naive choice is (D=Usum_i n_{i\uparrow}n_{i\downarrow}) on a sector
where at least one doublon supplies coercivity. The hopping remainder contains
all ten coupled edges. The elementary CAR word triangle bound gives an
extensive estimate of order (2|E|t=20t), so the scalar criterion gives

\[
 \rho\lesssim20/8=2.5>1.
\]

Therefore this un-dressed criterion cannot certify convergence on the target
ladder. This is a failure of the chosen preconditioner/bound, not a proof that
the resolvent cannot be compressed. A successful dressing must absorb enough
hopping into (D), certify a relative form bound rather than a raw triangle
bound, or use a smaller eliminated sector with a stronger excitation count.

The implementation seam is clear: a local block builder emits (D,V,B), an
exact checker verifies local PSD and relative-form certificates, and a bounded
word engine emits (A_r^{-1}B) and the residual allowance. The retained Schur
operator still needs its own positivity certificate. This supplies a standard sufficient criterion. The bound above is inconclusive;
no new ladder convergence certificate or general compression theorem was obtained.
