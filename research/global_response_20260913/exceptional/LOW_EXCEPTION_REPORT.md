# Low-occupation exception with its full complement

The meaningful partition is now explicit.  Let (N_h) count electrons in the
last two spatial orbitals and define

\[
P=\{N_h\le 1\},\qquad R=I-P=\{N_h\ge2\}.
\]

The R certificate in this directory was independently proposed and replayed
with the joint density-square/filling construction for the whole complement;
it is not a certificate for a selected subset of R.  The construction uses a
4-orbital one-particle projector and a 6-by-6 (H6) or 8-by-8 (H8) PSD check.
No determinant basis is built.  Determinant dimensions below are combinatorial
oracle counts only: H6 has |R|=672 and |P|=252; H8 has |R|=9207 and
|P|=3663.

For H6, replay proves

\[
R(H-b)R\succeq -6.29381987977925\,R,
\]

which is 39.24 mHa above the frozen target b≈−6.333059 Ha and therefore gives
a positive complement denominator.  A new local CAR decomposition uses source
labels `[0,1,2,4,8]` and every output label with at least two occupied local
modes.  The corrected exact replay finds a nonzero P–R coupling.  Its squared
norm majorant is \(\nu\le7.04204811188\,\mathrm{Ha}^2\).  With the H6 target
\(b=U-0.001\), the complement gap is 0.04023874645 Ha, so the plain scalar
Schur penalty is approximately 175.007 Ha.  This does not close the global
bound; a response-preserving coupling treatment is required.

For H8 the same full-complement lower bound is approximately −9.261457 Ha,
below its frozen target, so this partition does not even provide a positive
denominator there.

Thus H6 supplies a useful, genuinely partition-complete intermediate result:
the compact low-occupation exception has a certified entire complement, while
the remaining obstruction is quantitatively isolated in the coupling Schur
term.  H8 is an honest failure of the same orbital pattern.  Neither result is
a full ground-energy certificate.

Replay receipts are in `results/global_response_20260913/exceptional/`; the
 exact joint checker reports zero many-body states and matrix entries constructed.
The earlier zero-coupling receipt is invalidated under
`invalidated_zero_coupling/`.
