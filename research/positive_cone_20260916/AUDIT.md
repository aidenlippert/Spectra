# Audit of the positive-cone amplitude certificate

After the down-spin particle-hole transformation, the amplitude map is

\[
 L(C)=KC+CK^T-8\sum_i n_iCn_i+32C.
\]

The exact particle-hole replay checks all signed transitions on both local and
full balanced configurations. Five focused tests pass. The submitted rational
candidate has `C\succ0` and a positive residual, with reported width
`4.77409316436\times10^{-7}t`.

## Why the certificate is valid for this Hubbard map

For a generic linear map, `C\succ0` and `L(C)-\ell C\succeq0` do not imply a
spectral lower bound. A precise counterexample is the self-adjoint map on
`2\times2` matrices `L(X)=\operatorname{Tr}(X)I-X`. With `C=I` and `\ell=1`,
the residual is zero, while every traceless matrix has eigenvalue `-1`.

The Hubbard map has extra structure. Its Hilbert-Schmidt Rayleigh minimum may
be taken over Hermitian `X`. In an eigenbasis of `X`,

\[
 \operatorname{Tr}(KX^2)=\operatorname{Tr}(K|X|^2),
 \qquad
 \operatorname{Tr}(XN_iXN_i)=\sum_{ab}\lambda_a\lambda_b|(N_i)_{ab}|^2
 \le\operatorname{Tr}(|X|N_i|X|N_i).
\]

Since the latter term has coefficient `-U`, `|X|` cannot increase the
Rayleigh quotient. A lowest eigenmatrix `W` can therefore be chosen PSD. If
`D=L(C)-\ell C\succeq0` and `C\succ0`, then

\[
 0\le\operatorname{Tr}(WD)=(E_0-\ell)\operatorname{Tr}(WC),
\]

and `Tr(WC)>0`, proving `E_0\ge\ell`. Strict positivity is needed for this
overlap, not because an inverse is taken. No irreducibility or CP-semigroup
assumption is needed for this real Hermitian map.

## Scope

At fixed total particle number `N=8`, the balanced `M_S=0` sector equals the
unrestricted spin minimum: every even-`N` SU(2) irrep has integer spin and
contains `M_S=0`. Other particle numbers are outside the result. The replay
supports a fixed-sector rational lower certificate; it does not establish a
universal many-body solver or a physical-model error bound. A complex or
non-Hermitian extension would require a separate cone/Perron argument.
