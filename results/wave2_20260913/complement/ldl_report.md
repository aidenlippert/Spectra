# H4 selected-subspace exact-LDL Schur control

Using the actual Wave1 H4 supports (16 determinants from `step_01_upper.json`
and 20 from `step_03_upper.json`), I reconstructed the 70-dimensional CAR
Hamiltonian exactly. For each selected block `A` and complement `C`, a floating
eigenvalue only proposes a 1e-10 rational floor; exact fraction LDL inertia
replays `M-lI >= 0`. The coupling is the exact Frobenius square of `B`.

The general block inequality used is
`lambda_min(H) >= (a+mu-sqrt((a-mu)^2+4||B||^2))/2`, with the square root
rounded upward, so it remains valid even without `mu>a`. Both selected supports
actually passed the sharper gate `mu>a`.

For 16 states: `a=-36523947161/1e10`, `mu=-35752753711/1e10`,
`||B||_F^2=173333636329311272497151/250000000000000000000000`, endpoint
`-22236972579/5000000000` (about -4.44739). For 20 states: `a=-2282749833/625000000`,
same `mu`, `||B||_F^2=97728619254766533128611/125000000000000000000000`, endpoint
`-89977779211/20000000000` (about -4.49889). The numerical H4 ground oracle is
about -3.666999956, so both are valid but loose lower bounds. Complement costs
are 2,916 and 2,500 entries. Diagonal-only values remain diagnostics and are
never used as bounds.

This improves the Wave1 single-determinant analysis by producing accepted
whole-block bounds, but accuracy remains poor because Frobenius coupling and
LDL lower floors leave a large Schur penalty. No scalability claim follows;
H4 full complement is a diagnostic. The reusable obstruction is quantitative:
an accurate selected state does not suffice unless the complement lower bound
is close to the selected energy and the coupling norm is small.
