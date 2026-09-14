# Corrected original H4 selected-subspace result

The prior Wave2 run accidentally used the Boys fixture. This rerun uses the
original active-space fixture, SHA256
`8e64505deb4cc06a87f75e669b0f3073d15ebca72daa6d086e06eeea0b0f5120`, and loads
the actual Wave1 supports from `step_01_upper.json` (16 states) and
`step_02_upper.json` (20 states). The Boys result remains preserved separately
and is explicitly marked as a different fixture.

Exact LDL floors certify both projected and complement blocks. The 16-state
case gives `a=-18333059109/5000000000`, `mu=-35752753711/10000000000`, exact
coupling Frobenius square about `0.2506065`, and valid block endpoint
`-16494513953/4000000000` (about -4.12363). The 20-state support gives
`a=-916749989/250000000`, `mu=-35752753711/10000000000`, and exactly zero
coupling: the selected support is invariant under this Hamiltonian replay.
Thus the endpoint is exactly `-916749989/250000000` (about -3.666999956),
matching the saved H4 reference interval to displayed precision. This is an
actual successful finite H4 selected-subspace certificate, though identifying
an invariant 20-state support and proving its complement remains an
exponential-size diagnostic (2,500 complement entries).

The LDL checker rejects negative pivots and zero-pivot rows with nonzero
coupling; square-root ceilings are increased until their squared rational value
dominates the exact radicand. No H6 dense computation was performed.
