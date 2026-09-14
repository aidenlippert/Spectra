# H8 full-support FCI validation control

Generated from the frozen canonical H8 rational fixture with the existing validation-only FCI script, requesting `top=5000` (larger than the 4,900 alpha/beta determinant pairs). The resulting integer witness has 2,468 nonzero amplitudes after 1e10 rounding. This is a validation upper only and is prohibited as SOS discovery input.

The recorded numerical FCI electronic energy is -9.255304429088598 Ha; the exact CAR replay of the rounded witness is recorded in `h8/receipt.json` and `h8/upper.json`. Wall time was 13.58 s.
