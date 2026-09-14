# Transfer test: connected six-mode graph

The coefficient-space certificate engine was tested on the asymmetric six-mode
Hamiltonian already used in the project, with an additional hopping edge
`(1,2)` of weight `delta`. At `delta > 0` the graph is connected, so the
matched-flavor component decomposition is gone. The coupling is `t=1/5` and
there are three particles. Discovery uses the mixed cubic CAR dictionary and
the particle-number ideal; it constructs coefficient maps and PSD Gram blocks,
never a Fock-sector matrix.

The independently computed upper witness is a full 20-dimensional numerical
Rayleigh vector, rounded to integer amplitudes and evaluated by exact signed
fermion action. It is an upper bound only; it is not used to certify a lower
bound.

| delta | lower bound | exact replay upper | interval width |
|---:|---:|---:|---:|
| 0 | 0.540868577421760 | 0.540877180? | 8.60371e-6 |
| 0.01 | 0.540850847887070 | 0.540866357? | 1.55129e-5 |
| 0.1 | 0.539783608921510 | 0.539798969? | 1.53702e-5 |

The question tested here has a clean answer: tightness largely survives a
connected perturbation, but the residual doubles from the disconnected
control. The exact rational upper values are stored in each certificate JSON;
the question marks above avoid pretending decimal rounding is exact.

Run the sweep with `python -m experiments.marginal_transfer --delta 1/10`.
Replay without numerical libraries with `python -S -c '...'` calling
`marginal_symbolic.verify` on each JSON. Tests cover graph connectivity,
replay, and independent variational positivity.
