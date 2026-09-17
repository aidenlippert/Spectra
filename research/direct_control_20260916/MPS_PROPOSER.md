# Direct tensor-state investigation

The exact accepting path reuses research/correlated_pair_20260913/mps_exact.py. It validates the supplied integer charge-flow MPS and computes moments by integer transfer contractions. The H8 input has maximum bond 144 and 10,093 nonzero tensor entries. Its normalized initial contrast is 1.390133812147592... . The strong-control certificate needs only this moment; the conjugate quadrature has zero expectation because the tensors are real.

validated_tensor.py implements proposed tensor canonicalization, operator action and compression with a posteriori residual allowances. QR and SVD are proposers: actual reconstruction and isometry errors are charged. These diagnostic bounds assume the stated IEEE binary64/BLAS rounding model. They are distinct from the pure rational accepting path for strong control.

Neither tensor action experiment certified the original dynamics. The first bond-64 action bound was about 1.34 million; canonicalizing the operator and state reduced it to 192.9873, still unusable. This is failure of these sufficient bounds, not a theorem that the best tensor approximation fails.

The orbital-order and local-basis probes operate directly on the MPS. The local rotation uses 28 exactly orthogonal rational Givens rotations, implemented as local fermionic gates with separately bounded numerical application error. The one-particle rotation proposal is inherited from prior work; cold discovery was not demonstrated here.

Checked Schmidt-tail diagnostics give state-approximation lower bounds for a declared basis/order, not error bounds for one observable. Valid receipts have _fixed in their names. The first square-root-scale receipt is invalid and excluded.

An early mps_propose.py misread the CAR encoding and tensor normalization. That implementation and its self-confirming test were removed. Its reported moments must not be cited as evidence. The independently tested exact backend and final tensor diagnostics are retained.
