# Observable-only control: established scope

The accepted implementation is collective_control.py. It retains the actual molecular Hamiltonian and supplied MPS, but uses a short, strong constant pulse. Its exact CAR algebra, Pauli-product norm bounds and scalar MPS contractions do not construct a state trajectory or enumerate a particle sector.

On H8, v=64 Ha, u=0, and T=25/1024 atomic units give the robust interval [-1.98454768699, -0.79533721980]. The original maximum amplitude was 0.5 Ha. This is a separate regime, not a solution of the original four-piece protocol.

collective_obstruction.py proves that the global commutator-sum envelope cannot certify the target for any duration with constant W-only control, u=0, and |v|<=0.5. This limits that proof rule, not physical reachability or the original four-phase family. Initial moments and model/state bindings are checked again.

An earlier diagnostic divided already scaled coefficients by the Hamiltonian denominator again and incorrectly reported nearly vanishing commutators. That result is invalid. The repaired coefficient-only diagnostic gives CAR l1 upper bounds 63.72438317538 and 62.860446170876 Ha for [H,D] and [H,W]. The accepted construction uses tighter exact Pauli-product bounds for [H,D] and [H,A]: 19.128975101863 and 18.676738267037 Ha. An inconclusive upper bound is not an obstruction.

An adjoint-weighted residual remains a possible refinement. A scalar residual contraction alone is insufficient. For rho_dot=-i[H,rho], the backward-observable defect is O_dot+i[H,O]. Its expectation against an approximate state needs an additional allowance for that state's error. The complete forward/backward identity is in DERIVATION.md; no accepted weak-control result is claimed from it.
