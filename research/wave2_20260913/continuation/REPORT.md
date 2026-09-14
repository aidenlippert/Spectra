# Original H4 homotopy continuation audit

The prior H4 certificate supplies exact endpoint bounds
`-3.667000108966591573 <= E0 <= -3.6669999559641995`, width
`1.530023922e-7`. The requested path uses `lambda=0,1/4,1/2,3/4,1`.

The complete rational CAR fixture was reconstructed with the independent
70×70 oracle. Exact Gershgorin lower bounds and rational `e_0` trial uppers
were evaluated at all five points. The exact maximum absolute row sum of the
off-diagonal is `1706063013273/10^12 = 1.706063013273`, giving the valid global
fallback `L(lambda)>=L(0)-lambda*||V||_infty`; it is too loose to improve the
endpoint certificate. The relative PSD scan failed for small epsilon and
passed only the vacuous `eta=1, epsilon=1` case. LDL refuses the coupled
zero-pivot indefinite matrix.

This is a quantitative failure of the available evidence, not evidence that
the original molecular homotopy fails. A stronger nonlocal relative bound is
required for improvement.
