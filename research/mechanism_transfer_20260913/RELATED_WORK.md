# Direct comparison with the two relevant constructions

Rubin, Low and DePrince's paper is indeed at version 2, revised September 9,
2026. Its Eq. 28 uses linear and quadratic ladder-operator generators;
Eqs. 29--40 organize particle, hole and spin blocks. Section III.3 restricts
to spin-summed density generators, and Section IV connects weighted SOS with
v2RDM constraints. Its particle/spin constraints are part of the model, not
an optional numerical convention. [Paper, v2](https://arxiv.org/html/2602.05069v2)

The present code retains the quadratic baseline and adds restricted
degree-three generators Q[k,s,t]a[i]. Their tied anticommutator cancels
degree-six terms exactly. It uses a total-number ideal, without imposing a
fixed total-spin target, and applies a separately verified Hamiltonian tail.
Consequently, the paper's rank-2 numerical results are not a matched control
for this program. Its spin-free density algebra also differs from retaining
all four spin components inside cubic operators. This is a comparison of
the explicit algebras, not a claim that cubic positivity itself is new.
[Relevant definitions, Eqs. 28 and 52--79](https://arxiv.org/html/2602.05069v2#S3)

Ahmadi, Dash and Hall's Section 3.1 describes adding PSD atoms from negative
eigenvectors of a dual matrix; Section 4 discusses normalizing the violation
by a matrix norm. Full PSD feasibility, rather than a sampled absence of
violations, is the condition that closes their larger SDP comparison.
[Column-generation paper, Sections 3--4](https://arxiv.org/pdf/1512.05402)

Our earlier ordinary eigenvector pricing follows that established principle.
The fixed-N trace metric changes the geometry used to choose fermionic
operators. Coupled blocks then permit cross terms between selected operators;
they are more general than independent nonnegative weights on the same
rank-one atoms. This remains a restricted SOS optimization strategy. Neither
normalized eigenvectors nor adding useful directions alone establishes
novelty. The present pass instead tests a fixed, Hamiltonian-derived rule
without per-instance column generation.

The candidate's alpha=1 spin-summed contraction satisfies the ordinary
Heisenberg equation-of-motion algebra, verified here by exact CAR reduction.
No priority or novelty claim is made for that identity. A useful contribution
would require evidence that a compact rule preserves substantial certified
gain at lower total cost, and that the benefit transfers. Operator projection
coverage alone does not supply that evidence.
