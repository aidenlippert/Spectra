# H-only commutator dictionaries: corrected bounded experiments

The root implementation keeps the entire symmetry-partitioned quadratic
baseline and enriches linear charge blocks with cubic components of [H,a_i]
and their adjoints. The omitted linear component is already in the baseline.
A second variant resolves the cubic component by its creation-orbital index.
Those components sum to the first variant, so the exact enriched span includes
it. Power-of-two rational rescaling controls coefficients without changing
these spans. Discovery uses H alone, never an FCI state or a source Gram.

Every upper-triangle Gram coefficient is expanded exactly: p_i^dagger p_i,
and p_i^dagger p_j plus its adjoint for i<j. Every resulting CAR monomial,
including degree six, enters the equality map. The fixed-number ideal uses
all symmetry-compatible Hermitian body-two multipliers. All charge blocks are
exported by rational linear combination into existing word dictionaries and
a global integer factor denominator. The original exact verifier remains the
acceptance gate. The implementation adds no trusted polynomial-factor schema.

## Results

All six corrected lower/upper intervals were independently replayed with
python -S against the same frozen rational Hamiltonians.

| Fixture / dictionary | Gram entries | Map nonzeros | Build / solve / export seconds | Exact width Ha |
|---|---:|---:|---:|---:|
| H4 quadratic baseline | 1,136 | 1,280 | .080 / .178 / .064 | .00448024 |
| H4 summed commutator | 1,232 | 20,580 | .830 / .218 / .378 | .00448540 |
| H4 creator channels | 3,152 | 45,440 | 1.221 / .786 / .437 | .00356381 |
| H6 quadratic baseline | 5,724 | 6,180 | .469 / 1.885 / .333 | .01300618 |
| H6 summed commutator | 5,940 | 433,284 | 15.440 / 4.796 / 5.157 | .01300710 |
| H6 creator channels | 16,020 | 1,241,904 | 25.300 / 28.118 / 5.933 | .01220835 |

None reaches .0016 Ha. The summed commutator gives effectively the same
numerical objective as the quadratic baseline on these fixtures; this is
not a theorem of algebraic redundancy. Creator resolution improves the bound
modestly. Several Clarabel results are optimal_inaccurate, so no cone optimum
or universal impossibility is inferred from solver status.

For comparison the successful full H6 cubic dictionary built 217,268 Gram
entries and 113,614 map nonzeros. The creator-resolved dictionary therefore
reduces Gram entries about 13.6-fold while increasing map nonzeros about
10.9-fold. Its 818,628 exact polynomial word-pair products are explicitly
counted. This is a failed accuracy candidate and a warning against measuring
only PSD block dimension: Hamiltonian-weighted operators can densify the
coefficient map. No generic efficient-accuracy condition follows.

Exact artifacts, intervals, and costs are in
results/certificate_scaling/commutator_dictionary/h{4,6}_{corrected,creator_channels,quadratic_baseline}.

The initial agent draft under h4 omitted square coefficients from its SDP,
exported only one charge block, and lacked the requested baseline/ideal.
Its weak certificate can remain a valid lower bound, but its numerical result
is NOT a valid comparison for the declared dictionary. Only the corrected
root runs above count. Four controls test the density commutator, independent
polynomial square expansion, degree-six retention, exact factor export, and
creator-channel span/charge preservation. A read-only independent audit found
no algebraic error in the root implementation within this stated scope.
