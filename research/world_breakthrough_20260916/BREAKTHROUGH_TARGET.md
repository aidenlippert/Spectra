# A field-level target for the next research direction

Assessment date: 16 September 2026. This is a research recommendation, not an
achieved result or a prediction that a particular approach will work.

## Recommendation

Pursue a constructive, controlled theory of strong electronic correlations:
identify which collective responses must be retained, construct them without a
full-state teacher, and bound the effect of the remaining degrees of freedom
at a cost useful in difficult two-dimensional, low-temperature regimes.

The mathematical target is a new tractable physical class plus a constructive
algorithm. The algorithmic target is useful observables at demonstrably better
total cost and controlled error. The physics target is a previously unresolved
statement about a collective phase or response. Any one of these can be a
substantial scientific advance; a new material is not a prerequisite.

## The existing bar

Provably efficient construction of ground states already exists for gapped
one-dimensional local Hamiltonians. Rigorous renormalization algorithms also
cover extensions within one dimension. The new contribution must extend scope,
cost or guarantees beyond that established class, rather than merely rename
recursive reduction. [Landau, Vazirani and Vidick](https://arxiv.org/abs/1307.5143),
[Arad, Landau, Vazirani and Vidick](https://arxiv.org/abs/1602.08828).

Efficient algorithms for short-time local-observable dynamics also exist under
specified locality and initial-state assumptions. A short pulse that suppresses
the influence of the interacting drift therefore does not by itself establish
a new general mechanism for hard quantum dynamics.
[Wild and Alhambra](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.4.020340).

Small representation size alone is insufficient: general PEPS contraction is
computationally hard, and area laws alone do not imply efficient classical
descriptions in higher dimensions. These are worst-case results, not a claim
that all useful physical tensor networks are intractable.
[Schuch et al.](https://arxiv.org/abs/quant-ph/0611050),
[Ge and Eisert](https://arxiv.org/abs/1411.2995).

## A precise candidate mathematical problem

For a retained projector P, Q=I-P, and energy E below the spectrum of QHQ, the
exact eliminated response is

\[
\Sigma(E)=PHQ(QHQ-E)^{-1}QHP.
\]

The retained eigenvalue equation involves PHP minus this response. The formula
and its use in renormalization are established; their rediscovery would not
be the breakthrough. [Bach, Ballesteros and Fröhlich](https://arxiv.org/abs/1308.0911).

The proposed new result would supply efficiently constructible bounds

\[
\Sigma_-(E)\preceq\Sigma(E)\preceq\Sigma_+(E),
\qquad\|\Sigma_+(E)-\Sigma_-(E)\|\le\delta,
\]

over the required energy interval, while keeping the retained operators and
their certified responses manageable after further elimination steps. The
projectors, the eliminated-sector spectral premise, induced interactions and
all construction costs must be provided, not assumed after solving the full
problem. Observable-specific weighted bounds may be more economical than
uniform operator bounds, but must include their domain and omitted-state error.

A meaningful complexity theorem would tie these costs to independently
checkable physical structure. Defining a small complexity parameter only after
the exact solution is known would not address discovery. The criterion must
apply to a physically significant class beyond an already easy high-temperature,
weak-coupling, short-time, or fixed-width regime. This proposed theorem has not
been proved here, and there is no claim that one compact response form must
work near every critical point.

## A discriminating physics target

A strong candidate is controlled low-temperature, thermodynamic-limit phase
behavior of a declared doped repulsive square-lattice Hubbard model, including
the distinction between the nearest-neighbor and next-nearest-neighbor variants.
The question is whether pairing, stripe order, and competing states persist
when width, environment and approximation errors are controlled. In two
dimensions, a finite-temperature superconducting claim needs appropriate BKT
or equivalent scaling diagnostics, not just a nonzero finite-system pairing
correlator. This is a proposed target, not an assertion that every part of its
phase diagram remains unknown.

There are already strong competing results. Xu et al. reported ground-state
superconductivity coexisting with stripes in the model with next-nearest
hopping. Their calculation resolves small energy differences using several
sizes and boundary conditions. [Published in Science, 2024](https://arxiv.org/abs/2303.08376).

The 2025 preprint by Wang and Devereaux reports finite-temperature signatures
of underlying electron-doped d-wave superconductivity using determinant Monte
Carlo. This is a significant baseline, not proof that every finite-temperature
transition is settled. [Primary preprint](https://arxiv.org/abs/2510.16616).

The 2025 enhanced-XTRG preprint reports finite-system Hubbard calculations down
to approximately 0.004t. Consequently, merely reaching a low temperature or
showing pairing on a small cylinder would not define a new field-level result.
[Zhang and von Delft](https://arxiv.org/html/2510.25022v1).

An independently reproducible, genuinely new controlled phase-boundary result
or a demonstrated algorithmic change in accessible size/temperature/accuracy
would be meaningful. Independent numerical and quantum-simulator comparisons
strengthen the claim where available; lack of an immediate laboratory experiment
does not invalidate a mathematical breakthrough.

## Scope and decision

Fixed-basis electronic structure is QMA-complete in the worst case. This is
not a proof that useful broad physical classes lack efficient algorithms, nor
a proved separation of complexity classes. It does mean that a universal cheap
solver for every Hamiltonian is a much stronger claim than a new physical
tractability theorem. [O'Gorman et al.](https://arxiv.org/abs/2103.08215).

The most promising direction in this assessment is a mathematical result that
enables a new algorithm and then settles a consequential physics question.
H6/H8 can remain correctness tests. They should not determine the ambition or
serve as the sole evidence of novelty. A successful model result would still
need further modeling and experimental work before it predicts a real material
or a synthesis policy.
