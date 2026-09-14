# Separating invariant certificate feasibility from numerical discovery

The H8 SCS300 spin-invariant discovery run produced a weak exported certificate:
its independent spectral interval width is0.699480442455 Ha. A control now
transforms a previously passing H8 certificate into the same spin-invariant
SOS word family. It tests representational feasibility, explicitly reusing an
existing certificate. Its cost cannot count as discovering a certificate from
H alone.

## Exact representation identity

Let p_(a,r) be the divided-power spin-j descendants, with invariant integer
weights w_r proportional to1/binom(2j,r). For
p=sum_(a,r)c_(a,r)p_(a,r), Haar averaging gives

\[
\mathcal T(p^\dagger p)=
\sum_{a,b}\left[\sum_r\frac{c_{a,r}^*c_{b,r}}{(2j+1)w_r}\right]
\sum_s w_s p_{a,s}^\dagger p_{b,s}.
\]

The bracket is PSD. Multiplicity copies need not be orthogonal: spin acts
identically on the descendant index of every copy. A nonorthogonal change of
multiplicity coordinates commutes with that action. No inverse multiplicity
metric belongs in the formula. An initial agent objection claiming otherwise
was retracted after the action-matrix derivation and exact controls.

For independent polynomial checks, the adjoint Casimir is
C(P)=[Sz,[Sz,P]]+([S+,[S-,P]]+[S-,[S+,P]])/2.
On an even polynomial of degree at most six,

\[
\mathcal T=\prod_{j=1}^{3}\left(I-\frac{C}{j(j+1)}\right)
\]

projects onto the spin-invariant component; factors above half the actual
degree can be omitted. It also satisfies
T((Nhat-N)X)=(Nhat-N)T(X).
`spin_twirl.py` implements this exact rational polynomial projector.

## Executed transformation

`spin_twirl_control.py` matches each source block to a unique current
quadratic/cubic word family. It exactly inverts each local spin-pattern basis
(dimension at most8), then uses floating Gram contractions and the existing
rational exporter. The standard exact verifier re-expands the final SOS against
the ORIGINAL Hamiltonian, charging all numerical error. It twirls the number
multiplier exactly and retains compatible parity sectors. The original scalar
bound is retained.

The control saves proposed Gram matrices for diagnostics, reports
`source_certificate_used=true`, and records the source certificate hash.
It refuses unsupported/ambiguous source blocks. The particular prior H6
certificate contains explicit linear words outside this degree-2/3 control;
H4 and H8 source blocks fit unambiguously.

Five exact projector tests cover all four spin ranks, nonorthogonal copies,
Schur contraction, direct commutation with spin generators, idempotence, wrong
weights, degree refusal and number-ideal equivariance. A separate H4 integration
test compares the exported SOS with the independent exact Casimir projection
of the complete source SOS, bounds their rational coefficient difference, and
checks preservation of the original H. The full focused suite now has58 passes.

H4 transformed independently passes at2.51996453e-7 Ha, with0.905s transformation
plus export/replay in the source process (source discovery remains additional).
H8 transformed output took106.274s: coordinate/Gram transform2.572s,
export11.631s, exact replay92.071s. It has320,192 invariant Gram entries,
6,255 factor rows,1,252,050 nonzero factor entries and7,403,714 certificate bytes.
Its coefficient-l1 lower bound is-9.258278219641182 Ha. A separate spectral
residual proof gives-9.256055015126634 Ha and costs95.480s; independent interval
replay is pending as of this note.

This control does not establish a chemical-accuracy optimum for the solver's
coefficient-l1 objective, nor fast discovery in the invariant cone. The spectral
residual proof can be substantially sharper than that objective. Source
certificate discovery and FCI upper-witness discovery remain separate costs.
Artifacts are under `results/certificate_scaling/spin_irrep/h4_twirl_control`,
`h8_twirl`, and `h8_twirl_spectral`.
