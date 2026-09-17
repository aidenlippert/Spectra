# Prior art: local response certificates and their limits

Date: 2026-09-16

## Local-cluster and marginal certificates

The standard Anderson/cluster lower bound decomposes a local Hamiltonian into patches and lower-bounds each patch by its smallest eigenvalue, with coefficients chosen so every original term is counted correctly. This is already used as a practical lower-bound companion to tensor-network variational upper bounds. A recent primary treatment is Eisert, *Lower bounds to variational problems with guarantees*, [Phys. Rev. A 113, 022214 (2026)](https://doi.org/10.1103/cz6k-y46r). It also develops a hierarchy of semidefinite relaxations inspired by quantum marginals. This means that “overlapping local blocks plus an SDP dual” is established territory by itself.

The quantum-marginal side is also mature. Yu, Simnacher, Wyderka, Nguyen, and Gühne give an SDP hierarchy whose levels converge to the compatible pure-state marginal set: [Nature Communications 12, 1012 (2021)](https://doi.org/10.1038/s41467-020-20799-5). The dual of a local-marginal relaxation is naturally a Hamiltonian lower-bound witness. Its limitation is equally important: fixed local marginal order is not generally exact for correlated systems, and higher levels grow in the marginal/extension size.

## Positivity-preserving tensor representations

Locally purified tensor networks represent a density operator as a partial trace of a pure tensor network, guaranteeing positivity and allowing trace-norm error control in suitable simulations: Werner et al., [A positive tensor network approach for simulating open quantum many-body systems](https://arxiv.org/abs/1412.5746). This is useful infrastructure for positive response objects or thermal states. It is not automatically a lower-bound certificate: positivity of a represented density matrix gives a variational upper expectation, while a certified lower bound still requires a dual witness or an operator inequality. An MPO being numerically positive, or an MPDO having a purification, does not certify that it encloses the exact eliminated response.

## Local Schrieffer--Wolff bounds

Bravyi, DiVincenzo, and Loss prove linked-cluster structure and truncation-error bounds for Schrieffer--Wolff effective Hamiltonians in local spin systems: [Schrieffer--Wolff transformation for quantum many-body systems](https://arxiv.org/abs/1105.0675). The theorem gives controlled locality/order dependence in a perturbative, gapped setting. It does not cover arbitrary strongly correlated fermions at a critical point or provide a nonperturbative two-sided Loewner enclosure. Any local response proposal using perturbative denominators must state its gap and coupling regime explicitly.

## Approximate Markov structure

Chen and Rouzé prove quasi-local recovery maps and conditional-mutual-information decay for bounded-degree quantum Gibbs states at arbitrary temperature: [Quantum Gibbs states are locally Markovian](https://arxiv.org/abs/2504.02208). This is a promising source of compact environment messages, but the theorem concerns local state recovery and preparation under its stated assumptions. It does not turn a local recovery map into a ground-energy lower bound or an exact Schur-response enclosure. Kato and Brandão's earlier one-dimensional result is [Quantum Approximate Markov Chains are Thermal](https://arxiv.org/abs/1609.06636); its one-dimensional scope is a warning against presenting the same argument as a 2D strongly correlated breakthrough.

## Why fixed block size cannot generally deliver arbitrary total-energy precision

Let (H_L=\sum_{i=1}^{L}h_i) be a translation-invariant finite-range chain and let a fixed patch relaxation produce an energy-density lower bound (e_r\le e_0) from patches of diameter (r). Composing the patch certificates gives

\[
E_0(H_L)\ge L e_r+O(r),
\]

so a nonzero density defect (e_0-e_r) becomes an (O(L)) total-energy error. A fixed (r) can achieve arbitrary *energy-density* accuracy only if one proves (e_r=e_0) or a special exactness property. Critical systems are especially unsuitable for assuming this: their long correlation length and finite-size corrections are not determined by a bounded patch alone. To demand a fixed absolute total-energy error as (L\to\infty), the certificate order must generally grow, or it must carry a nonlocal/collective boundary response.

This is why the proposed Schur response is materially different from a plain patch sum: it is intended to carry the boundary coupling and its induced energy. But the boundary response must be represented with a matrix inequality, not merely one scalar residual or a local patch ground energy.

## What could still be new

The defensible novelty target is a theorem or algorithm with all of these properties:

1. It starts from a declared interacting fermion Hamiltonian and constructs a retained-sector response directly from local/tensor contractions.
2. It certifies the omitted contribution in Loewner order, or proves an observable-specific bound sufficient for the target energy.
3. It allows charge transfer and collective cross-cluster correlations.
4. Its error and representation size compose through a second elimination.
5. It returns a precise failure certificate when locality, gap, or Markov assumptions are insufficient.

None of the cited cluster SDP, positive tensor-network, local-SW, or Markov results supplies this conjunction. That is a scoped research opportunity, not evidence that a general theorem exists. The most useful first test is a fully coupled non-frustration-free and then critical benchmark: compare a fixed-order patch bound, a patch-plus-boundary response, and the exact/block-Lanczos reference while measuring both total-energy error and representation cost.

## Implementation note

For a numerical prototype, Clarabel documents the cone representation used for positive-semidefinite constraints, including the scaled upper-triangular vectorization expected by its PSD cone: [Clarabel supported cone types](https://clarabel.org/stable/api_cone_types/). This can make a local SDP implementation cleaner, but a floating-point solver output still requires exact or interval post-verification before it becomes a Spectra certificate.

## Sources and inspection limits

I inspected the primary article records/abstracts and the statements available at the linked journal or arXiv pages. I did not reproduce the proofs, search patents, or perform an exhaustive literature review. The conclusions are therefore a scoped novelty and feasibility audit.
