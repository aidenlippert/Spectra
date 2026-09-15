"""Run the bounded keystone experiment and emit reproducible evidence."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform
from time import perf_counter
import numpy as np
from experiments import discovery
from experiments.pauli import closure, commutator_i, generator_matrix
from experiments.certificates import PauliTerm, certify_hamiltonian, check_certificate
from experiments.stabilizer import diagnose_css


@lru_cache(maxsize=4096)
def dense_pauli(label):
    mats = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
            'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
    if not label or len(label) > 6 or any(c not in mats for c in label):
        raise ValueError('dense validator limited to 1..6 qubits')
    out = np.ones((1, 1), dtype=complex)
    for letter in label:
        out = np.kron(out, mats[letter])
    return out


def dense_hamiltonian(h):
    return sum(float(c) * dense_pauli(p) for p, c in h.items())


def trajectory_check(h, basis, seed, oracle_h=None, model_error=Fraction(0)):
    """Independent dense Schrödinger validation of reduced Heisenberg dynamics."""
    rng = np.random.default_rng(seed); dim = 2 ** len(next(iter(h)))
    psi = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    psi /= np.linalg.norm(psi)
    z = np.array([np.vdot(psi, dense_pauli(p) @ psi).real for p in basis])
    largest = 0.; largest_skew = 0.; epsilon_integral = 0.
    # Piecewise changed coefficient controls withheld from the identification data.
    for _ in range(4):
        factors = {p: Fraction(int(rng.integers(5, 17)), 10) for p in sorted(set(h) | set(oracle_h or {}))}
        controlled = {p: c * factors[p] for p, c in h.items()}
        physical = controlled if oracle_h is None else {p: c * factors[p] for p, c in oracle_h.items()}
        _, matrix, residual = generator_matrix(controlled, basis)
        a = np.array(matrix, dtype=float)
        largest_skew = max(largest_skew, float(np.max(np.abs(a + a.T))))
        dt = 0.2
        epsilon_integral += dt * float(np.linalg.norm([float(x) for x in residual.values()]))
        # All undisplayed dictionary coefficients also obey a factor <=1.6.
        epsilon_integral += dt * 2 * float(np.sqrt(len(basis))) * 1.6 * float(model_error)
        w, v = np.linalg.eigh(dense_hamiltonian(physical))
        psi = v @ (np.exp(-1j * w * dt) * (v.conj().T @ psi))
        wa, va = np.linalg.eigh(1j * a)
        z = (va @ (np.exp(-1j * wa * dt) * (va.conj().T @ z))).real
        truth = np.array([np.vdot(psi, dense_pauli(p) @ psi).real for p in basis])
        largest = max(largest, float(np.linalg.norm(truth - z)))
    return {'max_observable_vector_error': largest, 'skew_symmetry_error': largest_skew,
            'analytic_error_bound': epsilon_integral,
            'floating_validation_tolerance': 1e-10, 'pieces': 4, 'total_time': 0.8,
            'numerically_within_bound': bool(largest <= epsilon_integral + 1e-10),
            'model_norm_error': str(model_error),
            'scope': ('oracle dynamics with conditional Hamiltonian error bound' if oracle_h is not None else 'nominal dense numerical validation; exact closure algebra is the certificate')}


def representation(h, seed, previous=None):
    controls = [{p: Fraction(1)} for p in h]
    basis, complete = closure(controls, h.keys(), budget=64)
    basis = sorted(basis)
    if not complete:
        raise RuntimeError('unexpected complete-basis budget failure')
    residual = {}
    for p in h:
        residual[p] = {o: str(c) for o, c in generator_matrix({p: 1}, basis)[2].items() if c}
    if any(residual.values()):
        raise RuntimeError('closure verification failed')
    broken = {}
    if previous:
        for p in h:
            for o, c in generator_matrix({p: 1}, previous)[2].items():
                if c:
                    broken.setdefault(p, {})[o] = str(c)
    cap_basis, cap_complete = closure(controls, h.keys(), budget=2)
    return {'basis': basis, 'basis_size': len(basis), 'full_nonidentity_dimension': 63,
            'exactly_closed_for_each_control': True,
            'coefficient_transfer': 'same operator support permits arbitrary piecewise real coefficients',
            'prior_basis_failure_witnesses': broken,
            'budget_two_abstained': not cap_complete,
            'budget_two_returned_basis': sorted(cap_basis),
            'verification_commutator_calls': len(h) * len(basis),
            'trajectory': trajectory_check(h, basis, seed)}


def energy_evidence(h, model_error, oracle_h):
    start = perf_counter()
    terms = [PauliTerm(p, c) for p, c in sorted(h.items())]
    width = len(terms[0].pauli)
    def value(bits):
        return sum((t.coeff * (-1) ** sum(bit for p, bit in zip(t.pauli, bits) if p == 'Z')
                    for t in terms if all(p in 'IZ' for p in t.pauli)), Fraction(0))
    states = list(product((0, 1), repeat=width))
    state = min(states, key=value)
    cert = certify_hamiltonian(terms, state=state, omitted_coefficients=[model_error])
    valid = check_certificate(terms, cert)
    if not valid:
        raise RuntimeError('independent energy certificate rejected')
    lower, upper = Fraction(cert['lower']), Fraction(cert['upper'])
    reference = float(np.linalg.eigvalsh(dense_hamiltonian(oracle_h))[0])
    baseline = -sum(map(abs, h.values()), Fraction(0)) - model_error
    numerical_bracket = float(lower) - 1e-10 <= reference <= float(upper) + 1e-10
    if not numerical_bracket:
        raise RuntimeError('oracle spectrum contradicts claimed interval')
    return {'certificate': cert, 'independent_rational_check': valid,
            'basis_states_evaluated': len(states), 'termwise_lower': str(baseline),
            'lower_improvement': str(lower - baseline), 'absolute_interval_width': str(upper - lower),
            'oracle_ground_energy_numeric': reference, 'oracle_numeric_bracket_check': numerical_bracket,
            'wall_seconds': perf_counter() - start,
            'scope': 'bound on physical H conditional on finite dictionary and exact sensor error intervals'}


def integrate(seed):
    result = discovery.benchmark(seed)
    previous = None
    for stage in result['stages']:
        stage_rep = None
        for ci, ctx in enumerate(stage['contexts']):
            h = {p: Fraction(c) for p, c in ctx['learned_hamiltonian'].items()}
            true_h = {p: Fraction(c) for p, c in ctx['oracle_hamiltonian'].items()}
            rep = representation(h, seed + 100 + stage['stage'] * 10 + ci, previous)
            ctx['representation'] = rep
            if ctx['conditional_model_norm_error'] is None:
                ctx['energy'] = {'abstained': True, 'reason': ctx['interval_status']}
            else:
                eta = Fraction(ctx['conditional_model_norm_error'])
                ctx['energy'] = energy_evidence(h, eta, true_h)
                ctx['oracle_trajectory'] = trajectory_check(h, rep['basis'], seed + 500 + stage['stage'] * 10 + ci, true_h, eta)
                if not ctx['oracle_trajectory']['numerically_within_bound']:
                    raise RuntimeError('oracle dynamics exceeded uncertainty-aware bound')
            if ci == 0:
                stage_rep = rep['basis']
        previous = stage_rep
    return result


def source_hashes(root):
    return {str(p.relative_to(root)): sha256(p.read_bytes()).hexdigest()
            for folder in ('experiments', 'tests', 'research')
            for p in sorted((root / folder).glob('*')) if p.is_file() and p.suffix in ('.py', '.md')}


def summarize(data):
    contexts = [c for run in data['runs'] for stage in run['stages'] for c in stage['contexts']]
    seq = [c['methods']['sequential'] for c in contexts]
    scratch = [c['methods']['from_scratch'] for c in contexts]
    return {'hamiltonian_contexts': len(contexts),
            'successful_sequential_support_recoveries': sum(x['oracle_support_correct'] and not x['abstained'] for x in seq),
            'rational_energy_certificates_verified': sum(c['energy'].get('independent_rational_check', False) for c in contexts),
            'max_heldout_derivative_error': max(x['holdout_max_error'] for x in seq),
            'max_dense_dynamics_error': max(c['representation']['trajectory']['max_observable_vector_error'] for c in contexts),
            'max_oracle_dynamics_error': max(c['oracle_trajectory']['max_observable_vector_error'] for c in contexts),
            'all_oracle_dynamics_within_certified_bound': all(c['oracle_trajectory']['numerically_within_bound'] for c in contexts),
            'learned_basis_sizes': sorted(set(c['representation']['basis_size'] for c in contexts)),
            'sequential_correlation_multiply_adds': sum(x['cost']['correlation_multiply_adds'] for x in seq),
            'scratch_correlation_multiply_adds': sum(x['cost']['correlation_multiply_adds'] for x in scratch),
            'equal_sensor_record_budget': all(c['methods']['sequential']['sensor_records_used'] == c['methods']['from_scratch']['sensor_records_used'] for c in contexts),
            'out_of_family_rejections': sum(r['misspecified_status'] == 'out_of_family_or_noise_bound_violated' for r in data['runs']),
            'unidentifiable_rejections': sum(r['unidentifiable_abstained'] for r in data['runs']),
            'compounding_scientific_capability_demonstrated': False,
            'interpretation': 'support reuse saves some greedy search; fixed probes, conventional full-dictionary solve and retrieval on shared support explain performance; no novel scientific operation demonstrated'}


def write_report(data, path):
    s = data['summary']
    lines = ['# Bounded keystone experiment — measured result', '',
             'This run exercises support identification, exact operator closure, conditional energy certificates and finite-patch failure detection. It does not establish compounding scientific intelligence.', '',
             f"Generated: {data['generated_at']}. Seeds: {data['seeds']}. NumPy {data['numpy']}; Python {data['python']}.", '',
             '| Check | Result |', '|---|---|']
    for label, key in [('Evaluated Hamiltonian contexts', 'hamiltonian_contexts'),
                       ('Sequential support recoveries', 'successful_sequential_support_recoveries'),
                       ('Rational energy certificates verified', 'rational_energy_certificates_verified'),
                       ('Worst held-out derivative error', 'max_heldout_derivative_error'),
                       ('Worst dense observable-vector discrepancy', 'max_dense_dynamics_error'),
                       ('Worst simulator dynamics error with learned parameters', 'max_oracle_dynamics_error'),
                       ('Simulator dynamics inside conditional bound', 'all_oracle_dynamics_within_certified_bound'),
                       ('Reduced basis sizes', 'learned_basis_sizes'),
                       ('Out-of-family rejection cases', 'out_of_family_rejections')]:
        lines.append(f'| {label} | {s[key]} |')
    lines += ['', '## What the cost comparison says', '',
              f"Greedy correlation multiply-add counts: sequential {s['sequential_correlation_multiply_adds']:,}; scratch {s['scratch_correlation_multiply_adds']:,}. These count one part of the computation, not total speedup. Sensor data budgets are equal. Full-dictionary least squares avoids greedy search, and retrieval can refit already-known supports. Physical shots, derivative estimation and apparatus costs are absent from the ideal sensor model.", '',
              'The first stage identifies two coefficients. Two later stages introduce one interaction each. Each stage also tests a coefficient vector generated by a declared deterministic rescaling and disjoint measurement-observable probes. Base Hamiltonians repeat across seeds; these are 18 evaluations, not 18 statistically independent physical systems. This is instance transfer within a public two-body Pauli family; it does not discover a new vocabulary or choose an informative experiment.', '',
              '## Certificates and boundaries', '',
              'For each coefficient, rational sensor inequalities enclose its value conditional on the declared family and bounded sensor error. Their summed worst errors bound the Hamiltonian operator norm. The energy certificate combines anticommuting Pauli groups with an exact computational-basis trial state and widens both ends by that norm bound. A separate checker verifies the supplied partition and rational square-root inequalities. Dense diagonalization is only a numerical cross-check.', '',
              'Closure is checked algebraically for every independently controlled Hamiltonian term. The reduced generator is skew-symmetric in normalized Pauli coordinates, so residual error accumulates linearly in time. New couplings can invalidate the previous basis; the run records omitted-operator witnesses. This exact statement covers the learned model. The actual simulator is also checked against an uncertainty-aware bound: an omitted Hamiltonian of norm eta contributes at most 2 sqrt(k) eta to the vector residual, multiplied by the declared control scale (at most 1.6). A zero nominal closure residual does not erase model uncertainty.', '',
              'The CSS diagnostic detects the manufactured X1X2/X2X3 bridge obstruction and passes when the buffer includes the bridge. These are diagnostic examples, not instances of the 2026 recursive memory code.', '',
              '## Verdict', '',
              ('**Pipeline checks passed; compounding scientific capability was not demonstrated.**' if s['pipeline_gate_passed'] else '**Pipeline gate failed; examine the recorded failures.**') + ' No superiority over a strong conventional method is claimed. The next research gate is a learned representation or experiment rule that beats equally equipped cache and conventional baselines on new coupled systems after acquisition and verification costs are charged.', '',
              'Full rational certificates, errors, baselines, source hashes and cost records are in [keystone_results.json](keystone_results.json).']
    path.write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seeds', type=int, nargs='+', default=[7, 11, 23])
    parser.add_argument('--output', type=Path, default=Path('results'))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    start = perf_counter()
    bad = diagnose_css([(1, 1, 0), (0, 1, 1)], [], 3, {0, 2}, {0, 2})
    good = diagnose_css([(1, 1, 0), (0, 1, 1)], [], 3, {0, 2}, {0, 1, 2})
    data = {'generated_at': datetime.now(timezone.utc).isoformat(), 'seeds': args.seeds,
            'python': platform.python_version(), 'numpy': np.__version__,
            'runs': [integrate(seed) for seed in args.seeds],
            'stabilizer_diagnostic': {'missing_bridge_rejected': not bad[0].passed,
                                     'witness_bits': bad[0].witness, 'restored_bridge_passed': good[0].passed}}
    data['summary'] = summarize(data)
    s = data['summary']
    s['pipeline_gate_passed'] = (s['successful_sequential_support_recoveries'] == s['hamiltonian_contexts']
        and s['rational_energy_certificates_verified'] == s['hamiltonian_contexts']
        and s['max_heldout_derivative_error'] <= 0.006
        and s['all_oracle_dynamics_within_certified_bound']
        and s['out_of_family_rejections'] == len(args.seeds)
        and s['unidentifiable_rejections'] == len(args.seeds)
        and data['stabilizer_diagnostic']['missing_bridge_rejected']
        and data['stabilizer_diagnostic']['restored_bridge_passed'])
    data['elapsed_seconds'] = perf_counter() - start
    data['source_sha256'] = source_hashes(root)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'keystone_results.json').write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    write_report(data, args.output / 'REPORT.md')
    print(json.dumps(data['summary'], indent=2))
    if not s['pipeline_gate_passed']:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
