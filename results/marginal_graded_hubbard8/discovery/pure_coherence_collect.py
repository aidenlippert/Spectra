"""Collect source-bound pure-coherence energy, family and extension evidence."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/pure_coherence'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path):
    return str(path.relative_to(ROOT))


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def accepted(path):
    data = read(path)
    if data.get('accepted') is not True:
        raise ValueError('Required accepted receipt: ' + rel(path))
    return data


def make_report(cases, transfer, receipts, entries, historical, tests, subtests, diagnostics):
    lines = ['# Joint pure-coherence energy certificates and a transfer limitation', '',
             'ENERGYv19 and FAMILYv15 integrate the two independent pure-coherence constraints identified after complete diagonal spin-word closure. Reoptimizing older coefficients jointly with the new terms produces stronger exact energy lower bounds. Both new projector moments close, but full quantum overlap still requires a separate test. The frozen held-out comparison reveals a transfer limitation.', '',
             '## Accepted matched results', '',
             '| Case | Periodic lower/site | Million-site open lower/site | Enlarged family ceiling | Family gap |',
             '|---|---:|---:|---:|---:|---:|']
    for case, data in cases.items():
        lines.append(f"| {case} | {float(F(data['periodic_lower'])):.17g} | {float(F(data['open_lower_per_site'])):.17g} | {float(F(data['family_ceiling'])):.17g} | {data['family_gap_float']:.12g} |")
    lines += ['', 'The family ceilings limit attainable lower certificates in this specified relaxation. They are not physical ground-energy upper bounds. Physical uppers remain -0.6106763470511881 for W=0 and -0.6184244823693281 for W=1.', '']
    for case, data in cases.items():
        lines.append(f"- {case}: the lower improves by {float(F(data['lower_improvement'])):.12g}/site over the preceding certificate. Signed separation from the preceding full-family ceiling is {float(F(data['signed_previous_family_separation'])):.12g}; strict full-family separation proved: {data['strict_previous_family_separation']}. The selected enlarged ceiling uses {data['mixture_sources']} positive sources, maximum weight length {data['max_weight_characters']} characters. Selected evidence: `{data['directory']}`.")
    lines += ['', 'The preceding ENERGYv18 and FAMILYv14 proofs were freshly replayed under the current implementation. A stronger certificate alone does not prove the new directions outperform every reoptimized recipe in the preceding family.', '',
              '## What joint reoptimization establishes', '',
              'The earlier exact one-source and three-source caps allowed arbitrary real coefficients on the two new terms while keeping all older fields frozen. Those physical caps were reconstructed again here under current code. Both jointly optimized lower bounds exceed those caps:', '',
              '| Case | Joint lower minus frozen-old-recipe ceiling | Removing only new terms from current fixed recipe loses |',
              '|---|---:|---:|']
    for case, data in cases.items():
        lines.append(f"| {case} | {float(F(data['joint_excess_above_frozen_ceiling'])):.12g} | {float(F(data['matched_fixed_recipe_contribution'])):.12g} |")
    lines += ['', 'These comparisons have different scopes. Exceeding the frozen-old-recipe ceiling proves that allowing older coefficients to respond mattered. The matched ablation measures removal of the new terms from the resulting fixed recipe. Neither comparison establishes separation from the reoptimized preceding full family. Certificate byte identity binds the comparisons to the selected polished energy when a different mixture directory supplies the best ceiling.', '',
              '## Exact closure and residual quantum consistency', '',
              'The new sources are |346>+|409> and |314>+|614>. Removing the diagonal from their symmetry-averaged reflection-odd projector telescopes leaves five-site operators with 8 and 16 entries and six-site telescopes with 60 and 120 entries. Their norms are bounded by 1/2. The first carries four-bit occupation-dependent spin exchange, the second six-bit coherence. Prior exact annihilator functionals establish independence from all older affine terms.', '',
              'Each selected positive mixture has both new moments exactly zero, all 120 spin-word moments zero, and the entire preceding constraint hierarchy. An independent complete five-site density-matrix reconstruction also checks the two original positive-projector differences, rather than relying only on the production moment fields. Full diagonal prefix/suffix laws still agree and retain their stationary classical order-five extension.', '']
    for case, data in cases.items():
        witness = data['residual_overlap_witness']
        lines.append(f"- {case}: complete overlap has {data['residual_difference_diagonal_nonzeros']} diagonal mismatches and {data['residual_difference_upper_nonzeros']} nonzero upper-triangle differences. Quantum extension of this particular symmetry-averaged mixture refuted: {data['residual_full_overlap_refuted']}." + (f" Positive-projector vector {witness['five_site_vector']} detects mismatch {witness['difference_float']:.12g}." if witness else ''))
    lines += ['', 'Local positivity remains intact. A residual overlap refutes extension of the particular mixture, not all mixtures at its energy or the accepted ceiling. Classical occupation-law consistency does not extend the quantum coherences or impose fixed global particle number. The complete symmetry-compatible real test space remains 3960-dimensional: 120 diagonal and 3840 offdiagonal directions.', '',
              '## Frozen transfer', '',
              'The initial W=0 recipe before polishing was frozen and transferred to U=5, t=1, V=1/4, W=-1/5. All projector sources, penalties and auxiliary coefficients were preserved; the recorded physical profiles were rescaled or shifted and only the scalar threshold was recomputed. Geometry and filling were unchanged.', '',
              '| Million-site open-chain recipe | Lower/site |', '|---|---:|',
              f"| New recipe with pure-coherence terms | {float(F(transfer['with_pure_coherence_open_lower'])):.17g} |",
              f"| Only those two terms removed | {float(F(transfer['without_pure_coherence_open_lower'])):.17g} |",
              f"| Previous frozen spin-word recipe | {float(F(transfer['previous_frozen_open_lower'])):.17g} |", '',
              f"The signed contribution of the new terms is {transfer['ablation_loss_float']:.12g}/site. The signed change of the complete new recipe relative to the previous frozen recipe is {transfer['previous_recipe_improvement_float']:.12g}/site: a negative value means a worse lower bound. Thus helpful new terms do not establish robustness of the full jointly adapted recipe. All three energy replays and exact frozen-field comparisons pass. The physical upper remains {float(F(transfer['physical_upper'])):.17g}. No broader transfer claim is made.", '',
              '## Implementation and validation', '',
              'The new production module reuses the existing signed projector builder and removes diagonal entries, which commutes with its checked signed permutations and embeddings. Both new operators enter every local matrix before PSD replay. ENERGYv19 requires the full previous hierarchy; older versions reject the new field. FAMILYv15 requires both new zero moments and rejects fixed-coefficient fields. Its conservative source cap is 209, while the numerical system has 200 independent rows. The 4096-character weight limit, 4096-state coverage, 94 blocks and maximum local PSD dimension 200 are unchanged.', '',
              'The joint search varies 199 coefficients: 120 full spin-word directions, 75 offdiagonal directions, two hopping profiles and two penalties. Each run is bounded by 500 full-spectrum evaluations, with gradient/Hessian checks and fresh physical reconstruction. Initial W0/W1 searches used 460/500 evaluations; polishing used 461/500. Family pricing uses bounded physical integer vectors, then exact fraction-free reconstruction and independent physical replay. Initial pricing has 40 two-eigenvector rounds and refinement has up to 80 one-eigenvector rounds. Iteration limits and remaining negative reduced eigenvalues do not prove convergence. Timing records are not controlled performance benchmarks.', '']
    for name, diagnostic in diagnostics.items():
        if diagnostic.get('exact_basis_failure'):
            lines.append(f"- Rejected exact basis in `{name}`: {diagnostic['exact_basis_failure']}. This attempt contributes no accepted ceiling.")
    for name, diagnostic in diagnostics.items():
        if diagnostic.get('constraint_scale') and diagnostic.get('proposal_written'):
            lines.append(f"- Zero-pricing export retry in `{name}` scaled the numerical constraint rows by {diagnostic['constraint_scale']} and used {diagnostic['export_method']}. It produced an exact proposal; independent physical replay and the original 4096-character weight gate decide acceptance. No additional pricing or gate relaxation was used.")
    lines += ['', f"All 74 focused integration tests, 3 fraction-free tests, and the full {tests} tests plus {subtests} subtests pass. The existing calibration return-value warning remains. A copied fraction-free test initially retained the old dimension boundary; it was corrected to exercise 200 accepted rows and 201 refused rows, then passed. The failed test log is preserved. No production source changed after full-suite collection. The collector audits {receipts} current receipts and {entries} proof hash entries, plus {historical} historical receipts against unchanged sources or preserved snapshots. Both previous provenance manifests are accounted for. No agents, GPU or paid resources were used.", '',
              'Full-family numerical attainment, general quantum representability, generic molecular/long-range/higher-dimensional transfer, and scalability at requested accuracy remain unproved. The goal remains active.']
    return lines


def main():
    log = (BASE / 'full_validation.log').read_text()
    count = re.search(r'(\d+) passed', log)
    subtests = re.search(r'(\d+) subtests passed', log)
    suites = list(ET.parse(BASE / 'full_validation.xml').getroot().iter('testsuite'))
    if not count or not subtests or ' failed' in log or not suites or any(
        int(s.attrib.get(k, 0)) for s in suites for k in ('errors', 'failures', 'skipped')
    ):
        raise ValueError('Completed passing full suite required')
    for filename, expected in [('focused_validation.log', '74 passed'), ('fraction_free_validation.log', '3 passed')]:
        if expected not in (BASE / filename).read_text():
            raise ValueError('Focused gate missing: ' + filename)
    current, files = [], set()
    entries = construction = 0
    for path in BASE.rglob('*.json'):
        if path.name in ('summary.json', 'provenance.json'):
            continue
        data = read(path)
        if not isinstance(data, dict):
            continue
        for name, digest in data.get('source_sha256', {}).items():
            source = ROOT / name
            if sha(source) != digest:
                raise ValueError('Stale current proof/construction source: ' + name)
            files.add(source)
            construction += 1
        if data.get('accepted') is True and 'source_sha256' in data:
            current.append(path)
            entries += len(data['source_sha256'])
    cases = {}
    for case in ('W_zero', 'W_plus_1'):
        initial = BASE / case / 'final'
        choices = [initial]
        for name in ('refined', 'scaled'):
            candidate = BASE / case / name
            receipt = candidate / 'range_two_family_limit_replay.json'
            if receipt.exists() and read(receipt).get('accepted') is True:
                choices.append(candidate)
        folder = min(choices, key=lambda p: F(read(p / 'range_two_family_limit_replay.json')['periodic_family_upper']))
        energy = accepted(folder / 'range_two_replay.json')
        family = accepted(folder / 'range_two_family_limit_replay.json')
        closure = accepted(folder / 'pure_coherence_closure.json')
        overlap = accepted(folder / 'overlap/full_overlap_replay.json')
        comparison = accepted(initial / 'previous_family_comparison.json')
        ablation = accepted(initial / 'ablation_comparison.json')
        initial_family = accepted(initial / 'range_two_family_limit_replay.json')
        if sha(folder / 'profile_joint_r1_2_certificate.json') != sha(initial / 'profile_joint_r1_2_certificate.json'):
            raise ValueError('Comparisons do not bind selected energy')
        if not family.get('pure_coherence') or not family.get('full_spin_word') or not closure.get('full_spin_word_overlap_exactly_zero') or not closure.get('classical_markov_extension'):
            raise ValueError('Full diagonal closure missing')
        if overlap.get('difference_diagonal_nonzeros', 0):
            raise ValueError('Independent full-RDM diagonal closure disagrees')
        if overlap.get('pure_coherence_projector_moments') != {'346,409,1':'0','314,614,1':'0'}:
            raise ValueError('Independent complete RDM coherence closure missing')
        frozen_limit = accepted(initial / 'frozen_limit_comparison.json')
        nearest = family['family_replay']['nearest_constraint_replay']
        for key, size in [('pure_coherence_moments', 2), ('spin_word_moments', 120), ('coherent_projector_moments', 2),
                          ('three_spectator_hopping_moments', 18), ('two_spectator_hopping_moments', 30),
                          ('pair_transfer_moments', 4), ('spin_telescope_moments', 4), ('spectator_hopping_moments', 14)]:
            if len(nearest[key]) != size or any(map(F, nearest[key].values())):
                raise ValueError('Incomplete exact hierarchy: ' + key)
        lower = F(energy['lower_replay']['periodic_lower_density'])
        cap = F(family['periodic_family_upper'])
        initial_cap = F(initial_family['periodic_family_upper'])
        gap = cap - lower
        if gap < 0 or gap != F(family['family_gap']) or cap > initial_cap:
            raise ValueError('Inconsistent or worsened selected family bracket')
        certificate = read(folder / 'range_two_family_limit_certificate.json')
        max_chars = max(len(item['weight']) for item in certificate['mixture'])
        if max_chars > 4096:
            raise ValueError('Existing rational weight gate exceeded')
        cases[case] = {
            'directory': rel(folder), 'periodic_lower': str(lower), 'family_ceiling': str(cap),
            'family_gap': str(gap), 'family_gap_float': float(gap),
            'initial_family_ceiling': str(initial_cap), 'family_refinement': str(initial_cap - cap),
            'open_lower_per_site': energy['lower_per_site'], 'physical_upper_per_site': energy['upper_per_site'],
            'mixture_sources': nearest['mixture_sources'], 'max_weight_characters': max_chars,
            'lower_improvement': comparison['lower_improvement'],
            'strict_previous_family_separation': comparison['strict_family_separation'],
            'signed_previous_family_separation': comparison['exact_separation'],
            'matched_fixed_recipe_contribution': ablation['exact_fixed_recipe_contribution'],
            'joint_excess_above_frozen_ceiling': frozen_limit['joint_excess_above_frozen_ceiling'],
            'strict_joint_improvement_beyond_frozen_recipe': frozen_limit['strict_joint_improvement_beyond_frozen_recipe'],
            'full_diagonal_overlap_zero': True, 'classical_markov_extension': True,
            'markov_states_with_positive_mass': closure['markov_states_with_positive_mass'],
            'markov_edges_with_positive_flow': closure['markov_edges_with_positive_flow'],
            'residual_full_overlap_refuted': overlap.get('stationary_extension_refuted', False),
            'residual_overlap_witness': overlap.get('witness'),
            'residual_difference_upper_nonzeros': overlap.get('difference_upper_nonzeros', 0),
            'residual_difference_diagonal_nonzeros': overlap.get('difference_diagonal_nonzeros', 0),
        }
    transfer = accepted(BASE / 'held_out/frozen_transfer_comparison.json')
    space = accepted(BASE / 'overlap_space.json')
    old_root = BASE.parent / 'spin_word'
    obstruction_root = BASE.parent / 'spin_coherence'
    old = accepted(old_root / 'summary.json')
    historical = old['current_receipts'] + old['historical_audit']['receipt_paths'] + accepted(obstruction_root / 'summary.json')['current_receipts']
    if len(historical) != len(set(historical)):
        raise ValueError('Duplicate historical receipt')
    snapshots = {}
    for name in ('pair_transfer/source_before', 'pair_transfer/family_source_before',
                 'two_spectator/source_before', 'three_spectator/source_before',
                 'coherent_projector/source_before', 'spin_word/source_before', 'pure_coherence/source_before'):
        directory = BASE.parent / name
        for source, digest in read(directory / 'original_sha256.json').items():
            path = directory / Path(source).name
            if sha(path) != digest:
                raise ValueError('Changed preserved source snapshot')
            snapshots[source, digest] = path
            files.add(path)
    unchanged = preserved = 0
    for name in historical:
        receipt = accepted(ROOT / name)
        for source, digest in receipt['source_sha256'].items():
            if sha(ROOT / source) == digest:
                unchanged += 1
            elif (source, digest) in snapshots:
                preserved += 1
            else:
                raise ValueError('Unaccounted historical source change: ' + source)
    manifest_audits = {}
    for previous_root in (old_root, obstruction_root):
        same = archived = 0
        for name, digest in read(previous_root / 'provenance.json')['sha256'].items():
            if sha(ROOT / name) == digest:
                same += 1
            elif (name, digest) in snapshots:
                archived += 1
            else:
                raise ValueError('Prior manifest changed without snapshot: ' + name)
        manifest_audits[rel(previous_root / 'provenance.json')] = {'unchanged': same, 'preserved_source_snapshots': archived}
    diagnostics = {}
    for path in BASE.rglob('*.json'):
        if path.name not in ('thermal_proposal.json', 'diagonal_family_limit_diagnostic.json'):
            continue
        data = read(path)
        diagnostics[rel(path)] = {k: data[k] for k in (
            'accepted', 'proposal_written', 'matrix_evaluations', 'optimizer_success', 'optimizer_message',
            'seconds', 'states', 'mixture_sources', 'basis_size', 'numerical_family_upper',
            'constraint_scale', 'export_method', 'rational_family_upper_float',
            'proposed_periodic_lower_float', 'exact_basis_failure') if k in data}
        if 'pricing' in data:
            diagnostics[rel(path)].update(pricing_rounds=len(data['pricing']), last_pricing=data['pricing'][-1:])
    result = {
        'accepted': True, 'energy_version': 19, 'family_version': 15, 'source_cap': 209,
        'independent_numerical_rows': 200, 'goal_status': 'active', 'cases': cases,
        'frozen_transfer': transfer, 'full_overlap_operator_space': space, 'diagnostics': diagnostics,
        'current_receipts': list(map(rel, current)), 'current_receipts_audited': len(current),
        'current_source_hash_entries_checked': entries, 'construction_source_hash_entries_checked': construction,
        'historical_audit': {'receipts': len(historical), 'receipt_paths': historical,
                             'unchanged_hash_entries': unchanged, 'preserved_snapshot_entries': preserved,
                             'scope': 'Historical source-version audit, not fresh current-code replay.'},
        'prior_manifest_audits': manifest_audits,
        'focused_tests_passed': 74, 'fraction_free_tests_passed': 3,
        'full_validation': {'tests_passed': int(count[1]), 'subtests_passed': int(subtests[1]),
                            'junit_suites': [s.attrib for s in suites]},
        'report': 'research/marginal_pure_coherence_energy.md', 'provenance': rel(BASE / 'provenance.json'),
        'remaining': [
            'Both numerical family gaps remain nonzero; no attainment or convergence proof.',
            'Fixed-recipe ablation does not prove separation from the reoptimized previous family.',
            'A classical occupation-law extension does not extend quantum coherences or impose fixed global particle number.',
            'The complete frozen new recipe is worse at the held-out target despite helpful new-term ablation.',
            'Full quantum representability, generic molecular/long-range/higher-dimensional transfer and requested-accuracy scalability remain unproved.',
        ],
    }
    lines = make_report(cases, transfer, len(current), entries, len(historical), count[1], subtests[1], diagnostics)
    report = ROOT / result['report']
    report.write_text('\n'.join(lines) + '\n')
    files.update(path for path in BASE.rglob('*') if path.is_file() and path.name not in ('summary.json', 'provenance.json'))
    files.update((BASE.parent / 'discovery').glob('pure_coherence*.py'))
    files.update((ROOT / 'tests').glob('test_marginal_pure_coherence*.py'))
    files.update([Path(__file__).resolve(), report, old_root / 'provenance.json', old_root / 'summary.json', obstruction_root / 'provenance.json', obstruction_root / 'summary.json'])
    files.update(ROOT / name for name in historical)
    write(BASE / 'provenance.json', {'sha256': {rel(path): sha(path) for path in sorted(files)},
                                   'scope': 'Current proof/construction sources, tests, report, numerical attempts and historical snapshots. Numerical proposals remain nonaccepting.'})
    write(BASE / 'summary.json', result)
    central = ROOT / 'results/marginal_final_validation.json'
    data = read(central)
    data['pure_coherence_energy'] = result
    data['latest_validation_scope'] = 'ENERGYv19/FAMILYv15 add two pure-coherence constraints. Jointly optimized exact lower bounds and enlarged family ceilings accepted. Full diagonal and new projector moments close; complete quantum overlap independently checked. Matched/frozen ablations and signed previous-recipe transfer comparison accepted. Full regression passes. General representability, full-family attainment and scalability remain unproved; goal active.'
    write(central, data)
    print(json.dumps({'accepted': True, 'current_receipts': len(current), 'current_hash_entries': entries,
                      'historical_receipts': len(historical), 'provenance_files': len(files),
                      'full_tests': int(count[1]), 'subtests': int(subtests[1])}))


if __name__ == '__main__':
    main()
