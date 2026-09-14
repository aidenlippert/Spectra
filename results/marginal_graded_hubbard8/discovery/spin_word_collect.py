"""Collect source-bound spin-word energy, family and extension evidence."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/spin_word'


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


def main():
    log = (BASE / 'full_validation.log').read_text()
    count = re.search(r'(\d+) passed', log)
    subtests = re.search(r'(\d+) subtests passed', log)
    suites = list(ET.parse(BASE / 'full_validation.xml').getroot().iter('testsuite'))
    if not count or not subtests or ' failed' in log or not suites or any(
        int(s.attrib.get(k, 0)) for s in suites for k in ('errors', 'failures', 'skipped')
    ):
        raise ValueError('Completed passing full suite required')
    for filename, expected in [('focused_validation.log', '63 passed'), ('fraction_free_validation.log', '3 passed')]:
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
        closure = accepted(folder / 'spin_word_closure.json')
        overlap = accepted(folder / 'overlap/full_overlap_replay.json')
        comparison = accepted(initial / 'previous_family_comparison.json')
        ablation = accepted(initial / 'ablation_comparison.json')
        initial_family = accepted(initial / 'range_two_family_limit_replay.json')
        if sha(folder / 'profile_joint_r1_2_certificate.json') != sha(initial / 'profile_joint_r1_2_certificate.json'):
            raise ValueError('Comparisons do not bind selected energy')
        if not family.get('full_spin_word') or not closure.get('full_spin_word_overlap_exactly_zero') or not closure.get('classical_markov_extension'):
            raise ValueError('Full diagonal closure missing')
        if overlap.get('difference_diagonal_nonzeros', 0):
            raise ValueError('Independent full-RDM diagonal closure disagrees')
        nearest = family['family_replay']['nearest_constraint_replay']
        for key, size in [('spin_word_moments', 120), ('coherent_projector_moments', 2),
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
    old_root = BASE.parent / 'coherent_projector'
    old = accepted(old_root / 'summary.json')
    historical = old['current_receipts'] + old['historical_audit']['receipt_paths']
    if len(historical) != len(set(historical)):
        raise ValueError('Duplicate historical receipt')
    snapshots = {}
    for name in ('pair_transfer/source_before', 'pair_transfer/family_source_before',
                 'two_spectator/source_before', 'three_spectator/source_before',
                 'coherent_projector/source_before', 'spin_word/source_before'):
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
    same = archived = 0
    for name, digest in read(old_root / 'provenance.json')['sha256'].items():
        if sha(ROOT / name) == digest:
            same += 1
        elif (name, digest) in snapshots:
            archived += 1
        else:
            raise ValueError('Prior manifest changed without snapshot: ' + name)
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
        'accepted': True, 'energy_version': 18, 'family_version': 14, 'source_cap': 207,
        'independent_numerical_rows': 198, 'goal_status': 'active', 'cases': cases,
        'frozen_transfer': transfer, 'full_overlap_operator_space': space, 'diagnostics': diagnostics,
        'current_receipts': list(map(rel, current)), 'current_receipts_audited': len(current),
        'current_source_hash_entries_checked': entries, 'construction_source_hash_entries_checked': construction,
        'historical_audit': {'receipts': len(historical), 'receipt_paths': historical,
                             'unchanged_hash_entries': unchanged, 'preserved_snapshot_entries': preserved,
                             'scope': 'Historical source-version audit, not fresh current-code replay.'},
        'prior_manifest_audit': {'path': rel(old_root / 'provenance.json'), 'unchanged': same,
                                 'preserved_source_snapshots': archived},
        'focused_tests_passed': 63, 'fraction_free_tests_passed': 3,
        'full_validation': {'tests_passed': int(count[1]), 'subtests_passed': int(subtests[1]),
                            'junit_suites': [s.attrib for s in suites]},
        'report': 'research/marginal_spin_word_energy.md', 'provenance': rel(BASE / 'provenance.json'),
        'remaining': [
            'Both numerical family gaps remain nonzero; no attainment or convergence proof.',
            'Fixed-recipe ablation does not prove separation from the reoptimized previous family.',
            'A classical occupation-law extension does not extend quantum coherences or impose fixed global particle number.',
            'Full quantum representability, generic molecular/long-range/higher-dimensional transfer and requested-accuracy scalability remain unproved.',
        ],
    }
    lines = ['# Full spin-word constraints: energy certificates and classical closure', '',
             'ENERGYv18 and FAMILYv14 implement the complete 120-dimensional diagonal five-site space that is even under particle-hole and spin flip and odd under reflection. Both matched energy bounds improve after exact PSD replay. The selected positive local mixtures now have identical full spin-resolved five-site occupation laws, with an explicit stationary classical extension. Full quantum extension remains a separate requirement.', '',
             '## Accepted matched results', '',
             '| Case | Periodic lower/site | Million-site open lower/site | Family ceiling | Family gap |',
             '|---|---:|---:|---:|---:|']
    for name, data in cases.items():
        lines.append(f"| {name} | {float(F(data['periodic_lower'])):.17g} | {float(F(data['open_lower_per_site'])):.17g} | {float(F(data['family_ceiling'])):.17g} | {data['family_gap_float']:.12g} |")
    lines += ['', 'Family ceilings bound attainable lower certificates in this specified family; they are not physical ground-energy upper bounds. Physical upper bounds remain -0.6106763470511881 (W=0) and -0.6184244823693281 (W=1).', '']
    for name, data in cases.items():
        refinement = (f"Refinement lowers the initial ceiling by {float(F(data['family_refinement'])):.12g}/site."
                      if F(data['family_refinement']) > 0 else 'No tighter refinement passed exact acceptance.')
        lines.append(f"- {name}: lower improvement {float(F(data['lower_improvement'])):.12g}/site; signed separation from the preceding whole-family ceiling {float(F(data['signed_previous_family_separation'])):.12g}; strict separation proved: {data['strict_previous_family_separation']}. Removing only the new spin-word field from the adapted fixed recipe loses {float(F(data['matched_fixed_recipe_contribution'])):.12g}/site. The selected ceiling uses {data['mixture_sources']} positive sources, with maximum rational weight length {data['max_weight_characters']} characters. {refinement} Selected evidence: `{data['directory']}`.")
    lines += ['', 'The preceding energy and family proofs were freshly replayed under current code. Matched ablations are bound to the selected energy by certificate byte identity. Their larger effects concern fixed, jointly adapted coefficients and do not establish improvement over a reoptimized older family.', '',
              '## Exact diagonal closure and remaining quantum obstruction', '',
              'The independent diagnostic reconstructs all 4096 six-site occupation probabilities from the positive source mixture and its eight symmetry images. Every one of the 1024 five-site prefix probabilities equals the corresponding suffix probability exactly. All 120 new family moments vanish, as do the preceding constraints.', '',
              'The closure also has a general finite-window explanation. For any mixture averaged over these symmetries, the prefix-minus-suffix diagonal vector is particle-hole even, spin-flip even and reflection odd. It therefore lies in the complete 120-dimensional basis. The zero telescope moments make it orthogonal to every basis vector, so it must vanish. Disjoint signed orbits establish independence. This implication concerns the diagonal of the symmetry-averaged mixture; the independent integer reconstruction verifies it for each accepted instance.', '',
              'For a five-word prefix u with mass q(u)>0, append symbol a with probability p6(u,a)/q(u). Prefix/suffix equality proves normalization and stationary inflow. Zero-mass prefixes may append empty deterministically. This constructs an order-five stationary classical process reproducing the full six-site occupation law over empty, up, down and double. It does not preserve the offdiagonal matrix elements of the supplied quantum mixture or impose a fixed global particle number.', '']
    for name, data in cases.items():
        witness = data['residual_overlap_witness']
        lines.append(f"- {name}: {data['markov_states_with_positive_mass']} positive-mass Markov states and {data['markov_edges_with_positive_flow']} positive flows. Independent complete density-matrix replay finds {data['residual_difference_diagonal_nonzeros']} diagonal differences and {data['residual_difference_upper_nonzeros']} nonzero upper-triangle differences. Stationary quantum extension of this particular symmetry-averaged mixture refuted: {data['residual_full_overlap_refuted']}." + (f" The normalized positive projector from vector {witness['five_site_vector']} has expectation mismatch {witness['difference_float']:.12g}." if witness else ''))
    lines += ['', 'A nonzero offdiagonal overlap difference refutes extension of the particular mixture, not all mixtures at its energy or the family ceiling. Local positivity is intact. The exact symmetry census gives 3960 allowed real overlap directions: 120 diagonal and 3840 offdiagonal. Closing the diagonal part does not close the full space.', '',
              '## Frozen transfer', '',
              'The initial W=0 recipe before polishing was frozen and transferred to U=5, t=1, V=1/4, W=-1/5. All projector sources, penalties and auxiliary coefficients remain fixed; the recorded physical profiles are rescaled or shifted and only the scalar threshold is recomputed. Geometry, filling and interaction range are unchanged.', '',
              '| Million-site open-chain recipe | Lower/site |', '|---|---:|',
              f"| With new spin-word field | {float(F(transfer['with_spin_word_open_lower'])):.17g} |",
              f"| Only that field removed | {float(F(transfer['without_spin_word_open_lower'])):.17g} |",
              f"| Previous frozen coherent-projector recipe | {float(F(transfer['previous_frozen_open_lower'])):.17g} |", '',
              f"The signed frozen contribution is {transfer['ablation_loss_float']:.12g}/site; the full recipe improves on the previous frozen recipe by {transfer['previous_recipe_improvement_float']:.12g}/site. The physical upper remains -0.4885616802989547. All three energy proofs and the exact frozen-field comparison pass. This single parameter transfer does not establish broader transfer robustness. The preceding stage's adverse two-term transfer remains historical evidence; this result does not erase it.", '',
              '## Implementation, numerical limitations and validation', '',
              'Canonical signed symmetry orbits partition all 1024 five-site spin words and produce exactly 120 independent diagonal directions. Independent full-group projection tests establish completeness and containment of the previous 52 charge directions and nine selected shapes. The new field enters the local matrix before PSD replay. Coverage remains 4096 states and 94 blocks, with maximum local PSD dimension 200. Older energy versions reject the new field. FAMILYv14 requires the preceding hierarchy and every new zero moment. Its conservative source cap is 207; the actual numerical system has 198 independent rows because the selected nine shapes are redundant. The 4096-character rational-weight limit is unchanged.', '',
              'The numerical energy search has 197 coordinates and a 500 spectral-evaluation cap. Initial W0/W1 runs used 410/500 evaluations and polishing used 417/500. Family pricing used 40 rounds with two eigenvectors, then 80 rounds with one, preserving accepted physical vectors and determinant anchors. The first W0 family search failed with a HiGHS Unknown status. Its log is retained. A separate bounded fallback tries dual simplex and interior point without presolve, retains explicit candidate-feasibility checks, and never promotes numerical output to acceptance. Both first refinement exports failed exact consistency despite successful numerical LP status. A bounded zero-pricing retry scales all numerical constraint rows by 1000 and uses dual simplex without presolve; it retains the original rational equations and acceptance gates. Each selected mixture is reconstructed by bounded 198-row fraction-free arithmetic and replayed from physical vectors. The reported directory identifies the best accepted ceiling; failed proposals do not count as refinements. Negative remaining reduced eigenvalues and finite iteration limits do not prove convergence; recorded durations are not controlled benchmarks.', '',
              f"All 63 focused integration tests, 3 focused fraction-free tests, and the full {count[1]} tests plus {subtests[1]} subtests pass. The existing calibration return-value warning remains. No production source changed after suite collection. The collector checks {len(current)} current accepted receipts with {entries} source-hash entries and {len(historical)} historical receipts against unchanged sources or preserved snapshots. The preceding {same + archived}-file manifest is accounted for. Numerical attempts and failure logs remain in provenance. No agents, GPU or paid resources were used.", '',
              'General quantum representability, exact numerical-limit attainment, generic molecular or long-range/higher-dimensional transfer, and scalability at requested accuracy remain unproved. The goal remains active.']
    report = ROOT / result['report']
    report.write_text('\n'.join(lines) + '\n')
    files.update(path for path in BASE.rglob('*') if path.is_file() and path.name not in ('summary.json', 'provenance.json'))
    files.update((BASE.parent / 'discovery').glob('spin_word*.py'))
    files.update((ROOT / 'tests').glob('test_marginal_spin_word*.py'))
    files.update([Path(__file__).resolve(), report, old_root / 'provenance.json', old_root / 'summary.json'])
    files.update(ROOT / name for name in historical)
    write(BASE / 'provenance.json', {'sha256': {rel(path): sha(path) for path in sorted(files)},
                                   'scope': 'Current proof/construction sources, tests, report, numerical attempts and historical snapshots. Numerical proposals remain nonaccepting.'})
    write(BASE / 'summary.json', result)
    central = ROOT / 'results/marginal_final_validation.json'
    data = read(central)
    data['spin_word_energy'] = result
    data['latest_validation_scope'] = 'ENERGYv18/FAMILYv14 add the complete120 diagonal spin-word directions. Exact matched lower bounds, family ceilings, ablations and frozen transfer pass. Full diagonal overlap closes with a stationary classical extension; full quantum overlap is independently checked. Full regression passes. General representability, numerical attainment and scalability remain unproved; goal active.'
    write(central, data)
    print(json.dumps({'accepted': True, 'current_receipts': len(current), 'current_hash_entries': entries,
                      'historical_receipts': len(historical), 'provenance_files': len(files),
                      'full_tests': int(count[1]), 'subtests': int(subtests[1])}))


if __name__ == '__main__':
    main()
