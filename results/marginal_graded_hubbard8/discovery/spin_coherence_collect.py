"""Collect compact coherence obstructions and exact fixed-recipe limitations."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/spin_coherence'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path):
    return str(path.relative_to(ROOT))


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def main():
    prior_root = BASE.parent / 'spin_word'
    prior = read(prior_root / 'summary.json')
    manifest = read(prior_root / 'provenance.json')['sha256']
    for name, expected in manifest.items():
        if sha(ROOT / name) != expected:
            raise ValueError('Previous proof or provenance changed: ' + name)
    validations = {}
    for stem, count in [('validation', 30), ('fixed_limit_validation', 9)]:
        suites = list(ET.parse(BASE / (stem + '.xml')).getroot().iter('testsuite'))
        if len(suites) != 1 or int(suites[0].attrib['tests']) != count or any(int(suites[0].attrib.get(k, 0)) for k in ('errors', 'failures', 'skipped')):
            raise ValueError('Focused test gate failed')
        if f'{count} passed' not in (BASE / (stem + '.log')).read_text():
            raise ValueError('Completed test log required')
        validations[stem] = suites[0].attrib
    files, current = set(), []
    entries = construction = 0
    for path in BASE.rglob('*.json'):
        if path.name in ('summary.json', 'provenance.json'):
            continue
        data = read(path)
        if not isinstance(data, dict):
            continue
        for name, expected in data.get('source_sha256', {}).items():
            source = ROOT / name
            if sha(source) != expected:
                raise ValueError('Stale current source: ' + name)
            files.add(source)
            construction += 1
        if data.get('accepted') is True:
            current.append(path)
            entries += len(data['source_sha256'])
    if len(current) != 5:
        raise ValueError('Exactly five current proof receipts expected')
    support = read(BASE / 'support_replay.json')
    if not support['accepted'] or support['prior_telescope_count'] != 73 or F(support['new_direction_minor_determinant']) != F(1, 128):
        raise ValueError('Two independent new directions not proved')
    cases = {}
    for case in ('W_zero', 'W_plus_1'):
        telescope = read(BASE / case / 'telescope_replay.json')
        cap = read(BASE / case / 'fixed_limit/fixed_family_replay.json')
        numeric = read(BASE / case / 'numeric/proposal.json')
        diagnostic = read(BASE / case / 'fixed_limit/family_diagnostic.json')
        if not telescope['accepted'] or not cap['accepted'] or numeric['accepted'] or diagnostic['accepted']:
            raise ValueError('Proof/proposal separation failed')
        if not telescope['zero_diagonal'] or F(telescope['removed_diagonal_moment']) or not F(telescope['exact_moment']):
            raise ValueError('Pure coherence obstruction failed')
        if cap['new_coherence_moments'] != ['0', '0'] or not cap['all_old_fields_frozen'] or not cap['new_coefficients_unrestricted_real']:
            raise ValueError('Exact fixed-recipe scope failed')
        lower, ceiling = F(cap['accepted_seed_lower']), F(cap['fixed_recipe_family_ceiling'])
        if lower != F(prior['cases'][case]['periodic_lower']) or ceiling-lower != F(cap['maximum_gain_over_seed']) or ceiling < lower:
            raise ValueError('Fixed-recipe bracket differs from accepted seed')
        numerical_discrepancy = abs(float(ceiling)-diagnostic['numerical_upper'])
        if numerical_discrepancy > 1e-9:
            raise ValueError('Independent physical ceiling disagrees with numerical chart')
        cases[case] = {'telescope': telescope, 'fixed_recipe_cap': cap, 'numerical_probe': numeric,
                       'family_discovery': diagnostic, 'numeric_exact_ceiling_discrepancy': numerical_discrepancy,
                       'unchanged_full_family_gap': prior['cases'][case]['family_gap']}
    result = {'accepted': True, 'goal_status': 'active', 'new_energy_certificate': False,
              'independent_new_directions': 2, 'cases': cases, 'support': support,
              'current_receipts': list(map(rel, current)), 'current_receipts_audited': len(current),
              'current_source_hash_entries_checked': entries, 'construction_source_hash_entries_checked': construction,
              'prior_provenance_files_checked': len(manifest), 'prior_summary_sha256': sha(prior_root / 'summary.json'),
              'prior_provenance_sha256': sha(prior_root / 'provenance.json'),
              'focused_tests_passed': 39, 'focused_validation': validations,
              'production_sources_changed': False, 'full_suite_rerun': False,
              'prior_full_suite': prior['full_validation'],
              'report': 'research/marginal_spin_coherence_obstruction.md', 'provenance': rel(BASE / 'provenance.json'),
              'remaining': ['Joint reoptimization of older coefficients with new constraints has not been implemented or certified.',
                            'The full ENERGYv18 family gaps remain unchanged; the new ceiling applies only when all older fields are frozen.',
                            'No new energy lower certificate or coupling/geometry transfer is claimed.',
                            'General quantum representability, full-family attainment and requested-accuracy scalability remain unproved.']}
    lines = ['# Compact coherence obstructions and a fixed-recipe limitation', '',
             'Two exact, independent stationary constraints remain after full diagonal spin-word closure. They can be made purely offdiagonal and sparse. However, exact local-mixture ceilings show that changing only these two coefficients, while freezing every older coefficient, cannot materially improve the current energy recipes.', '',
             '## Exact new constraints', '',
             'Start from the normalized positive projectors onto |346>+|409> and |314>+|614>. Average under particle-hole and spin flip, take the reflection-odd part Y, and form T=Y_left−Y_right. Remove Y’s diagonal: the complete spin-word constraints make its expectation zero in the selected symmetry-averaged mixtures. Direct contractions of the resulting pure-coherence T still equal the independently reconstructed positive-projector overlap mismatch.', '',
             '| Source case | Five-site entries | Six-site entries | Changed fermion bits | Exact mismatch, decimal approximation | Norm bound on T |',
             '|---|---:|---:|---:|---:|---:|']
    for case, data in cases.items():
        t = data['telescope']
        lines.append(f"| {case} | {t['five_site_nonzeros']} | {t['six_site_nonzeros']} | {t['changed_fermion_bit_counts']} | {t['moment_float']:.12g} | {t['six_site_and_open_boundary_norm_bound']} |")
    lines += ['', 'Each Y has maximum absolute row sum 1/4, giving spectral norm at most 1/4. Both T and any open translated telescoping sum therefore have norm at most 1/2. Exact signed-permutation checks establish Hermiticity, spin conservation, all required symmetries, fermionic left/right translation and zero periodic sum. The original witnesses are positive projectors; their nonzero overlap difference is a consistency violation, not negative local positivity. The pure-coherence differences themselves need not be positive.', '',
              'The four-bit constraint is an occupation-dependent spin exchange and shares support with an old spin term. Merely counting changed modes would not establish independence. The exact linear functionals used instead are:', '',
              '- L0(M)=M[1370,1433]−M[350,413].',
              '- L1(M)=M[1337,1637].', '',
              'Both functionals annihilate all 73 existing nondiagonal operators. All selected entries lie in spin sector (4,2), excluding the actual fixed half and charged projectors. Four or six changed bits exclude physical one-body hopping, and offdiagonality excludes every diagonal term, including all 120 spin-word directions. On the two new T operators, the functional matrix is diag(1/8,1/16), with determinant 1/128. This proves two independent directions modulo the entire ENERGYv18/FAMILYv14 affine family.', '',
              '## Exact limitation when older coefficients are frozen', '',
              'Let A be the complete accepted local matrix, including every frozen old coefficient and both fixed projector penalties. For a positive trace-one local mixture rho with Tr(rho T0)=Tr(rho T1)=0,', '',
              '`lambda_min(A + gamma0 T0 + gamma1 T1) <= Tr(rho A)`', '',
              'for every real pair gamma. Subtracting the unchanged projector penalty cost and dividing by five gives the following exact ceilings. The witnesses do not require projector fidelity inequalities because both penalties remain fixed. All source vectors, weights, moments and local energies are replayed with standard-library rational arithmetic and physical CAR actions.', '',
              '| Case | Accepted periodic seed lower/site | Fixed-old-recipe ceiling | Maximum possible gain | Positive sources |',
              '|---|---:|---:|---:|---:|---:|']
    for case, data in cases.items():
        cap = data['fixed_recipe_cap']
        lines.append(f"| {case} | {float(F(cap['accepted_seed_lower'])):.17g} | {cap['ceiling_float']:.17g} | {cap['maximum_gain_over_seed_float']:.12g} | {cap['mixture_sources']} |")
    lines += ['', 'These bounds cover unrestricted real new coefficients. They prove a narrow limitation of the fixed old recipes, not of joint reoptimization. The full previous-family gaps remain about 5.75e-5 and 3.80e-4 per site; neither full-family numerical limit is resolved here.', '',
              'The bounded two-variable spectral probes found zero coefficients and no gain for W=0. W=1 proposed coefficients −93/500000 and 179/500000, with a proposed gain of 1.4e-7/site. This remains an unaccepted numerical energy proposal: no new energy schema or exact PSD lower replay was introduced. Both probes check an active spectral derivative and fresh physical matrix reconstruction. The exact ceilings above independently bound what any two-coordinate search could achieve.', '',
              'A separate three-row LP search proposed the positive mixtures. It used at most 20 pricing rounds, bounded integer amplitudes and an exact fraction-free reconstruction before the independent physical replay. Numerical success is not acceptance. The initial replay implementation misread the existing action dictionary as a list and failed; focused tests caught the same error. That version and its failed logs are preserved. The corrected dictionary indexing passed the independent occupation-energy tests and both exact replays.', '',
              '## Validation and remaining scope', '',
              f"All 39 focused tests passed: 30 overlap/coherence tests and 9 fixed-recipe cap tests. They include independent partial traces, norm checks, diagonal-removal accounting, amplitude rescaling, direct vacuum occupation energies, nonzero-moment refusal, malformed mixtures and scope refusals. The collector audits five accepted receipts and {entries} proof source-hash entries; all {len(manifest)} prior provenance files are unchanged. No production source changed and no full suite was rerun. The preceding full result remains 1074 tests plus 102 subtests.", '',
              'No energy lower bound or held-out transfer result changes in this turn. The new constraints are ready for a joint reoptimization experiment, where older diagonal, offdiagonal, profile and penalty coefficients can respond. Other coherence directions may also matter. General quantum representability, full-family numerical attainment, broader molecular/long-range/higher-dimensional transfer, and scalability at requested accuracy remain unproved. No agents, GPU or paid resources were used; the goal stays active.']
    report = ROOT / result['report']
    report.write_text('\n'.join(lines) + '\n')
    files.update(path for path in BASE.rglob('*') if path.is_file() and path.name not in ('summary.json', 'provenance.json'))
    files.update((BASE.parent / 'discovery').glob('spin_coherence*.py'))
    files.update((ROOT / 'tests').glob('test_marginal_spin_coherence*.py'))
    files.update([Path(__file__).resolve(), report, prior_root / 'summary.json', prior_root / 'provenance.json', BASE.parent / 'spin_coherence_replay.log'])
    write(BASE / 'provenance.json', {'sha256': {rel(path): sha(path) for path in sorted(files)},
                                   'scope': 'Five exact receipts, numerical proposals, tests, report and preserved failed replay version. Prior provenance verified unchanged.'})
    write(BASE / 'summary.json', result)
    central = ROOT / 'results/marginal_final_validation.json'
    data = read(central)
    data['spin_coherence_obstruction'] = result
    data['latest_validation_scope'] = 'Two compact pure-coherence constraints are exactly independent beyond ENERGYv18/FAMILYv14. Exact one/three-source ceilings restrict improvement from varying only those two coefficients to3.77e-8/1.73e-7 per site. All older fields must remain frozen for these ceilings.39 focused tests pass; production and energy bounds unchanged. Full-family attainment and representability remain unproved; goal active.'
    write(central, data)
    print(json.dumps({'accepted': True, 'current_receipts': len(current), 'source_hash_entries': entries,
                      'provenance_files': len(files), 'prior_provenance_files': len(manifest), 'focused_tests': 39}))


if __name__ == '__main__':
    main()
