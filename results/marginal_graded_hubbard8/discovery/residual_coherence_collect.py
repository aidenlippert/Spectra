"""Collect compact coherence obstructions and exact fixed-recipe limitations."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/residual_coherence'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path):
    return str(path.relative_to(ROOT))


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def main():
    prior_root = BASE.parent / 'pure_coherence'
    prior = read(prior_root / 'summary.json')
    manifest = read(prior_root / 'provenance.json')['sha256']
    for name, expected in manifest.items():
        if sha(ROOT / name) != expected:
            raise ValueError('Previous proof or provenance changed: ' + name)
    validations = {}
    for stem, count in [('validation', 40)]:
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
    if not support['accepted'] or support['prior_telescope_count'] != 75 or F(support['new_direction_minor_determinant']) != -1:
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
              'focused_tests_passed': 40, 'focused_validation': validations,
              'production_sources_changed': False, 'full_suite_rerun': False,
              'prior_full_suite': prior['full_validation'],
              'report': 'research/marginal_residual_coherence_obstruction.md', 'provenance': rel(BASE / 'provenance.json'),
              'remaining': ['Joint reoptimization of older coefficients with new constraints has not been implemented or certified.',
                            'The full ENERGYv19 family gaps remain unchanged; the new ceiling applies only when all older fields are frozen.',
                            'No new energy lower certificate or coupling/geometry transfer is claimed.',
                            'General quantum representability, full-family attainment and requested-accuracy scalability remain unproved.']}
    lines = ['# Residual coherence directions and exact frozen-recipe limits', '',
        'Two remaining positive-projector overlap witnesses yield independent stationary constraints beyond ENERGYv19/FAMILYv15. Exact three-source local mixtures limit the benefit of changing only these two coefficients to less than 2e-7 per site for W=0 and 8e-8 for W=1. These are bounds on a frozen-recipe subproblem; the full family gaps remain open.', '',
        '## Independent constraints', '',
        'The sources are |358>+|409> and |103>+|358>. Symmetry averaging and taking the reflection-odd part produce Y; T is the difference of its left and right embeddings. Removing the diagonal commutes with these operations. The full spin-word closure of the selected mixtures makes the removed moment exactly zero. Direct physical contractions reproduce the independently accepted full-density-matrix mismatches.', '',
        '| Source case | Five-site entries | Six-site entries | Changed fermion bits | Overlap mismatch | Bound on T and open boundary sum |',
        '|---|---:|---:|---:|---:|---:|']
    for case, data in cases.items():
        t = data['telescope']
        lines.append(f"| {case} | {t['five_site_nonzeros']} | {t['six_site_nonzeros']} | {t['changed_fermion_bit_counts']} | {t['moment_float']:.12g} | {t['six_site_and_open_boundary_norm_bound']} |")
    lines += ['', 'The first source changes eight fermionic occupation bits. The second is a same-spin hop across four sites whose amplitude depends on the occupations of the other modes. It changes only two bits; a generic argument excluding all one-body support by bit count would be invalid. The physical Hamiltonian and its variable hopping profiles contain nearest-neighbor hopping, which these selected long-range entries exclude.', '',
        'Exact row elimination finds particularly compact separating functionals:', '',
        '- L0(M) = -16 M[413,1433].',
        '- L1(M) = 8 M[1382,1433].', '',
        'Both annihilate each of the 75 preceding offdiagonal operators exactly, including the two pure-coherence terms integrated in v19. Offdiagonal entries exclude all diagonal terms. Both lie in spin sector (4,2), outside the actual fixed half and charged projector supports, and neither is a nearest-neighbor transition. Their matrix on the two new operators is [[0,1],[1,0]], with determinant -1. Thus neither constraint is redundant modulo the full preceding affine family.', '',
        'All Hermiticity, spin conservation, required signed symmetries, fermionic embedding translations and the six-site cyclic cancellation are checked exactly. Maximum absolute row sums bound Y by 1/4 and 1/8 respectively; the corresponding telescopes and open translated sums are bounded by 1/2 and 1/4. A positive-projector overlap mismatch violates stationarity of the particular local mixture. It does not violate that mixture’s local positivity or refute every mixture at its energy.', '',
        '## Certified limits with the previous recipe frozen', '',
        'For the complete fixed local matrix A, a positive trace-one mixture rho satisfying Tr(rho T0)=Tr(rho T1)=0 proves', '',
        '`lambda_min(A + gamma0 T0 + gamma1 T1) <= Tr(rho A)`', '',
        'for arbitrary real coefficients gamma0 and gamma1. Every older coefficient, physical profile and projector penalty remains fixed. Subtract the same penalty offset and divide by five. No fidelity inequality is required because the penalties are fixed. The existing accepted recipe at gamma=0 supplies the lower endpoint of this subproblem’s rigorous bracket.', '',
        '| Case | Accepted periodic lower/site | Exact frozen-recipe ceiling | Maximum possible gain/site | Positive sources |',
        '|---|---:|---:|---:|---:|']
    for case, data in cases.items():
        cap = data['fixed_recipe_cap']
        lines.append(f"| {case} | {float(F(cap['accepted_seed_lower'])):.17g} | {cap['ceiling_float']:.17g} | {cap['maximum_gain_over_seed_float']:.12g} | {cap['mixture_sources']} |")
    lines += ['', 'The three-source witnesses have exact positive rational weights, exact trace one and both new moments zero. Their energies are recomputed from physical CAR actions, all old diagonal and offdiagonal corrections, and both fixed projectors using standard-library rational arithmetic. This includes the preceding pure-coherence coefficients; a separate nonzero offdiagonal test checks their contribution.', '',
        'These ceilings do not bound joint reoptimization of older fields. The full ENERGYv19/FAMILYv15 gaps remain 4.49049126240775e-5 and 1.6702121170760045e-4 per site. Independence and a violated overlap alone therefore do not justify expecting useful improvement from a search that freezes the old recipe.', '',
        '## Numerical probes and validation', '',
        '| Case | Proposed new coefficients | Unaccepted proposed gain/site | Spectral evaluations |',
        '|---|---|---:|---:|']
    for case, data in cases.items():
        n=data['numerical_probe']
        lines.append(f"| {case} | {n['new_coefficients']} | {n['proposed_gain_over_accepted_seed_float']:.12g} | {n['matrix_evaluations']} |")
    lines += ['', 'These numerical energy proposals have no new production schema or exact PSD acceptance. The spectral probes check the active derivative and fresh physical matrix reconstruction and are capped at 250 optimization evaluations. Separate three-row LP searches use at most 20 pricing rounds and exact basis reconstruction. Independent physical replay, not floating-point convergence, establishes the ceilings.', '',
        f"All 40 focused tests passed, including 22 new tests and the 18 preceding coherence/cap tests. Coverage includes independent partial traces, exact diagonal-removal accounting, norm checks, dependence and injected-old-term refusal, physical nearest-neighbor exclusion, frozen old-coherence contributions, amplitude rescaling, malformed mixtures, nonzero moments and wrong version/scope refusals. One new test initially supplied profiles with invalid declared means and reflection; the existing validation correctly refused them. The fixture was corrected to valid nonuniform profiles and the failed version/log retained. No production source changed; the preceding full result remains 1128 tests plus 102 subtests and was not rerun for this diagnostic-only change.", '',
        f"The collector checks five current accepted receipts, {entries} proof source-hash entries, and all {len(manifest)} prior provenance files unchanged. It binds the current operators and ceilings to the accepted v19 energy and v15 mixture sources. Numerical proposals remain explicitly nonaccepting. No agents, GPU or paid resources were used.", '',
        'No energy certificate or held-out transfer result changes in this turn. The preceding whole frozen recipe remains worse at its held-out target despite helpful individual terms. Joint reoptimization, other coherence directions, full-family attainment, general quantum representability, broader transfer and requested-accuracy scalability remain unresolved. The goal stays active.']
    report = ROOT / result['report']
    report.write_text('\n'.join(lines)+'\n')
    files.update(path for path in BASE.rglob('*') if path.is_file() and path.name not in ('summary.json','provenance.json'))
    files.update((BASE.parent/'discovery').glob('residual_coherence*.py'))
    files.update((ROOT/'tests').glob('test_marginal_residual_coherence*.py'))
    files.update((ROOT/'tests').glob('test_marginal_spin_coherence*.py'))
    files.update([Path(__file__).resolve(), report, prior_root/'summary.json', prior_root/'provenance.json', BASE.parent/'residual_coherence_obstruction.log'])
    write(BASE/'provenance.json', {'sha256':{rel(p):sha(p) for p in sorted(files)},
        'scope':'Five exact receipts, numerical proposals, source versions, tests and report. Previous607-file provenance unchanged.'})
    write(BASE/'summary.json', result)
    central=ROOT/'results/marginal_final_validation.json'
    data=read(central)
    data['residual_coherence_obstruction']=result
    data['latest_validation_scope']='Two independent residual coherence directions beyond ENERGYv19/FAMILYv15. Exact three-source ceilings cap gains from only these two coefficients at1.993e-7/7.147e-8 per site with every old field frozen.40 focused tests pass. No production or energy changes. Full-family attainment, representability and scalability remain unproved; goal active.'
    write(central,data)
    print(json.dumps({'accepted':True,'current_receipts':len(current),'source_hash_entries':entries,'construction_hash_entries':construction,'provenance_files':len(files),'prior_provenance_files':len(manifest),'focused_tests':40}))


if __name__ == '__main__':
    main()
