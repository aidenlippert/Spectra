"""Collect exact energy, family, ablation, transfer and residual-overlap evidence."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,re,xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/marginal_graded_hubbard8/coherent_projector'

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return str(p.relative_to(ROOT))
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n')

def main():
    log=(BASE/'full_validation.log').read_text();count=re.search(r'(\d+) passed',log);subtests=re.search(r'(\d+) subtests passed',log)
    suites=list(ET.parse(BASE/'full_validation.xml').getroot().iter('testsuite'))
    if not count or not subtests or ' failed' in log or not suites or any(int(s.attrib.get(k,0)) for s in suites for k in ('errors','failures','skipped')):
        raise ValueError('Completed passing full suite required')
    if '81 passed' not in (BASE/'focused_validation.log').read_text() or '3 passed' not in (BASE/'fraction_free_validation.log').read_text():raise ValueError('Focused gates missing')
    current=[];files=set();entries=construction=0
    for p in BASE.rglob('*.json'):
        if p.name in ('summary.json','provenance.json'):continue
        d=read(p)
        if not isinstance(d,dict):continue
        for name,h in d.get('source_sha256',{}).items():
            source=ROOT/name
            if sha(source)!=h:raise ValueError('Stale current construction/proof source: '+name)
            files.add(source);construction+=1
        if d.get('accepted') is True and 'source_sha256' in d:
            current.append(p);entries+=len(d['source_sha256'])
    cases={}
    for case in ('W_zero','W_plus_1'):
        folder=BASE/case/'final'
        names=('range_two_replay','range_two_family_limit_replay','coherent_projector_closure','previous_family_comparison','ablation_comparison')
        for name in names:
            if read(folder/(name+'.json')).get('accepted') is not True:raise ValueError('Missing final proof '+case+'/'+name)
        energy,family,closure,comparison,ablation=[read(folder/(name+'.json')) for name in names]
        overlap=read(folder/'overlap/full_overlap_replay.json')
        diagonal_overlap=read(folder/'overlap/diagonal_overlap_replay.json')
        if not diagonal_overlap.get('accepted') or not diagonal_overlap.get('coarse_charge_word_overlap_exactly_zero'):raise ValueError('Spin-word/charge-word comparison missing')
        if not family.get('coherent_projector') or not overlap.get('accepted'):raise ValueError('Enlarged family/overlap proof missing')
        if any(map(F,closure['coherent_projector_moments'].values())):raise ValueError('Independent coherent closure failed')
        lower=F(energy['lower_replay']['periodic_lower_density']);cap=F(family['periodic_family_upper']);gap=cap-lower
        if gap<0 or gap!=F(family['family_gap']):raise ValueError('Inconsistent family interval')
        nearest=family['family_replay']['nearest_constraint_replay']
        for key,n in [('coherent_projector_moments',2),('three_spectator_hopping_moments',18),('two_spectator_hopping_moments',30),('pair_transfer_moments',4),('spin_telescope_moments',4),('spectator_hopping_moments',14)]:
            if len(nearest[key])!=n or any(map(F,nearest[key].values())):raise ValueError('Incomplete exact closure: '+key)
        cases[case]={'directory':rel(folder),'periodic_lower':str(lower),'family_ceiling':str(cap),'family_gap':str(gap),'family_gap_float':float(gap),
                     'open_lower_per_site':energy['lower_per_site'],'physical_upper_per_site':energy['upper_per_site'],
                     'mixture_sources':nearest['mixture_sources'],'lower_improvement':comparison['lower_improvement'],
                     'strict_previous_family_separation':comparison['strict_family_separation'],'signed_previous_family_separation':comparison['exact_separation'],
                     'matched_fixed_recipe_contribution':ablation['exact_fixed_recipe_contribution'],
                     'residual_full_overlap_refuted':overlap.get('stationary_extension_refuted',False),
                     'residual_overlap_witness':overlap.get('witness'),
                     'residual_difference_upper_nonzeros':overlap.get('difference_upper_nonzeros',0),
                     'spin_word_diagonal_total_variation':diagonal_overlap['exact_diagonal_total_variation'],
                     'coarse_charge_word_overlap_zero':True,
                     'positive_diagonal_projector_states':len(diagonal_overlap['positive_projector_states'])}
    transfer=read(BASE/'held_out/frozen_transfer_comparison.json');space=read(BASE/'overlap_space.json')
    if not transfer.get('accepted') or not space.get('accepted'):raise ValueError('Transfer or space census missing')
    # Historical proofs remain proofs of their recorded source versions.
    old_root=BASE.parent/'three_spectator';old=read(old_root/'summary.json');overlap_root=BASE.parent/'full_overlap';old_overlap=read(overlap_root/'summary.json')
    historical=old['current_receipts']+old['historical_audit']['receipt_paths']+old_overlap['current_receipts']
    if len(historical)!=len(set(historical)):raise ValueError('Duplicate historical receipt')
    snapshots={}
    dirs=[BASE.parent/'pair_transfer/source_before',BASE.parent/'pair_transfer/family_source_before',BASE.parent/'two_spectator/source_before',old_root/'source_before',BASE/'source_before']
    for directory in dirs:
        for name,h in read(directory/'original_sha256.json').items():
            p=directory/Path(name).name
            if sha(p)!=h:raise ValueError('Changed preserved source snapshot')
            snapshots[name,h]=p;files.add(p)
    unchanged=preserved=0
    for name in historical:
        receipt=read(ROOT/name)
        if not receipt.get('accepted'):raise ValueError('Unaccepted historical receipt')
        for source,h in receipt['source_sha256'].items():
            if sha(ROOT/source)==h:unchanged+=1
            elif (source,h) in snapshots:preserved+=1
            else:raise ValueError('Unaccounted historical source change '+source)
    manifest_audits={}
    for root in (old_root,overlap_root):
        manifest=read(root/'provenance.json')['sha256'];same=archived=0
        for name,h in manifest.items():
            if sha(ROOT/name)==h:same+=1
            elif (name,h) in snapshots:archived+=1
            else:raise ValueError('Prior manifest changed without snapshot: '+name)
        manifest_audits[rel(root/'provenance.json')]={'unchanged':same,'preserved_source_snapshots':archived}
        files.update([root/'provenance.json',root/'summary.json'])
    diagnostics={}
    for p in BASE.rglob('*.json'):
        if p.name not in ('thermal_proposal.json','diagonal_family_limit_diagnostic.json'):continue
        d=read(p);diagnostics[rel(p)]={k:d[k] for k in ('accepted','proposal_written','matrix_evaluations','optimizer_success','optimizer_message','seconds','states','mixture_sources','rational_family_upper_float','proposed_periodic_lower_float','exact_basis_failure') if k in d}
        if 'pricing' in d:diagnostics[rel(p)].update(pricing_rounds=len(d['pricing']),last_pricing=d['pricing'][-1:])
    result={'accepted':True,'energy_version':17,'family_version':13,'source_cap':139,'goal_status':'active','cases':cases,
            'frozen_transfer':transfer,'full_overlap_operator_space':space,'diagnostics':diagnostics,
            'current_receipts':list(map(rel,current)),'current_receipts_audited':len(current),'current_source_hash_entries_checked':entries,
            'construction_source_hash_entries_checked':construction,
            'historical_audit':{'receipts':len(historical),'receipt_paths':historical,'unchanged_hash_entries':unchanged,'preserved_snapshot_entries':preserved,
                                'scope':'Historical source-version audit, not fresh current-code replay.'},
            'prior_manifest_audits':manifest_audits,'focused_tests_passed':81,'fraction_free_tests_passed':3,
            'full_validation':{'tests_passed':int(count[1]),'subtests_passed':int(subtests[1]),'junit_suites':[s.attrib for s in suites]},
            'report':'research/marginal_coherent_projector_energy.md','provenance':rel(BASE/'provenance.json'),
            'remaining':['Enlarged numerical limits and attainment remain unproved; both family gaps are nonzero.',
                         'Fixed-recipe matched ablation is not separation from the reoptimized preceding family.',
                         'Frozen new operators harm the held-out recipe; the jointly changed recipe improves on the preceding frozen recipe.',
                         'Both implemented overlap moments close, but complete five-site quantum overlap still requires separate examination.',
                         'General representability, generic molecular/long-range/higher-dimensional transfer and requested-accuracy scalability remain unproved.']}
    lines=['# Coherent-projector energy certificates and a transfer limitation','',
           'ENERGYv17 and FAMILYv13 integrate the two exact five-site projector-overlap constraints found in the preceding diagnostic. Both matched-model lower bounds improve after fresh exact PSD replay. The new terms help the matched fixed recipes, but their frozen contribution is harmful at the held-out couplings. Neither enlarged numerical limit is resolved.','',
           '## Accepted matched results','',
           '| Case | Periodic lower/site | Million-site open lower/site | Family ceiling | Remaining family gap |','|---|---:|---:|---:|---:|']
    for name,d in cases.items():lines.append(f"| {name} | {float(F(d['periodic_lower'])):.17g} | {float(F(d['open_lower_per_site'])):.17g} | {float(F(d['family_ceiling'])):.17g} | {d['family_gap_float']:.12g} |")
    lines+=['','Family ceilings bound attainable LOWER certificates in the specified family; they are not physical ground-energy uppers. Physical upper bounds remain -0.6106763470511881 (W=0) and -0.6184244823693281 (W=1). Each selected family mixture has 139 positive physical sources and satisfies both new moments and the preceding hierarchy exactly. Independent projector contractions also vanish.','']
    for name,d in cases.items():lines.append(f"- {name}: lower improvement {float(F(d['lower_improvement'])):.12g}/site; signed separation from the preceding whole-family ceiling {float(F(d['signed_previous_family_separation'])):.12g}; strict whole-family separation proved: {d['strict_previous_family_separation']}. Removing only the two new terms from this fixed recipe loses {float(F(d['matched_fixed_recipe_contribution'])):.12g}/site.")
    lines+=['','The larger matched ablation effects concern jointly adapted fixed coefficients. They do not imply the same improvement over a reoptimized older family. The preceding energy and family certificates were freshly replayed under the current implementation for these comparisons.','',
            '## Frozen transfer and signed ablation','',
            'The INITIAL W=0 recipe, before polishing, was frozen and transferred to U=5, t=1, V=1/4, W=-1/5. Projector sources, penalties and all auxiliary coefficients were preserved. Physical profiles were rescaled or shifted by the recorded recipe; only the scalar spectral threshold was recomputed. Geometry, filling and interaction range stayed fixed.','',
            '| Million-site open-chain recipe | Lower/site |','|---|---:|',
            f"| With both coherent-projector terms | {float(F(transfer['with_coherent_open_lower'])):.17g} |",
            f"| Only those two terms removed | {float(F(transfer['without_coherent_open_lower'])):.17g} |",
            f"| Previous frozen three-spectator recipe | {float(F(transfer['previous_frozen_open_lower'])):.17g} |",'',
            f"The signed contribution of the frozen new terms is **{transfer['ablation_loss_float']:.12g}/site**. Removing them improves the bound: their matched benefit does not transfer with these coefficients. The full new frozen recipe nevertheless improves on the previous frozen recipe by {transfer['previous_recipe_improvement_float']:.12g}/site. Its other jointly adapted coefficients matter. The physical upper is -0.4885616802989547. All three energy replays and exact frozen-field comparisons are accepted.",
            '', 'An initial use of the matched-target upper-state driver correctly refused all three held-out cases because that target has no matching stored filter recipe. Those refusal logs are retained. The dedicated transfer driver then evaluated the fixed W=0 physical filter at the actual held-out Hamiltonian and passed exact lower/upper replay; the source-selection gate was not weakened.','',
            '## What remains outside the implemented overlap constraints','',
            'An exact signed-permutation orbit census finds 32,264 real Hermitian five-site matrix directions preserving both spin numbers. Imposing particle-hole evenness, spin-flip evenness and reflection oddness leaves **3,960 directions: 120 diagonal and 3,840 offdiagonal**. Every matrix unit belongs to one checked orbit; disjoint consistent orbits give independent basis vectors and sign-conflicted orbits give none. This is an operator-space dimension statement, not a performance or general representability result.','']
    for name,d in cases.items():
        w=d['residual_overlap_witness']
        lines.append(f"- {name}: stationary extension of the new selected symmetry-averaged mixture refuted by complete overlap replay: {d['residual_full_overlap_refuted']}. Nonzero upper-triangle differences: {d['residual_difference_upper_nonzeros']}."+(f" Positive-projector witness {w['five_site_vector']} has exact mismatch summarized by {w['difference_float']:.12g}." if w else ''))
        lines.append(f"  The spin-resolved diagonal word laws have exact total variation {float(F(d['spin_word_diagonal_total_variation'])):.12g}, detected by a positive diagonal projector onto {d['positive_diagonal_projector_states']} states. Every coarser charge-word overlap is exactly zero. Thus even full diagonal spin-word consistency remains open, despite charge-law closure; the signed indicator telescope has norm at most one.")
    lines+=['','Closing the two added moments does not close full quantum consistency. A failed overlap test refutes the extension of that particular mixture, not all mixtures at its energy or the accepted family ceiling. The census motivates a broader operator-space treatment rather than assuming two more scalar tests establish representability.','',
            '## Implementation and validation','',
            'The production implementation fixes two canonical projector sources. It checks signed orbit closure, Hermiticity, spin numbers, reflection/particle-hole/spin-flip invariance and fermionic telescoping on the full 4096-state space. Both operators enter the local matrix before exact PSD acceptance. Local coverage remains 94 blocks, with maximum local PSD dimension 200. Older energy versions reject the new field. FAMILYv13 requires the full preceding hierarchy and two exact zero moments, allows at most 139 sources, and forbids fixed-coefficient fields. Older source caps and the 4096-character rational-weight limit are preserved.','',
            'Numerical proposals optimize 138 coefficients with at most 500 full-spectrum evaluations per run. Initial W0/W1 runs used 412/500 evaluations; polish runs used 411/500. Fresh physical reconstruction and derivative checks accompany the nonaccepting proposals. Initial family pricing ran 40 rounds with two eigenvectors per block; refinement ran 80 rounds with one. Exact bounded 139-row fraction-free reconstruction is separate from numerical selection. Remaining negative reduced eigenvalues and iteration/evaluation limits do not establish convergence.','',
            f"All **81 focused integration tests**, **3 focused fraction-free tests**, and the full **{count[1]} tests plus {subtests[1]} subtests** pass. The existing calibration return-value warning remains. No production source changed after full-suite collection. {len(current)} current accepted receipts and {entries} source-hash entries were checked; {len(historical)} historical receipts are audited against unchanged sources or preserved snapshots. Failed precondition attempts and untrusted proposals are retained. No GPU, paid resources or new agents were used.",
            '', 'General quantum representability, exact numerical-limit attainment, generic molecular or long-range/higher-dimensional transfer, and scalability at requested accuracy remain unproved. The goal stays active.']
    report=ROOT/result['report'];report.write_text('\n'.join(lines)+'\n')
    files.update(p for p in BASE.rglob('*') if p.is_file() and p.name not in ('summary.json','provenance.json'))
    files.update((BASE.parent/'discovery').glob('coherent_projector*.py'))
    files.update((ROOT/'tests').glob('test_marginal_coherent_projector*.py'))
    files.update([Path(__file__).resolve(),report]);files.update(ROOT/name for name in historical)
    write(BASE/'provenance.json',{'sha256':{rel(p):sha(p) for p in sorted(files)},'scope':'Current proof/construction sources, tests, report and historical source snapshots. Numerical proposals remain nonaccepting.'})
    write(BASE/'summary.json',result)
    central=ROOT/'results/marginal_final_validation.json';d=read(central);d['coherent_projector_energy']=result
    d['latest_validation_scope']='ENERGYv17/FAMILYv13 add two coherent-projector constraints. Stronger exact matched lower bounds, matching family ceilings and signed matched/held-out ablations accepted. New moments close but full quantum overlap remains separately tested. Full regression passes. General representability, numerical attainment and scalability remain unproved; goal active.'
    write(central,d)
    print(json.dumps({'accepted':True,'current_receipts':len(current),'current_hash_entries':entries,'historical_receipts':len(historical),'provenance_files':len(files),'full_tests':int(count[1]),'subtests':int(subtests[1])}))

if __name__=='__main__':main()
