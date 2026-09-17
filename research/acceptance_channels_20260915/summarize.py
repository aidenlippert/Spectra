"""Assemble the continuation's evidence without relabeling it a cold calculation."""
from datetime import datetime,timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from research.acceptance_channels_20260915.campaign import ROOT,OUT,dump
from research.molecular_collective_20260913.core import digest


def main():
    runs=[]
    for path in sorted((OUT/'runs').glob('*.json')):
        record=json.loads(path.read_text())
        if 'wall_seconds' in record:runs.append(record)
        elif record.get('status')=='starting':raise ValueError(('An experiment is still running',path.name))
    expected=digest(json.loads((ROOT/'results/interacting_scaling_20260915/models/h12_heldout/fixture.json').read_text()))
    exact={}
    for path in sorted((OUT/'replays').glob('*/original_interval.json')):
        record=json.loads(path.read_text());width=Fraction(record['upper_Ha'])-Fraction(record['lower_Ha'])
        if width!=Fraction(record['width_Ha']) or width<0 or record['original_fixture_sha256']!=expected:
            raise ValueError('Exact interval failed endpoint/model consistency')
        if record['fresh_complete_replay_in_this_invocation'] is not True:
            raise ValueError('A fresh full replay is required')
        lower=json.loads((path.parent/'exact/lower.json').read_text())
        exact[path.parent.name]={'width_mHa':float(1000*width),'target_met':width<=Fraction(1,625),
            'receipt':str(path.relative_to(ROOT)),'complete_replay_seconds':record['complete_seconds'],
            'singlet_residual_mHa':float(1000*Fraction(lower['singlet']['residual_l1'])),
            'upper_Ha':record['upper_Ha'],'lower_Ha':record['lower_Ha']}
    if not exact:raise ValueError('No new complete exact replay is available')
    proposals={}
    for case in ('h12','h12_cached'):
        for path in sorted((OUT/case).glob('*/discovery.json')):
            if path.parent.is_symlink():continue
            d=json.loads(path.read_text())
            width=d.get('predicted_width_mHa',d.get('best',{}).get('unverified_width_mHa'))
            proposals[f'{case}/{path.parent.name}']={'predicted_width_mHa':width,
                'exact_acceptance_required':True,'record':str(path.relative_to(ROOT))}
    old=json.loads((OUT/'baseline.json').read_text())
    iterations=json.loads((OUT/'iteration_comparison.json').read_text())
    residual_diagnosis=json.loads((OUT/'fixed_identity_residual_diagnosis.json').read_text())
    trace_path=OUT/'trace_diagnostic.json'
    trace=json.loads(trace_path.read_text()) if trace_path.exists() else None
    null_channels=json.loads((OUT/'singlet_null_channels.json').read_text())
    paired_path=OUT/'paired_replay_comparison.json'
    paired=json.loads(paired_path.read_text()) if paired_path.exists() else None
    original_component_path=OUT/'original_paired_component.json'
    component=None
    if original_component_path.exists():
        original_component=json.loads(original_component_path.read_text())
        contracted_component=json.loads((OUT/'paired_component.json').read_text())
        for key in ('source_certificate_sha256','selected_block_count','polynomial_terms','maximum_degree','polynomial_sha256','stats'):
            if original_component[key]!=contracted_component[key]:raise ValueError(('Paired component mismatch',key))
        contracted_seconds=contracted_component.get('evaluation_seconds',contracted_component.get('contraction_seconds'))
        component={'identical_exact_polynomial_and_factor_statistics':True,
            'original_evaluation_seconds':original_component['evaluation_seconds'],
            'contracted_evaluation_seconds':contracted_seconds,
            'ratio':original_component['evaluation_seconds']/contracted_seconds,
            'scope':'One isolated component evaluation each; this ratio is not a complete replay or solver speedup.'}
    best_name,best=min(exact.items(),key=lambda item:item[1]['width_mHa'])
    matched=None
    dense=exact.get('dense_t2',exact.get('dense_t2_long'))
    if dense and 'matched_continuation' in exact:
        a,b=exact['matched_continuation'],dense
        if Fraction(a['upper_Ha'])!=Fraction(b['upper_Ha']):raise ValueError('Matched upper endpoints differ')
        matched={'unchanged_width_mHa':a['width_mHa'],'four_direction_width_mHa':b['width_mHa'],
            'gain_mHa':a['width_mHa']-b['width_mHa'],
            'one_conditional_warm_comparison':True,'cold_speedup_claimed':False}
    result={'assembled_UTC':datetime.now(timezone.utc).isoformat(),'previous_exact_width_mHa':old['width_mHa'],
        'best_new_exact':{'name':best_name,**best},'all_exact_replays':exact,'numerical_proposals':proposals,
        'matched_comparison':matched,'completed_metered_stages':len(runs),
        'matched_iteration_diagnostic':iterations,'paired_replay_comparison':paired,
        'paired_component_comparison':component,
        'fixed_identity_residual_diagnosis':residual_diagnosis,
        'numerical_trace_diagnostic':None if trace is None else {k:v for k,v in trace.items() if k not in ('blocks','lift')},
        'exact_singlet_null_channel_receipt':'results/acceptance_channels_20260915/singlet_null_channels.json',
        'metered_stage_wall_seconds':sum(r['wall_seconds'] for r in runs),
        'peak_metered_child_RSS_bytes':max(r['peak_child_RSS_bytes'] for r in runs),
        'process_failures':[r for r in runs if r['status']!='passed'],
        'cost_scope':'New metered subprocess stages, including failed searches, retries, input caching and complete replay. Prior discovery/preparation and unmetered early unit-test/import time and lightweight receipt analyses are additional. Not a cold time-to-certificate.',
        'prior_cost_record':'results/interacting_scaling_20260915/accounting.json',
        'source_dual_from_attachment_replayed':False,'current_labelled_H8_channel_membership_proved':True,
        'external_orbital_basis_match_established':False,
        'exact_family_obstruction_proved':False,'new_transfer_test':False,
        'focused_regressions':{'tests':28,'passed':all(any(r['name']==name and r['status']=='passed' for r in runs)
            for name in ('final_focused_regressions','exact_paired_independent_oracle','exact_singlet_null_test')),
            'scope':'21 numerical/algebraic tests, six exact paired-contraction tests and one independent singlet-null test.'},
        'new_external_spending_USD':0,'github_operations':False,
        'run_receipt_sha256':{r['name']:hashlib.sha256((OUT/'runs'/(r['name']+'.json')).read_bytes()).hexdigest() for r in runs}}
    if not result['focused_regressions']['passed']:raise ValueError('Required focused tests did not pass')
    dump(OUT/'accounting.json',result)
    target='met' if best['target_met'] else 'not met'
    lines=['# H12 acceptance and collective-channel continuation','',
        f'The best new independently replayed original-model interval is **{best["width_mHa"]:.9f} mHa**, compared with the preserved **{old["width_mHa"]:.9f} mHa**. The 1.6 mHa target is **{target}**. No exact family-limit theorem was obtained.','',
        'The complete standard-library replay checks the actual MPS upper, singlet and nonsinglet lower witnesses, exact orbital transformation, and its charged coefficient allowance. The Hamiltonian, physical sector and upper are unchanged. This is a continuation from frozen inputs, not a fresh calculation from integrals.','',
        '| Complete replay | Width (mHa) | Target met |','|---|---:|---|']
    lines += [f'| {name} | {row["width_mHa"]:.9f} | {row["target_met"]} |' for name,row in exact.items()]
    if matched:
        lines += ['',f'The two 300-second continuations started from the same checkpoint. The four-direction export has an exactly checked width **{matched["gain_mHa"]:.9f} mHa** smaller than the unchanged-family export. Their selector, construction and replay costs are additional. This is one conditional warm comparison, not evidence of an end-to-end speedup.']
    common=iterations['last_common_logged_iterate']
    lines += ['',f'**That comparison does not isolate a benefit from the new correlations.** The four-direction run completed 1,408 iterations and the control 1,201. At every positive common logged iteration, the four-direction predicted width was worse. At iteration {common["iteration"]}, the control predicted {common["control_predicted_width_mHa"]:.9f} mHa and the four-direction run {common["four_predicted_width_mHa"]:.9f} mHa. The equal-time outcome includes different throughput; its reproducibility has not been established. These intermediate scores are numerical diagnostics, not separately replayed certificates.']
    lines += ['',f'The new exact singlet residual allowance is **{dense["singlet_residual_mHa"]:.9f} mHa**, or {residual_diagnosis["residual_share_percent"]:.2f}% of the complete width. Keep this exported scalar, positive squares, upper and other allowances fixed. The existing norm certificate bounds R between -eta and +eta, so replacing only its scalar lower can improve the result by at most 2 eta. Even that optimistic substitution leaves **{residual_diagnosis["optimistic_width_after_residual_only_replacement_mHa"]:.9f} mHa**. Residual-only repair of this fixed identity cannot meet 1.6 mHa. This does not limit reoptimization or the full family.']
    lines += ['',
        'Interpreted in the current orbital labels, the supplied four-term H8 polynomial is already contained in both current H8 and H12 cones. Exact rational coordinates reproduce it in existing highest-weight blocks, and an exact spin-average identity handles its adjoint. The external source bundle was not supplied: its negative dual evaluation, the physical orbital-basis match, and any required transport were not checked. The H12 searches use their own frame-bound inputs.','',
        '| Numerical experiment | Best predicted width (mHa) |','|---|---:|']
    lines += [f'| {name} | {row["predicted_width_mHa"]:.9f} |' for name,row in proposals.items() if row['predicted_width_mHa'] is not None]
    lines += ['',
        'The existing solver already selected exports by the residual-penalized lower. The new fixed-Gram ideal LP changed its internal coefficient fit; it improved the warmed control by about 0.107 mHa. The block-scale LP exhausted its time limit and retained its source. A second numerical algorithm and a full accepted-L1 conic formulation were both exercised. Their time-limited outputs did not beat the best candidate and are not optimality proofs.','',
        'The selected-channel prototype constructs paired cubic maps through three quartic contractions and updates the normal operator with a Woodbury solve. Four directions add 16 Gram entries; each direction can contain 792 terms. The numerical map agrees with direct polynomial evaluation to approximately 3e-15. The reference accepting path independently expands the actual integer factors.','',
        'The later eight-direction and mixed-spin experiments preserved earlier directions. The mixed diagnostic used a 2,520-by-2,520 moment matrix and selected four further directions, each potentially containing 2,520 terms. A small final rank does not erase those discovery or verification costs. No complete global cubic coefficient map was built for these additions, but the existing 320,543 coefficient rows and all original prepared maps remain dependencies.','',
        f'New metered subprocess stages total **{result["metered_stage_wall_seconds"]/60:.2f} minutes**, with peak child RSS **{result["peak_metered_child_RSS_bytes"]/1e9:.3f} GB**. This includes failed searches, retries and byte-verified input caching after macOS offloading caused read delays. It excludes inherited preparation/discovery, unmetered early unit-test/import time, and lightweight receipt analyses. Those exclusions prevent a complete cold-cost or speedup claim. The prior campaign ledger is linked in the machine-readable accounting. No external compute was purchased. GitHub operations stayed stopped.','',
        'All 28 focused regressions passed. These include exact comparisons between the contraction and original CAR algebra, original refusal paths, arbitrary-size integers, an independent literal-ladder oracle on all 256 local occupation states, and a complete two-site/two-electron singlet-basis check. The numerical selector’s negative eigenvalues are not energy gains or exact counterexamples to the present family. The trace diagnostic does not perform exact affine/nullspace repair or exact PSD verification. A stalled search is not an obstruction.','',
        'The remaining scientific items are H12 accuracy or an exact family obstruction; useful selection at complete measured cost; reducing the retained global bookkeeping and accepting cost; frozen-rule transfer; independent fresh-environment reproduction; and physical-model validation beyond finite-Hamiltonian solver precision. The existing H8/H10, water, orbital-expansion and coupling milestones remain preserved.','',
        'See [the mathematics](MATHEMATICS.md), [declared experiments](PROTOCOL.md), [run commands](RUNBOOK.md), and the machine-readable [accounting](../../results/acceptance_channels_20260915/accounting.json).']
    if paired:
        position=next(i for i,line in enumerate(lines) if line.startswith('The remaining scientific items'))
        lines[position:position]=[f'The experimental exact paired-contraction replay produced identical input witnesses, exact endpoints and singlet residual allowance. Its complete time was **{paired["contracted_complete_seconds"]:.2f} seconds**, against **{paired["original_complete_seconds"]:.2f} seconds** for the original replay ({paired["original_time_divided_by_contracted_time"]:.3f}× ratio). This is one replay of each path on this host; it is not a fresh solve or an end-to-end speedup. All unmatched factors, sector checks, orbital allowances and original validations remain on the reference path. The paired contraction is an exactly equivalent component with demonstrated scope, not a solution to the 320,543-row global bookkeeping problem.','']
    if component:
        position=next(i for i,line in enumerate(lines) if line.startswith('The remaining scientific items'))
        lines[position:position]=[f'The isolated molecular paired component also produced the identical exact polynomial and factor statistics. Original CAR evaluation took **{component["original_evaluation_seconds"]:.3f} seconds**; direct integer contraction took **{component["contracted_evaluation_seconds"]:.3f} seconds**, a **{component["ratio"]:.2f}×** component ratio. This local ratio must not be applied to the whole verifier. The full replay comparison also includes timing variation in unchanged steps.','']
    if trace:
        position=next(i for i,line in enumerate(lines) if line.startswith('The remaining scientific items'))
        lines[position:position]=[f'The retried numerical singlet-trace diagnostic completed. Its trace normalization/ideal defect was {trace["trace_ideal_defect"]:.3g}, but the candidate still acted on the trace kernel at up to {trace["maximum_null_action"]:.3g}. Even the positive-face-only mixture would weaken the proposed dual energy by {1000*(trace["face_only_energy_after_mixing_Ha"]-trace["approximate_dual_energy_Ha"]):.6f} mHa, and it does not repair that kernel defect. The suggested {trace["unverified_family_floor_mHa"]:.6f} mHa family floor is **unverified** and is not an accepted obstruction. Exact affine/nullspace repair and PSD proof remain missing.','']
    position=next(i for i,line in enumerate(lines) if line.startswith('The remaining scientific items'))
    lines[position:position]=[f'A small exact follow-up identified **{null_channels["independent_null_vectors_total"]} physical singlet-null directions** in the two actual 216-dimensional dictionaries: S_plus a_p_beta and S_plus a_p_alpha_dagger. Their commutators, dictionary coordinates and independence are checked exactly. This identifies part of the trace kernel without relying on numerical eigenvectors. It does not repair the dual, establish which null constraints are forced by the current ideal map, or prove a family obstruction.','']
    (ROOT/'research/acceptance_channels_20260915/REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:result[k] for k in ('best_new_exact','matched_comparison','metered_stage_wall_seconds','completed_metered_stages')}),flush=True)


if __name__=='__main__':main()
