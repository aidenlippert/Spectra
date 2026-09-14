"""Collect completed evidence and preserve the distinction between the goals."""
from fractions import Fraction as F
import json
from pathlib import Path
from research.reconstruction_compression_20260914.inputs import ROOT,OUT,sha,dump
from research.reconstruction_compression_20260914.transfer import hashes

def read(p):return json.loads(Path(p).read_text())

def collect():
    inputs=read(OUT/'frozen_inputs.json');runs=[read(p) for p in sorted((OUT/'runs').glob('*.json'))]
    pending=[r['name'] for r in runs if r['status']=='starting' and not r['name'].startswith('assemble_report')]
    if pending:raise ValueError(('Unfinished runs',pending))
    frozen=read(OUT/'transfer_freeze.json')
    if hashes()!=frozen['source_sha256']:raise ValueError('Transfer source freeze changed')
    preservation=read(OUT/'preservation_before.json')['files']
    changed=[p for p,v in preservation.items() if not (ROOT/p).is_file() or sha(ROOT/p)!=v['sha256']]
    if changed:raise ValueError(('Preserved bytes changed',changed))
    candidates=[]
    for path in sorted((OUT/'candidates').glob('*/certificate.json')):
        d=read(path.parent/'discovery.json');check=read(path.parent/'interval.json');cert=read(path)
        if check['certificate_sha256']!=sha(path):raise ValueError('Candidate replay hash mismatch')
        case=d['case'];U=F(check['upper_Ha']);L=F(check['lower']);method='coefficient_L1'
        spectral=path.parent/'spectral_replay'/'receipt.json'
        spectral_path=None;best_replay_seconds=check['replay_seconds'];proof_bytes=path.stat().st_size
        if spectral.exists():
            sp=read(spectral)
            if sp['certificate_sha256']!=sha(path):raise ValueError('Spectral replay hash mismatch')
            if F(sp['lower'])>L:
                L=F(sp['lower']);method='wedge_spectral';spectral_path=str(spectral)
                best_replay_seconds=sp['wall_seconds']
                proof_bytes+=(path.parent/'spectral'/'witness.json').stat().st_size
        source='h6' if case=='h6_actual' else 'h8' if case=='h8' else None
        solver_run=None
        for run in runs:
            command=run.get('command',[])
            if 'solve' not in command:continue
            index=command.index('solve')
            if command[index+1:index+4]!=[case,d['mode'],str(d['rank_cap'])]:continue
            tag=command[command.index('--tag')+1] if '--tag' in command else ''
            if path.parent.name==f"{case}_{d['mode']}_{d['rank_cap']}{tag}":
                solver_run=run;break
        row={'name':path.parent.name,'case':case,'mode':d['mode'],'rank':d['rank_cap'],'teacher_assisted':d['teacher_used'],
             'certificate':str(path),'certificate_sha256':sha(path),'certificate_bytes':path.stat().st_size,
             'lower_Ha':str(L),'upper_Ha':str(U),'width_mHa':float(1000*(U-L)),
             'target_1p6mHa_met':U-L<=F(1,625),'best_residual_method':method,'spectral_replay':spectral_path,
             'coefficient_interval_mHa':check['width_mHa'],'bare_b_Ha':cert['b'],
             'coefficient_remainder_allowance_Ha':check['residual_l1'],'best_residual_loss_Ha':str(F(cert['b'])-L),
             'required_lower_Ha':str(U-F(1,625)),'lower_shortfall_Ha':str(U-F(1,625)-L),
             'gram_entries':d['reduced_gram_entries'],'full_dictionary_entries':d['full_gram_entries'],
             'optimization_entry_reduction':d['full_gram_entries']/d['reduced_gram_entries'],
             'factor_rows':check['factor_rows'],'projected_coefficient_map_nonzeros':d['projected_map_nonzeros'],
             'search_and_export_seconds':d['total_seconds'],'coefficient_lower_replay_seconds':check['replay_seconds'],
             'best_lower_replay_seconds':best_replay_seconds,'total_lower_proof_bytes':proof_bytes,
             'peak_solver_process_RSS_bytes':solver_run.get('peak_child_RSS_bytes') if solver_run else None,
             'exact_equality_elimination':d.get('exact_quartic_ideal_elimination',False),
             'all_coefficient_rows':d['all_CAR_rows'],'active_equations':d.get('active_equations_after_exact_elimination',d['all_CAR_rows']),
             'uses_full_dictionary_map_and_moment_preparation':True,'full_many_body_space_enumerated':False}
        if source:row['deterioration_from_preserved_lower_Ha']=str(F(inputs[source]['lower_Ha'])-L)
        candidates.append(row)
    best={case:min((r for r in candidates if r['case']==key and not r['teacher_assisted']),key=lambda r:r['width_mHa'])
          for case,key in [('h6','h6_actual'),('h8','h8')]}
    transfer=next(r for r in candidates if r['case']=='fresh_h6_2p07')
    uppers={case:read(OUT/'upper'/case/'comparison.json') for case in ('h6','h8')}
    bundles={case:read(OUT/'bundles'/f'{case}.json') for case in ('h6','h8')}
    accounting={'completed_process_wall_seconds':sum(r.get('wall_seconds',0) for r in runs),
                'completed_child_CPU_seconds':sum(r.get('child_user_seconds',0)+r.get('child_system_seconds',0) for r in runs),
                'initial_aggregate_budget_seconds':2400,'peak_process_RSS_bytes':max(r.get('peak_child_RSS_bytes',0) for r in runs),
                'over_1GiB_target':[{'name':r['name'],'bytes':r['peak_child_RSS_bytes']} for r in runs if r.get('peak_child_RSS_bytes',0)>2**30],
                'failed_or_timed_out':[{'name':r['name'],'status':r['status'],'wall_seconds':r.get('wall_seconds')} for r in runs if r['status'] in ('failed','timeout','runner_error')],
                'completed_run_count':sum(r['status']!='starting' for r in runs),'pending_runs':pending,
                'scope':'Includes bounded failed attempts, tests, preparation, numerical discovery and replay. Excludes editing, file browsing, final inventory, and inherited historical discovery. Some independent processes overlapped; wall sums are not elapsed calendar time.'}
    result={'campaign':'reconstruction_compression_20260914','main_compressed_lower_target_achieved':all(r['target_1p6mHa_met'] for r in best.values()),
            'original_milestones_preserved':True,'preserved_file_count':len(preservation),'preserved_file_changes':changed,
            'transfer_rule_hashes_unchanged':True,'frozen_inputs':inputs,'candidates':candidates,'best_independent_compressed':best,
            'fresh_transfer':transfer,'upper_comparison':uppers,'fast_complete_bundles':bundles,'accounting':accounting,
            'remaining_obstruction':'No accuracy-preserving compressed lower found in the tested families and budgets. Dense projected coefficient maps and full moment/map preparation remain. This is not an impossibility theorem.',
            'historical_discovery_seconds':{'full_cubic_lower_h6':259.309,'full_cubic_lower_h8':721.563,
                'accepted_h6_MPS_final_run_internal':8.5489,'accepted_h8_MPS_successful_warm_path_processes':414.36},
            'historical_cost_limit':'These are retained historical measurements; they do not account for every earlier unsuccessful attempt. Teacher-assisted rows inherit the full-cubic discovery cost.',
            'focused_tests':'12 passed; primitive rational enclosure tests, exact CAR/congruence/degree-six tests, modular rank/refusal tests, Hamiltonian/sector/MPS mutations and endpoint refusal/fallback.',
            'scope':'Supplied rational finite-basis electronic Hamiltonians. No experimental-accuracy or general-scalability claim; CH2 model result preserved.'}
    if accounting['completed_process_wall_seconds']>2400:raise ValueError('Campaign process budget exceeded')
    return result

def write_report(r):
    rows=r['candidates'];best=r['best_independent_compressed'];a=r['accounting'];u=r['upper_comparison'];b=r['fast_complete_bundles'];fresh=r['fresh_transfer']
    lines=['# Reconstruction-preserving proof compression — completed bounded pass','',
      f"The main target remains open. The strongest independently constructed compressed intervals are **{best['h6']['width_mHa']:.6f} mHa on H6** and **{best['h8']['width_mHa']:.6f} mHa on corrected H8**. Both exceed 1.6 mHa. The preserved 0.158560 and 1.176860 mHa milestones and the CH2 model result remain unchanged.",'',
      f"The useful accepting-path result is a rigorous interval MPS verifier: H8 takes **{u['h8']['interval_replay_seconds']:.3f} s**, compared with **{u['h8']['integer_replay_seconds']:.3f} s** for the integer checker, an observed **{u['h8']['speed_ratio']:.2f}×** ratio on the same rounded state and declared endpoint. The endpoint allowance is **{u['h8']['charged_endpoint_allowance_mHa']:.8f} mHa**. Exact equality remains ambiguous and requires refusal or integer fallback.",'',
      '## Complete accepting results','',
      '| Case | Preserved exact interval (mHa) | Fast-upper + preserved-lower interval (mHa) | New complete accepting time (s) |','|---|---:|---:|---:|']
    for case in ('h6','h8'):lines.append(f"| {case.upper()} | {r['frozen_inputs'][case]['width_mHa']:.9f} | {b[case]['width_mHa']:.9f} | {b[case]['seconds']:.3f} |")
    lines+=['','These new complete bundles recheck the inherited full-cubic lower proof. They establish a cheaper verifier with a charged, tiny upper allowance; they do not establish a cheaper lower constructor. Both upper paths were compared with the same endpoint ceil(2^24 U)/2^24. The main lower searches used the original exact U throughout.','',
      '## Full lower-plus-frozen-upper comparisons','',
      'Every row below has an independently replayed rational lower bound. “Passed” in a process receipt only means the program exited successfully; the interval target is assessed separately. Search time includes numerical construction/solve/export and excludes the shared preparation listed below. Rank is a cap per odd block; quadratic repair blocks stay full.','',
      '| Case / rule | Rank | Gram entries | Certified interval (mHa) | Search/export (s) | Accepting lower replay (s) | Peak solver MiB |','|---|---:|---:|---:|---:|---:|---:|']
    names={'teacher':'teacher, unpaired','unpaired':'MPS, unpaired','simple':'simple, paired','low':'MPS, paired','cuts':'nested cuts, paired'}
    for row in rows:
        if row['case']=='fresh_h6_2p07':continue
        label=('H6' if row['case']=='h6_actual' else 'H8')+' / '+names[row['mode']]
        if row['exact_equality_elimination']:label+='; exact elimination'
        if row['name'].endswith('_eliminated'):label+='; indirect solve failed numerically'
        peak=row['peak_solver_process_RSS_bytes']
        lines.append(f"| {label} | {row['rank']} | {row['gram_entries']:,} | {row['width_mHa']:.6f} | {row['search_and_export_seconds']:.3f} | {row['best_lower_replay_seconds']:.3f} | {peak/2**20:.1f} |" if peak else f"| {label} | {row['rank']} | {row['gram_entries']:,} | {row['width_mHa']:.6f} | {row['search_and_export_seconds']:.3f} | {row['best_lower_replay_seconds']:.3f} | unavailable |")
    prep={case:read(OUT/'prepared'/key/'receipt.json') for case,key in [('h6','h6_actual'),('h8','h8')]}
    lines+=['',f"Shared Hamiltonian-map and MPS-moment preparation costs {prep['h6']['map_seconds']+prep['h6']['moment_seconds']:.3f} internal seconds for H6 and {prep['h8']['map_seconds']+prep['h8']['moment_seconds']:.3f} for H8, plus interpreter/import overhead in the run ledger. Every standalone candidate must pay that preparation. The initial H6 quotient-frame trial, replaced to match the teacher convention, is also charged in the campaign total.",'',
      'The two useful rank-32 H6 rows use the existing wedge spectral residual checker; other rows use the full CAR coefficient remainder. Their separate spectral-proposal costs remain in the ledger, and the accepting column includes the selected checker’s complete lower replay. The result index includes proof bytes and exact endpoints. No omitted degree-six terms or floating-point feasibility claims substitute for acceptance.','',
      f"The best searches use approximately {best['h6']['optimization_entry_reduction']:.1f}× fewer Gram entries on H6 and {best['h8']['optimization_entry_reduction']:.1f}× fewer on H8 than their uncompressed dictionaries. That does not give a matched-accuracy speedup. In particular, H8’s rank-24 projected coefficient map has {best['h8']['projected_coefficient_map_nonzeros']:,} nonzeros and exceeded the memory target.",'',
      'At rank eight, MPS-guided paired directions improve the interval relative to simple paired directions on both cases. The cut maps improve less. The tiny unpaired subspaces return almost the same bounds as the quadratic/simple control. This supports the relevance of state information, while showing that the tested small spaces lose too much reconstruction capability. It does not prove that larger or differently structured subspaces must fail.','',
      '## Exact equality reduction and the cut experiment','',
      'The paired anticommutator construction cancels degree-six terms before dense projection. An exact modular full-column-rank calculation proves that the quartic sector multipliers must then be zero in an exact reconstruction. H6’s active equations drop from 4,859 to 499 and its multiplier dimension from 499 to 25; H8 drops from 28,461 to 1,525 equations and from 1,525 to 41 multipliers. Equivalence applies to the paired model’s zero-remainder feasible set, not every residual-corrected certificate. It did not resolve its accuracy or convergence limitations.','',
      'The cut experiment constructs nested fixed maps at two orbital partitions, with a per-node dimension cap of two and an independent root rank. All root positivity constraints share global CAR moments and the full molecular Hamiltonian. Unknown states are never restricted to the MPS. However, this prototype still constructs full dictionary moments/maps and exports expanded factors. It is a finite operator-map restriction, not an efficient implementation of a molecular RDM consistency hierarchy. That higher-upside part remains unresolved.','',
      '## Frozen transfer','',
      f"The generation and search rules were hashed before a new H6 geometry at **2.07 Å** was generated. The independent bond-48 MPS, paired rank-32 lower, exact equality reduction and 65-second solver budget produced a complete interval of **{fresh['width_mHa']:.6f} mHa**. The target {'was met' if fresh['target_1p6mHa_met'] else 'was not met'}. The frozen transfer uses coefficient-L1 acceptance, which gave 4.011289 mHa for the corresponding original-fixture construction; it does not add the later spectral-residual improvement. No FCI trial state or full-cubic teacher solve was used on this geometry. The exact interval concerns the rational electronic Hamiltonian; nuclear/geometry/model errors are separate.",'',
      '## Accounting, validation, and preservation','',
      f"{a['completed_process_wall_seconds']:.3f} seconds of bounded process wall time are recorded against the initial 2,400-second campaign budget, including failures, preparation, tests and replay. Independent jobs sometimes overlapped, so this sum is not calendar elapsed time. Editing, file inspection, final inventory and inherited discovery are outside that sum. The ledger contains {a['completed_run_count']} completed process receipts.",'',
      'The teacher rows also inherit the old full-cubic discovery (259.309 s for H6 and 721.563 s for H8). Frozen MPS construction is another inherited dependency: the retained H8 warm path took 414.36 recorded seconds. These historical figures do not account for every earlier failed attempt. None is silently counted as ordinary Hamiltonian input or free preprocessing.','',
      'Two numerical builds exceeded the 1 GiB working target: '+', '.join(f"{x['name']} ({x['bytes']/2**30:.3f} GiB)" for x in a['over_1GiB_target'])+'. Their costs are retained; subsequent smaller cases and the frozen transfer stayed below that target. The indirect linear solver trial failed to converge to useful positivity and is retained as an 80,368 mHa failure, not a successful compression.','',
      'Twelve focused tests pass. They cover exact CAR reconstruction and degree-six retention, adjoint pairing, numerical moments against an explicit tiny oracle, positive congruences, modular rank and its refusal case, MPS/sector/Hamiltonian mutations, rigorous primitive enclosures, bad upper endpoints, and equality fallback. Initial implementation/test and report-assembly failures remain in the prospective logs.','',
      f"All **{r['preserved_file_count']:,} inherited files** match their starting hashes. The fresh rule hashes are unchanged. No external resources were provisioned or published, and the CH2 model result was preserved.",'',
      'The result supports keeping the faster rigorous upper checker. The lower-bound research still needs a representation that carries enough cross terms and constructs the coefficient map compactly. Increasing dense reduced rank already grows memory before H8 approaches the accuracy target. No exact energy-family obstruction was proved in this pass.','',
      f"[Mathematical derivation and limitations]({ROOT/'research/reconstruction_compression_20260914/MATHEMATICS.md'}) · [Exact result index]({OUT/'final_result.json'}) · [Accounting audit]({OUT/'audit.json'})",'']
    (ROOT/'research/reconstruction_compression_20260914/REPORT.md').write_text('\n'.join(lines))

def main():
    r=collect();dump(OUT/'final_result.json',r);dump(OUT/'audit.json',{'preserved_file_count':r['preserved_file_count'],
        'preserved_file_changes':r['preserved_file_changes'],'transfer_rule_hashes_unchanged':r['transfer_rule_hashes_unchanged'],**r['accounting']})
    write_report(r)
    print(json.dumps({'best_compressed_mHa':{k:v['width_mHa'] for k,v in r['best_independent_compressed'].items()},
                     'fresh_mHa':r['fresh_transfer']['width_mHa'],'completed_process_seconds':r['accounting']['completed_process_wall_seconds'],
                     'preserved_files':r['preserved_file_count']},indent=2))

if __name__=='__main__':main()
