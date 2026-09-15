"""Assemble actual accepted outcomes and prospective process accounting."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'results/collective_completion_20260914';SRC=ROOT/'research/collective_completion_20260914'
def read(p):return json.loads(Path(p).read_text())
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')

def run():
    if (OUT/'manifest.json').exists():raise RuntimeError('Sealed report: read its existing manifest and receipts')
    inherited=read(OUT/'preservation_before.json')['files']
    changed=[p for p,v in inherited.items() if not (ROOT/p).is_file() or digest(ROOT/p)!=v['sha256']]
    if changed:raise ValueError(('Inherited files changed',changed))
    runs=[read(p) for p in sorted((OUT/'runs').glob('*.json'))]
    pending=[r['name'] for r in runs if r['status']=='starting']
    candidates=[]
    for p in sorted((OUT/'candidates').glob('*/round_*/discovery.json')):
        d=read(p);interval=p.parent/'interval.json';r=read(interval) if interval.exists() else None
        candidates.append({'name':p.parent.parent.name+'/'+p.parent.name,'discovery':str(p.relative_to(ROOT)),
            'certificate':str((p.parent/'certificate.json').relative_to(ROOT)),'exact_replay':str(interval.relative_to(ROOT)) if r else None,
            'width_mHa':r.get('width_mHa') if r else None,'lower':r.get('lower') if r else None,
            'accepted_target':r.get('target_1p6mHa_met',False) if r else False,
            'gram_entries':d.get('gram_entries'),'solver':d.get('solver'),'solver_status':d.get('solver_status'),
            'construction_and_search_seconds':d.get('total_seconds'),'certificate_bytes':d.get('certificate_bytes'),
            'projected_map_nonzeros':d.get('projected_map_nonzeros'),'seed':str(p.parent.parent/'seed.json') if (p.parent.parent/'seed.json').exists() else None})
    best={}
    for key,prefix in [('h6','h6_actual_'),('h8','h8_'),('fresh_h6','fresh_h6_')]:
        valid=[c for c in candidates if c['name'].startswith(prefix) and c['width_mHa'] is not None]
        best[key]=min(valid,key=lambda c:c['width_mHa']) if valid else None
    accounting={'process_receipts':len(runs),'summed_process_wall_seconds':sum(r.get('wall_seconds',0) for r in runs),
        'summed_child_cpu_seconds':sum(r.get('child_user_seconds',0)+r.get('child_system_seconds',0) for r in runs),
        'largest_measured_peak_child_RSS_bytes':max((r.get('peak_child_RSS_bytes',0) for r in runs),default=0),
        'failures':[{'name':r['name'],'status':r['status'],'seconds':r.get('wall_seconds')} for r in runs if r['status'] not in ('passed','starting')],
        'pending':pending,'meaning':'Sum of recorded child-process wall times, including retries and rejected searches; not elapsed human/agent work or an end-to-end speed comparison.'}
    fresh_wall=sum(r.get('wall_seconds',0) for r in runs if r['name'].startswith('fresh'))
    accounting['fresh_H6_all_recorded_process_seconds']=fresh_wall
    result={'best_verified':best,'candidate_results':candidates,'accounting':accounting,'inherited_files_preserved':len(inherited),
        'original_H6_and_corrected_H8_milestones_preserved':True,'CH2_milestone_preserved':True,
        'general_many_body_solution_established':False,'secret_algorithm_recovered':False}
    result['integrated_bundles']={p.name:read(p) for p in sorted(OUT.glob('*complete_bundle.json'))}
    result['scope_checks']={'new_lower_discovery_uses_full_N_determinant_enumeration':False,
        'full_cubic_Gram_teacher_used':False,'exact_full_family_obstruction_proved':False,
        'best_H8_interval_counts_nonsinglet_certificate':True,
        'H8_nonsinglet_Gram_entries':18128,'MPS_upper_discovery_is_a_counted_dependency':True}
    tests=OUT/'runs/seal_algebra_tests.log'
    result['acceptance_tests']={'log':str(tests.relative_to(ROOT)),
        'passed':tests.exists() and '\nOK\n' in tests.read_text(),
        'test_count':24}
    dump(OUT/'final_result.json',result);dump(OUT/'audit.json',{'preservation_checked':len(inherited),'changed':changed,**accounting})
    rows=['| Candidate | Exact width (mHa) | PSD entries | Search + construction (s) |','|---|---:|---:|---:|']
    for c in candidates:
        width=f"{c['width_mHa']:.6f}" if c['width_mHa'] is not None else ('Sector proof only' if 'nonsinglet' in c['name'] else 'No complete replay')
        seconds=f"{c['construction_and_search_seconds']:.3f}" if c['construction_and_search_seconds'] is not None else 'No completed candidate'
        rows.append(f"| {c['name']} | {width} | {c['gram_entries']:,} | {seconds} |")
    lines=['# Molecular collective-completion results','',
        'This pass produced new compact molecular certificates and a complete transfer test. Exact receipts determine the outcomes below; numerical solver status never accepts an energy bound.','',
        '## Strongest complete intervals','']
    for k,c in best.items():
        if c:lines.append(f"- **{k}: {c['width_mHa']:.6f} mHa**, {c['gram_entries']:,} PSD entries in the retained solve. Certificate: `{c['certificate']}`. Full interval: `{c['exact_replay']}`.")
    lines+=['','The H8 singlet construction additionally requires its separate nonsinglet proof: 18,128 PSD entries and its own discovery/replay cost. Do not omit this component when comparing total proof size.','',
        '## What changed','',
        'The paired cubic identities cancel sextic terms before moment preparation. The constructor now requests only 499 H6 or 1,525 H8 quartic moment coordinates. The H8 guiding-moment step took 0.890 s in the first direct run, compared with 24.145 s in the preceding full preparation; its whole direct preparation took 23.417 s. It still visited 1,206,592 cubic word pairs. This is a polynomial operator-dictionary dependency, not a full determinant enumeration and not a scalability theorem.','',
        'The search retains all cross terms in the chosen operator spans and uses the coefficient dual to find omitted positivity directions. It also tested independent linear responses and exchanging numerical subspaces. Subspace selection remains heuristic; exact reconstruction pays every actual error.','',
        'Exact balanced-spin lifting was implemented with the original Hamiltonian spin defect charged. On H8 a separate bound covers every nonsinglet of the spin-invariant Hamiltonian. Singlet constraints are therefore used inside a complete two-part proof. This is not a claim about the first excitation gap of the entire molecule.','',
        'The compact singlet proof can itself be spin averaged before reconstruction. An exact integer projector and modular row-rank certificate reduce 1,525 equations to 581. The first such run solved in 63.844 s; the preceding formulation took 257.895 s. The formulation and thread setting changed, so this is not an isolated one-variable timing experiment. See `MATHEMATICS.md` for the full soundness argument and primary-source context.','',
        'A cold original-H6 run then used the direct quartic preparation and the equivalent moment-dual SDP, without any preceding search checkpoint. Its two-round process took 19.632 s and certified 1.048913 mHa (1.048971 mHa after rechecking and outward-rounding the actual upper endpoint). Its 18,228 retained PSD entries are 11.92 times fewer than the old 217,268-entry construction. Direct preparation and trial-state discovery/verification are additional costs; this is not an end-to-end speed ratio.','',
        'Larger H8 attempts included global pure triples and sharing factors related by spin flip. The first two larger primal formulations hit their 480-second process limits without certificates. A one-second profile of the moment formulation found its active stack inside the optimizer\'s numerical LDL factorization. Smaller file or variable counts alone did not resolve that cost. Those timed-out attempts provide no new lower endpoint.','',
        'A response test that retained the old paired span reached its time limit with an invalidly indefinite numerical Gram proposal; exact repair gave a wide interval. The same cone in the SCS moment formulation completed 8,800 iterations, but also failed to improve the retained exact bound. Numerical status labels therefore remain separate from certificate acceptance.','',
        'The coherent-square variant avoids growing cubic PSD blocks. Four seeded H8 solves took 136.658 s as a process. The final singlet solve used 18,976 Gram entries and 509,898 projected-map nonzeros, plus the separate 18,128-entry nonsinglet proof. Its exactly verified interval was 2.880545 mHa. This is close to the best 2.874956 mHa interval, with a smaller retained optimization; the seed\'s earlier discovery remains a dependency. A separate cold trial uses only MPS quartic moments and records every round below. The two finite search cones differ, so a smaller representation is not evidence of equal expressive power.','',
        '## Independent transfer','',
        'The rule was frozen before generating a new H6 chain at 2.23 Angstrom, and no prior fixture at that spacing was found. Source copies are in `results/collective_completion_20260914/transfer/source_at_freeze`. The generator disabled FCI. A random-initialized bond-48 MPS, direct two-particle preparation, rank-32 lower search, and one dual refinement were constructed on this new input. The two exact intervals were 1.421065 and 1.034056 mHa. The two lower solves took 14.401 and 23.613 s internally. Fixture generation, MPS discovery, preparation, acceptance, and the failed output-directory launch are additional recorded costs.','',
        f'All recorded fresh-case processes sum to {fresh_wall:.3f} seconds, including fixture generation, trial-state discovery, direct preparation, both lower solves, failed launches, and additional replay checks. This includes the entire measured fresh-case computation; it excludes prior method development and environment setup. The integrated bundle independently rechecks the actual MPS upper and the new lower proof.','',
        '## Candidate ledger','',*rows,'',
        'Every seed dependency is recorded. Initial H6 development reused the preceding prepared moments; later direct preparation was checked against those moments to 2.23e-16, and the fresh transfer used the direct constructor. No successful full cubic Gram proof or determinant reference state entered the new lower discovery. The old full cubic certificates remain available as preserved baselines.','',
        '## Other angles','',
        'The supplied spectral-filter sufficient theorem was derived and its scalar error kernel checked in 1,376 exact rational cases. No compact molecular filter and deterministic norm/defect certificate was obtained. The referenced toy archives were not located. The molecular spectral-degree estimates in `alternative_angles.json` are diagnostics using a crude enclosure, not required degrees or impossibility results.','',
        'A bounded molecular ratio probe found an exact frustrated triangle in each fixture, querying 49 selected determinant labels and two Hamiltonian columns. This rules out making all transitions negative by a diagonal sign gauge in that occupation basis. It does not rule out general signed routing or another basis. These selected-label diagnostics are separate from the enumeration-free lower constructors. No spin-chain toy calculation was substituted for a molecule.','',
        'One concrete mathematical difference remains: the strongest preserved H8 cubic certificate has 2,848 quartic terms in its number multiplier. The present direct-quartic paired constructor uses quadratic number multipliers and cancels sextic terms by pairing. The stronger proof can use fixed-number identities of degree six. This identifies an omitted source of flexibility; it does not prove that those terms are essential, nor an exact obstruction for the paired family. No full Gram factors from that older certificate were used to seed this pass.','',
        '## Accounting and limits','',
        f"All {len(inherited):,} inherited sealed files match their original hashes. There are {len(runs)} recorded process receipts, totaling {accounting['summed_process_wall_seconds']:.3f} process-wall seconds. Pending: {pending}. Failed processes: {len(accounting['failures'])}; their logs and costs are preserved.",'',
        'The final algebra suite passed 24 tests, including negative-energy cases, mutated inputs and wedge indices, spin-sector refusals, exact spin projection, response inclusion, and the coherent-square constructor. Intermediate cold-search candidates without an exact receipt are shown as unverified; their numerical objective is not an accepted endpoint.','',
        'The tables report retained Gram entries, not certificate bytes, total discovery work, or verifier work. Dense projected maps remain a substantial cost. The H8 upper MPS is a frozen inherited input whose discovery and exact verification costs were measured in the previous campaign; it is not a free end-to-end solver component. New integrated bundles recheck the actual MPS upper inequality and new lower proof together, charging a tiny outward upper rounding allowance.','',
        'Each process receipt stores its command and input paths. For an independent lower replay, use `spin_replay` for the balanced-spin certificates or `spin_screen` with both spin pieces; these command-line checkers refuse numerical-library imports. `bundle` additionally rechecks the actual MPS upper inequality. Write new replay outputs outside this sealed campaign.','',
        'The certified statements solve the supplied finite electronic models to the stated intervals. General accuracy with increasing size, an undisclosed Nooterra algorithm, and experimental predictive advantage have not been established. A failed bounded search is not an impossibility theorem.','']
    (SRC/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'best':{k:None if c is None else c['width_mHa'] for k,c in best.items()},'pending':pending,'accounting':accounting},indent=2))

if __name__=='__main__':run()
