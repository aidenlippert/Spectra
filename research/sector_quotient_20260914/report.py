"""Assemble measured outcomes and seal this research pass after all runs stop."""
from fractions import Fraction as F
import datetime
import json
import math
from pathlib import Path
from research.reconstruction_compression_20260914.inputs import sha, dump
from research.sector_quotient_20260914.budget import ROOT, OUT

SOURCE = ROOT/'research/sector_quotient_20260914'


def run():
    if (OUT/'manifest.json').exists(): raise RuntimeError('Sealed campaign')
    runs = [json.loads(p.read_text()) for p in sorted((OUT/'runs').glob('*.json'))]
    pending = [r['name'] for r in runs if r['status'] == 'starting']
    if pending: raise ValueError(('Numerical processes still pending', pending))
    inherited = json.loads((OUT/'preservation_before.json').read_text())['files']
    changed = [p for p, rec in inherited.items() if not (ROOT/p).is_file() or sha(ROOT/p) != rec['sha256']]
    if changed: raise ValueError(('Changed inherited files', changed))
    requests = json.loads((OUT/'request_sources.json').read_text())
    if any(sha(p) != digest for p, digest in requests['files'].items()): raise ValueError('Request attachment changed')
    frozen = json.loads((OUT/'frozen_inputs.json').read_text())['h8']; U = F(frozen['upper_Ha'])
    labels = {'fixed_target_repaired': 'Restored freedom, fixed 1.4 mHa proposal',
              'unchanged_spans_repaired': 'Restored freedom, original compact spans',
              'linear_closure_repaired': 'Restored freedom + prescribed linear completion'}
    outcomes = []
    for p in sorted((OUT/'candidates').glob('**/interval.json')):
        rec = json.loads(p.read_text())
        if rec.get('kind') != 'complete_spin_screen_v1': continue
        if F(rec['upper_Ha']) != U or F(rec['width_Ha']) != U-F(rec['lower']):
            raise ValueError('Endpoint or frozen-upper mismatch')
        certificate = p.parent/'certificate.json'
        if not certificate.exists(): raise ValueError('Missing replayed certificate')
        case = str(p.parent.relative_to(OUT/'candidates'))
        if (rec['valid_on'] != 'Entire fixed-N sector' or rec['many_body_states_enumerated'] != 0 or
            not rec['singlet']['spin_twirl'] or
            F(rec['lower']) != min(F(rec['singlet']['lower']), F(rec['all_nonsinglets']['lower']))-F(59,250000000000)):
            raise ValueError('Incomplete sector or spin-defect accounting')
        cert = json.loads(certificate.read_text()); b = F(cert['core']['b'])
        repair = json.loads((p.parent/'repair.json').read_text())
        norms = {k: math.sqrt(float(F(repair[k+'_squared_Frobenius_norm'])))
                 for k in ('W', 'V', 'Z') if k+'_squared_Frobenius_norm' in repair}
        outcomes.append({'case': case, 'label': labels.get(case, case), 'interval': str(p.relative_to(ROOT)),
                         'certificate': str(certificate.relative_to(ROOT)), 'certificate_sha256': sha(certificate),
                         'certificate_bytes': certificate.stat().st_size, 'lower_Ha': rec['lower'],
                         'exported_b_Ha': str(b), 'full_lower_adjustment_Ha': str(b-F(rec['lower'])),
                         'coefficient_Frobenius_norms': norms,
                         'width_Ha': rec['width_Ha'], 'width_mHa': float(F(rec['width_Ha'])*1000),
                         'replay_seconds': rec['replay_seconds'], 'target_met': F(rec['width_Ha']) <= F(1, 625)})
    if not outcomes: raise ValueError('No independently replayed complete molecular result')
    best = min(outcomes, key=lambda r: F(r['width_Ha']))
    old_compact = json.loads((ROOT/'results/collective_completion_20260914/candidates/h8_spin_r32_invariant_dual48/round_0/interval.json').read_text())
    strong = json.loads((OUT/'strong_h8_diagnostic.json').read_text())
    if F(strong['strong_wedge_replay']['lower']) != F(frozen['lower_Ha']):
        raise ValueError('Strong inherited lower no longer matches its replay')
    wall = sum(r.get('wall_seconds', 0) for r in runs)
    cpu = sum(r.get('child_user_seconds', 0)+r.get('child_system_seconds', 0) for r in runs)
    failed = [{'name': r['name'], 'status': r['status'], 'seconds': r.get('wall_seconds')} for r in runs if r['status'] != 'passed']
    audit = {'completed_processes': len(runs), 'pending': pending, 'recorded_process_wall_seconds': wall,
             'recorded_CPU_seconds': cpu, 'peak_process_RSS_bytes': max(r.get('peak_child_RSS_bytes', 0) for r in runs),
             'failed_processes': failed, 'inherited_files_verified': len(inherited), 'changed_inherited_files': changed,
             'request_files_unchanged': True, 'outcomes': outcomes, 'best_new': best,
             'old_compact_width_mHa': old_compact['width_mHa'], 'old_strong_width_mHa': frozen['width_mHa'],
             'new_compact_H8_target_met': best['target_met'], 'general_many_body_solution_achieved': False,
             'exact_dual_obstruction_obtained': False}
    dump(OUT/'audit.json', audit)
    rows = '\n'.join(f"| {v['label']} | {v['width_mHa']:.6f} | {'Yes' if v['target_met'] else 'No'} |" for v in outcomes)
    norm_rows = '\n'.join(f"| {v['label']} | {v['coefficient_Frobenius_norms']['V']:.9g} | {v['coefficient_Frobenius_norms']['Z']:.9g} |"
                         for v in outcomes if 'V' in v['coefficient_Frobenius_norms'])
    improvement = old_compact['width_mHa']-best['width_mHa']
    best_raw_width = float((U-F(best['exported_b_Ha']))*1000)
    best_adjustment = float(F(best['full_lower_adjustment_Ha'])*1000)
    conclusion = (f"The best new complete interval is **{best['width_mHa']:.6f} mHa**. "
                  f"It {'meets' if best['target_met'] else 'does not meet'} the 1.6 mHa target. "
                  f"Its change relative to the retained paired interval is {improvement:+.6f} mHa of narrowing.")
    text = f'''# H8 fixed-number completion: measured outcome

{conclusion}

The supplied reduction is now implemented and checked against the actual strong
H8 proof. The new constructor permits independent cubic Gram blocks and collective
fixed-number cancellation. This is a completed test of the proposed algebraic
change; it is not a general many-body solution or a proof that every compact
representation has been exhausted.

## The actual strong proof confirms the missing freedom

The inherited proof contains 2,848 quartic number-multiplier coefficients. Exact
CAR expansion verifies **W+L2(X4)=-R6**, with the project's actual sign convention.
Projection verifies the corresponding V and contraction-free Z relationships
exactly. Coefficient Frobenius norms are 0.4608828142 for W, 0.1318644498 for V,
and 0.0000141065 for Z. The last quantity is nonzero: it remains covered by the
old exact residual witness, never silently discarded. These are coefficient norms,
not molecular excitation gaps or energy improvements.

The independently replayed inherited full proof remains at
**{frozen['width_mHa']:.6f} mHa** using the frozen upper. It was read for diagnosis
only; no full-proof factors or eigendirections entered the new compact spans.

## Complete H8 comparisons

| Construction | Full width, mHa | At most 1.6 mHa? |
|---|---:|---|
| Preserved strong full-cubic proof | {frozen['width_mHa']:.6f} | Yes |
| Preserved compact paired proof | {old_compact['width_mHa']:.6f} | No |
{rows}

Every new row has a separate standard-library exact replay using the sealed
spin/CAR verifier. It includes the S=0 proof, the separate MS=1 proof covering
all nonsinglets, and the original-H defect 59/250000000000 Ha charged once.
All rows use exactly the same frozen rational upper endpoint. The same valid
upper is reused here; its prior construction and verification are dependencies.

The best new candidate has a raw exported interval of {best_raw_width:.6f} mHa;
the complete residual and original-H spin adjustment adds {best_adjustment:.6f}
mHa. Its remaining distance from the target therefore cannot be attributed
mainly to replay precision. The final improvement belongs to the combined change
with the prescribed linear directions. A matched paired run with those appended
directions was not performed, so the energy gain cannot be assigned uniquely to
unpairing rather than the additional directions or further optimization.

The strict unchanged-span model has 91,752 independent Gram entries, plus
18,128 for the nonsinglet proof. The paired comparison has 54,932 plus 18,128.
The old paired construction is a feasible special case. Its numerical embedding
matches every enlarged coefficient equation to 4.29e-10. All retained cross terms
remain available. No higher MPS moments were filled in or treated as known.

After the unchanged-span test, an additional fixed-number issue was identified:
linear words had been removed using identities available in the full cubic space,
but their replacements were 94.83–97.59% outside the compact spans in relative
Euclidean distance. An exact completion appended 16 prescribed number-dressed
linear directions and their adjoints. This increases the singlet optimization to
96,904 entries, or 115,032 including the nonsinglet proof. It retains every old
column. The additions follow the number identity; no larger-rank ranking search
or full-proof teacher chooses them.

## What the numerical searches establish

The direct 1.4 mHa target did not survive its residual allowance: its complete
replay is 3.558913 mHa. A number-projection repair retained all contraction-free
remainder and used a separately verified residual wedge witness. Its raw target
was never presented as a certified energy interval.

Early optimization runs had poorly conditioned linear solves. A measured
preconditioner reduced the fixed-target solves from repeatedly hitting 120
iterations to usually converging in one. Nevertheless, its coefficient residual
leveled off. Additional bound-maximization and the prescribed linear completion
were tested. The table contains only accepted full intervals; other saved
iterates remain untrusted proposals, even when a numerical process exited normally.

There is **no exact dual obstruction**. A diagnostic dual had negative moment
eigenvalues and an unbounded-for-this-purpose functional coordinate around 29.77;
it fails the accepting conditions for a residual-compatible family ceiling. The
finite runs therefore do not prove that these spans, much less all compact
representations, can never reach the target.

The projected lift coefficient provides an additional diagnostic of whether the
search actually uses the restored freedom:

| Construction | Coefficient norm of V | Coefficient norm of Z |
|---|---:|---:|
| Strong inherited proof | 0.1318644498 | 0.0000141065 |
{norm_rows}

These are exact rational squared norms, displayed after taking a floating square
root. They are not energy errors or an objective to match: a different successful
certificate need not resemble the strong proof. Every nonzero Z is paid by the
complete replay. The norm comparison is only a diagnostic of the returned factors.

In the original compact spans, the returned V norm is about 0.00000223. After
adding the prescribed number-derived linear directions, it rises to about
0.00126809. Thus the additional directions do let the returned construction use
more of this freedom. The complete energy interval, rather than that increase,
determines whether the change is useful.

## Costs and dependencies

- {len(runs)} completed instrumented processes total **{wall:.3f} process-wall seconds**,
  with {cpu:.3f} recorded CPU seconds. This includes failed and unsuccessful
  attempts, preparations, tests, repairs, and exact replays. It excludes editing,
  reading, interactive diagnostic probes, final inventory, and prior discovery;
  it is not a fresh-problem end-to-end timing claim.
- Peak recorded process memory is **{audit['peak_process_RSS_bytes']/1e9:.3f} GB**.
  Each dense normal-matrix cache occupies 582,496,712 bytes. It is a numerical
  construction cost on coefficient equations, not a determinant matrix or a
  small accepting certificate. Two different span caches are retained.
- Fresh full-degree operator preparation took 178.657 process seconds. It
  visits 1,206,592 cubic dictionary word pairs and stores 2,009,792 sparse-map
  nonzeros. The exact spin reduction retains 8,533 independent equations from
  28,461 full coefficient rows. There is no N-particle determinant enumeration.
- The compact-span seed inherits earlier H8 searches totaling 380.135 internal
  seconds, plus direct quartic-guide preparation and the nonsinglet construction.
  The upper MPS inherits its 414.36-second recorded discovery path. Historical
  preparation and earlier failed attempts are not made free by using frozen files.
- The old full-cubic search's 721.563-second historical discovery is a dependency
  of the inherited strong baseline, not of the new constructor. Its diagnostic
  replay in this pass is separately charged. These timings are not a matched
  hardware speed comparison.

The first quotient repair took 177.721 seconds. Proposal-side spin projection
was then accelerated using exact canonical spin flips and tested against the
original projector. The accepting verifier was left unchanged.

## Validation and preserved milestones

The focused suite passes **15 tests**, including all 325 pair-matrix units at
m=5,6, sparse m=12,16 cases, the nonzero zero-contraction counterexample, the
signed-factor norm test, independent six-mode bit action, a nonzero completion
through the unchanged accepting checker, all molecular adjoint-pair maps, all
662 six-mode matrix units for spin projection, and the number-dressed linear
identities. The 64 enumerated bit inputs belong only to an isolated algebra test.

All **{len(inherited):,} inherited files** and both supplied request files match
their original hashes. The original H6 result and fresh-geometry 1.034056 mHa
milestone remain preserved; the latter is 1.034064 mHa when its actual MPS upper
is independently rechecked with outward rounding. This pass uses the same frozen
H8 upper as its paired comparison.

The mathematical construction and its exact scope are in [MATHEMATICS.md]({SOURCE/'MATHEMATICS.md'}).
Machine-readable results are in [audit.json]({OUT/'audit.json'}), with all process
receipts retained in the runs directory. The [best new certificate]({ROOT/best['certificate']})
has a separate [complete replay]({ROOT/best['interval']}).

The useful result of this pass is a tested constructor that restores the missing
fixed-number freedom, together with a measured molecular outcome. Further claims
of compact H8 accuracy or a limitation theorem must follow accepting certificates.
'''
    (SOURCE/'REPORT.md').write_text(text)
    files = {}
    for root in (SOURCE, OUT):
        for p in sorted(root.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and p.name not in ('manifest.json', 'seal.json'):
                files[str(p.relative_to(ROOT))] = {'bytes': p.stat().st_size, 'sha256': sha(p)}
    manifest = {'kind': 'sealed_sector_quotient_campaign_v1', 'created_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'files': files, 'file_count': len(files), 'total_bytes': sum(v['bytes'] for v in files.values()),
                'parent_manifest': 'results/collective_completion_20260914/manifest.json',
                'parent_sha256': requests['parent_manifest_sha256'], 'inherited_files_verified': len(inherited),
                'objective_status': {'fixed_number_freedom_restored_and_tested': True,
                                     'new_compact_H8_1p6mHa': best['target_met'], 'general_solution': False}}
    dump(OUT/'manifest.json', manifest)
    errors = [p for p, rec in files.items() if sha(ROOT/p) != rec['sha256']]
    if errors: raise ValueError(('Seal verification failed', errors))
    seal = {'manifest_sha256': sha(OUT/'manifest.json'), 'all_new_file_hashes_verified': True,
            'all_inherited_file_hashes_verified': True, 'files': len(files), 'bytes': manifest['total_bytes']}
    dump(OUT/'seal.json', seal)
    print(json.dumps({'best_new': best, 'old_compact_width_mHa': old_compact['width_mHa'],
                      'recorded_process_seconds': wall, 'seal': seal}, indent=2))


if __name__ == '__main__': run()
