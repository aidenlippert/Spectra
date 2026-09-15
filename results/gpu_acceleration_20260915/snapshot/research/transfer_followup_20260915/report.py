"""Summarize all recorded attempts and preserve the previous milestone chain."""
import datetime
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import subprocess
from research.transfer_solver_20260915.budget import ROOT, OUT, dump


def read(path, default=None):
    return json.loads(path.read_text()) if path.exists() else default


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        while block := stream.read(1 << 20):
            h.update(block)
    return h.hexdigest()


def run():
    if not (OUT/'followup_execution.json').exists():
        raise ValueError('Finish the declared follow-up work before reporting completion')
    if not (OUT/'independent_h4_oracle.json').exists() or not (OUT/'fragments/summary.json').exists():
        raise ValueError('Resolve the independent-oracle and composition checks before issuing this report')
    attempts = [read(p) for p in sorted((OUT/'runs').glob('*.json'))]
    byname = {row['name']: row for row in attempts}
    cases = []
    for name in ('h4_control', 'h8_cold', 'h6_asymmetric', 'water_asymmetric', 'h10_size'):
        p = OUT/'cases'/name
        interval = read(p/'interval.json')
        execution = read(p/'execution.json', {})
        stages = [byname[name+'_'+step] for step in ('generate', 'state', 'upper', 'nonsinglet', 'prepare', 'stage1', 'stage2', 'replay') if name+'_'+step in byname]
        cases.append({'case': name, 'specification': read(p/'specification.json'), 'execution': execution,
                      'interval': interval, 'forecast': read(p/'prepared/forecast.json'),
                      'fci_reference': read(p/'fci_reference.json'), 'peak_RSS_bytes': max((x['peak_child_RSS_bytes'] for x in stages), default=0),
                      'stage_process_wall_seconds': sum(x['wall_seconds'] for x in stages), 'stages': stages})
    baseline = []
    for c in cases:
        name = c['case']
        p = OUT/'cases'/name/'enumerated_baseline'
        receipt = read(p/'replay.json')
        if receipt:
            names = [name+'_enumerated_construct', name+'_enumerated_verify']
            runtime = sum(byname[n]['wall_seconds'] for n in names)+byname[name+'_generate']['wall_seconds']
            baseline.append({'case': name, 'receipt': receipt, 'seconds_including_shared_integral_generation': runtime,
                             'shared_integrals_charged_once_to_each_method_comparison': True})
    controls = []
    for p in sorted((OUT/'controls').glob('*')):
        source_name = 'h8_cold' if p.name.startswith('h8') else 'h4_control'
        control_steps = [x for x in attempts if x['name'].startswith(p.name+'_')]
        source_steps = [byname[source_name+'_'+step] for step in ('prepare', 'stage1', 'stage2', 'replay')]
        controls.append({'control': p.name, 'execution': read(p/'execution.json'), 'interval': read(p/'interval.json'),
            'forecast': read(p/'prepared/forecast.json'), 'control_process_seconds': sum(x['wall_seconds'] for x in control_steps),
            'source_same_stage_seconds': sum(x['wall_seconds'] for x in source_steps),
            'source_case': source_name, 'shared_state_integrals_upper_nonsinglet_excluded_from_both': True})
    physical = []
    for basis in ('cc-pvdz', 'cc-pvtz'):
        s = read(OUT/f'ch2_model_study/{basis}_S0/result.json')
        t = read(OUT/f'ch2_model_study/{basis}_S1/result.json')
        if s and t:
            physical.append({'basis': basis, 'singlet': s, 'triplet': t,
                'CASSCF_electronic_gap_kcal_mol': (s['CASSCF_total_Ha']-t['CASSCF_total_Ha'])*627.509474,
                'SC_NEVPT2_electronic_gap_kcal_mol': (s['SC_NEVPT2_total_Ha']-t['SC_NEVPT2_total_Ha'])*627.509474,
                'gap_convention': 'E_S-E_T', 'rigorous_physical_error_bar': False})
    fail = [x for x in attempts if x['status'] != 'passed']
    ledger = {'instrumented_attempts': attempts, 'count': len(attempts),
        'sum_process_wall_seconds': sum(x['wall_seconds'] for x in attempts),
        'sum_process_CPU_seconds': sum(x['child_user_seconds']+x['child_system_seconds'] for x in attempts),
        'peak_single_process_RSS_bytes': max(x['peak_child_RSS_bytes'] for x in attempts),
        'failed_attempts': fail, 'cases': cases, 'controls': controls, 'enumerated_baselines': baseline,
        'scope': 'All instrumented computation in this pass, including controls and failures. Editing, browsing, report generation and preservation I/O are not included. Earlier research costs remain in their earlier ledgers.',
        'physical_model_study': physical,
        'hardware': subprocess.run(['/usr/sbin/sysctl', '-n', 'hw.model', 'hw.memsize', 'hw.ncpu', 'machdep.cpu.brand_string'], capture_output=True, text=True).stdout.splitlines(),
        'platform': platform.platform(), 'one_heavy_process_at_a_time': True}
    dump(OUT/'cost_ledger.json', ledger)
    parent_path = ROOT/'results/h8_discovery_repeat_20260915/manifest.json'
    expected_parent = '0855be3e3d8c900cf3ab907106cc1526df0f465041e88daa1e56b8bd88fc8d81'
    if sha(parent_path) != expected_parent:
        raise ValueError('Prior conditional-repeat manifest changed')
    protected = read(ROOT/'results/h8_spin_import_20260914/preservation_before.json')['files']
    for parent in ('h8_spin_import_20260914', 'h8_discovery_repeat_20260915'):
        p = ROOT/'results'/parent/'manifest.json'
        protected.update(read(p)['files'])
        protected[str(p.relative_to(ROOT))] = {'sha256': sha(p)}
    changed = [name for name, row in protected.items() if not (ROOT/name).is_file() or sha(ROOT/name) != row['sha256']]
    if changed:
        raise ValueError(('Previously protected files changed', changed[:8]))
    protocol = read(OUT/'validation_set_protocol.json')
    if any(sha(ROOT/p) != h for p, h in protocol['source_hashes'].items()):
        raise ValueError('Frozen constructor source no longer matches')
    dump(OUT/'preservation.json', {'previous_files_verified': len(protected), 'changed_previous_files': [],
                                  'frozen_constructor_files_verified': len(protocol['source_hashes'])})
    table = []
    for c in cases:
        value = f"{c['interval']['width_mHa']:.9f}" if c['interval'] else 'No complete interval'
        target = ('yes' if c['interval']['target_met'] else 'no') if c['interval'] else 'not established'
        table.append(f"| {c['case']} | {value} | {target} | {c['execution'].get('elapsed_seconds',0):.2f} | {c['peak_RSS_bytes']/1e6:.1f} | {c['execution'].get('last_step','unknown')} |")
    comparison_rows = []
    for c in controls:
        interval = c['interval']
        width = f"{interval['width_mHa']:.9f}" if interval else 'unavailable'
        forecast = c['forecast'] or {}
        comparison_rows.append(f"| {c['control']} | {width} | {forecast.get('retained_Gram_entries','?')} | {c['control_process_seconds']:.2f} | {c['source_same_stage_seconds']:.2f} |")
    baseline_rows = [f"| {c['case']} | {c['receipt']['width_mHa']:.9f} | {c['seconds_including_shared_integral_generation']:.3f} | {c['receipt']['cost']['full_sector_labels_enumerated']} |" for c in baseline]
    fci_rows = []
    for c in cases:
        f = c['fci_reference']
        if f:
            fci_rows.append(f"| {c['case']} | {f['numerical_solve_seconds']:.4f} | {f['converged']} | {f.get('upper_minus_reference_mHa',float('nan')):.6f} | {f.get('reference_minus_lower_mHa',float('nan')):.6f} |")
    physical_rows = [f"| {x['basis']} | {x['CASSCF_electronic_gap_kcal_mol']:.5f} | {x['SC_NEVPT2_electronic_gap_kcal_mol']:.5f} |" for x in physical]
    fail_lines = [f"- {x['name']}: {x['status']}, {x['wall_seconds']:.3f} s; log in `runs/{x['name']}.log`." for x in fail]
    report = f'''# Cold construction, transfer, and measured controls

The prior **0.767448135463 mHa H8 milestone remains achieved and preserved**.
The new input-derived constructor produced a fresh H8 certificate of
**0.605110981693 mHa**, including fresh integral generation, MPS discovery,
nonsinglet construction, coefficient preparation, optimization and exact replay.
The original H8 rational Hamiltonian and the new one are exactly unitarily
equivalent under two paired orbital sign flips. This is a complete cold run,
not a repeat conditional on an old checkpoint.

## Frozen cases and complete cost

The changed-input cases were frozen together before their outcomes were known.
H4 was the development control. A frozen procedure adapts its algebra to the
input; it does not reuse H8 dimensions or parity labels.

| Case | Certified total width (mHa) | <=1.6 mHa | Complete campaign seconds | Peak process MB | Last step |
|---|---:|---|---:|---:|---|
{chr(10).join(table)}

All completed certificates cover the entire fixed-N sector, using separate S=0
and M_S=1 proofs plus one exact spin-defect allowance. Their upper Rayleigh
quotients and lower remainders use integer/rational arithmetic, without numerical
packages on the accepting path. The asymmetric H6 case has only the trivial
total-parity class; water changes the molecule and electron/orbital counts.
These are a few deliberately selected cases, not a statistical reliability test.

The per-case times start with new integrals. Historical state, coefficient maps,
optimization checkpoints, and nonsinglet certificates are absent from their
constructor inputs. The fresh MPS supplies moments to discovery. Numerical FCI
references below were run only after discovery and never supplied as teachers.

The H8 search retains **343,424 singlet Gram entries**, with largest block 260,
versus 335,168 and 256 in the earlier winner. The declared difference is that
mixed linear words are retained uniformly at every size. This small enlargement
means the fresh run is not an identical-cone reproduction of the previous winner.
Fewer exported factors are not used as a proxy for discovery cost.

## Equivalent and restricted representation controls

| Control | Exact width (mHa) | Singlet Gram entries | Control prep/solve/replay seconds | Same stages in full highest-weight case |
|---|---:|---:|---:|---:|
{chr(10).join(comparison_rows)}

The magnetic controls use the same attainable spin-averaged cone, physical zero
initial Gram matrices, new MPS functional, optimizer and stopping rule. Shared
integrals, state, upper and nonsinglet costs are excluded from both sides of
these representation-only comparisons. The restricted64 control instead retains
at most the first 64 deterministic rational columns per mixed highest-weight
block. It tests a declared smaller family; it is not the best possible restriction.
No failed finite run is promoted to an exact family obstruction.

The mathematical cone equivalence, exact kernel checks and acceptance rule are
in [DERIVATION.md]({ROOT}/research/transfer_solver_20260915/DERIVATION.md).
Symmetry reduction itself is established; see
[Gatermann and Parrilo](https://arxiv.org/abs/math/0211450).

## Comparisons at the same target and with numerical FCI

The enumerated baseline constructs its own integer upper and rational
Cholesky-remainder lower. It covers all fixed-N determinants, partitioned only
by alpha count, and independently replays the resulting certificate. It is
capped at 2,000 determinant labels. Its time includes the common integral
generation cost once for each method comparison.

| Case | Enumerated exact width (mHa) | Generation + construction + replay seconds | Enumerated labels |
|---|---:|---:|---:|
{chr(10).join(baseline_rows)}

Both methods target <=1.6 mHa; the table states their actual attained widths.
This is a small rigorous reference implementation, not a survey of the strongest
available certified solvers. It can directly reveal when eliminating enumeration
has not yet paid off. No overall competitive advantage follows from this campaign.

| Case | Numerical FCI solve seconds | Converged | Spectra upper minus FCI (mHa) | FCI minus Spectra lower (mHa) |
|---|---:|---|---:|---:|
{chr(10).join(fci_rows)}

FCI uses fixed tight numerical convergence, not a rigorous lower guarantee.
These energies are for the floating integrals before rational rounding; the
recorded rounding envelope and numerical comparison slack are explicit in the
receipts. Its timings and numerical agreement must not be presented as exact
error bars. See the [PySCF FCI documentation](https://pyscf.org/user/ci.html).

## Independent checks and fragment composition

`independent_h4_oracle.json` records a separately implemented full bit-state
Hamiltonian, rational LDL positivity test and explicit integer MPS expansion.
It imports neither the original CAR nor tensor checker. Its enumeration is
charged as a validation control, not hidden in construction. The raised-lower
mutation must fail. This strengthens validation on a small system; it is not
a formal verification of every implementation path.

`fragments/summary.json` records exact composition of 2, 4 and 8 independent
H4 fragments with fixed local electron counts, including direct exact product-MPS
uppers. Removing those local constraints or adding cross-fragment hopping is
refused. This establishes the compositional rule, not size consistency of a
monolithic restricted optimizer. The latter remains a separate question.

## CH2 physical-model control

The existing exact fixed-geometry STO-3G CAS(6,6) gap near 25.49 kcal/mol remains
a certificate for that model, not a reliable experimental prediction. The new
small control optimizes the A1 singlet and B1 triplet CASSCF(6,6) geometries,
checks total spin at every evaluation, and adds SC-NEVPT2 single-point corrections.

| Basis | CASSCF E_S-E_T (kcal/mol) | SC-NEVPT2 E_S-E_T (kcal/mol) |
|---|---:|---:|
{chr(10).join(physical_rows) if physical_rows else 'No completed pair; inspect the preserved failed jobs.'}

The experimental study reports 9.00 +/- 0.09 kcal/mol for the observed splitting
and an inferred electronic separation of 8.7 +/- 0.5 kcal/mol. The latter is the
appropriate rough electronic comparison here; no zero-point term was computed.
[Leopold, Murray, Miller and Lineberger](https://experts.umn.edu/en/publications/methylene-a-study-of-the-xsup3supbsub1sub-and-%C3%A3sup1supasub1sub-st/).
These new conventional calculations have no certified physical error bar.
Basis sensitivity is not an error enclosure, and agreement is retrospective.
The [CASSCF](https://pyscf.org/user/mcscf.html) and
[SC-NEVPT2](https://pyscf.org/user/mrpt.html) components are existing methods,
so the control does not establish a Spectra-specific predictive advantage.

## All attempts and preservation

The ledger contains **{len(attempts)} instrumented attempts**, summing to
**{ledger['sum_process_wall_seconds']:.3f} process-wall seconds** and
**{ledger['sum_process_CPU_seconds']:.3f} CPU seconds**. Maximum single-process
RSS was **{ledger['peak_single_process_RSS_bytes']/1e6:.3f} MB**. One heavy process
ran at a time. This includes unsuccessful attempts and the controls; it is not
the cost of one fresh solution. Editing, browsing, reporting and preservation I/O
are outside this compute ledger, and earlier R&D costs stay in their earlier
ledgers rather than being silently zeroed.

{chr(10).join(fail_lines) if fail_lines else 'Every instrumented job passed its process gate.'}

All **{len(protected):,} previously protected files** and all frozen constructor
source hashes match. No older achieved result has been overwritten.

Cold success and limited transfer are now evidenced. Broad reliability,
manageable larger-system scaling, the attainable-family ceiling, monolithic
fragment size consistency, competitive advantage over strong methods, and
prospective experimental value require their own evidence. A numerical miss
must be diagnosed using its upper, lower, residual and budget receipts rather
than being called a proof of impossibility.

- [Complete compute ledger]({OUT}/cost_ledger.json)
- [Frozen validation protocol]({OUT}/validation_set_protocol.json)
- [Exact original-H8 bridge]({OUT}/h8_original_model_bridge.json)
'''
    path = ROOT/'research/transfer_solver_20260915/REPORT.md'
    with path.open('x') as stream:
        stream.write(report)
    # Keep later analysis separate if a specific new result requires another pass.
    files = {}
    for folder in [OUT]+[ROOT/'research'/p for p in ('transfer_solver_20260915', 'transfer_followup_20260915', 'ch2_model_study_20260915')]:
        for p in sorted(folder.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and p.name not in ('manifest.json', 'seal.json'):
                files[str(p.relative_to(ROOT))] = {'bytes': p.stat().st_size, 'sha256': sha(p)}
    dump(OUT/'manifest.json', {'kind': 'sealed_cold_transfer_campaign_v1',
        'created_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'parent_manifest': str(parent_path.relative_to(ROOT)), 'parent_sha256': expected_parent,
        'files': files, 'file_count': len(files)})
    if any(sha(ROOT/p) != row['sha256'] for p, row in files.items()):
        raise ValueError('New pass files changed while sealing')
    dump(OUT/'seal.json', {'manifest_sha256': sha(OUT/'manifest.json'), 'new_files_verified': len(files),
                          'previous_files_verified': len(protected), 'all_hashes_match': True})
    print(path)


if __name__ == '__main__':
    run()
