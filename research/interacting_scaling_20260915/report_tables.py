"""Build compact report tables from accepted intervals and measured run receipts."""
from fractions import Fraction as F
import json
from pathlib import Path
from research.interacting_scaling_20260915.budget import ROOT, OUT
from research.interacting_scaling_20260915.account import belongs_to_cold


def read(path):
    return json.loads(path.read_text())


def link(path, label='receipt'):
    return f'[{label}](../../{path.relative_to(ROOT)})'


def width(interval):
    return f"{float(1000*(F(interval['upper_Ha'])-F(interval['lower_Ha']))):.9f}"


def build():
    ledger = read(OUT/'accounting.json')
    runs = [read(p) for p in sorted((OUT/'runs').glob('*.json'))]
    lines = ['# Measured campaign tables', '',
        'All energy widths below come from exact accepted endpoints. A process that',
        'finished successfully can still miss the 1.6 mHa accuracy target. Times are',
        'single measured runs on the local M1 host, without an estimated variance.', '',
        '| Fresh run | Width, mHa | Complete seconds | Target met | Selected Gram entries | Peak child RSS, MB |',
        '|---|---:|---:|---|---:|---:|']
    for name, d in ledger['cold_runs'].items():
        record = read(OUT/'cold'/name/'result.json')
        value = width(record['best']['interval']) if record.get('best') else 'No accepted interval'
        gram = f"{d['family']['Gram_entries']:,}" if d.get('family') else '—'
        lines.append(f"| {link(OUT/'cold'/name/'result.json', name)} | {value} | {d['complete_wall_seconds']:.3f} | {d['target_met']} | {gram} | {d['peak_child_RSS_bytes']/1e6:.1f} |")
    lines += ['', '## Representation and retained witnesses', '',
        'Selected Gram entries above describe the accepted representation. The sum',
        'below also counts earlier main-proof levels prepared during that cold run;',
        'it is not a peak allocation or arithmetic-work estimate. Magnetic-screen',
        'work is charged separately in the component ledger. Retained witness size',
        'includes the main and nonsinglet proof, MPS/upper and orbital transform,',
        'but excludes shared Hamiltonian descriptions, checker code and replay logs.', '',
        '| Run | Sum of prepared main Gram entries | Largest accepted block | Coefficient rows | Map nonzeros | Retained witness, MB | MPS nonzero integer entries |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for name, d in ledger['cold_runs'].items():
        f = d.get('family')
        if not f: continue
        witness = f.get('retained_witness_bytes_excluding_shared_input')
        size = f'{witness/1e6:.3f}' if witness is not None else '—'
        lines.append(f"| {name} | {d['sum_main_prepared_Gram_entries_including_rejected_levels']:,} | {f['largest_block']:,} | {f['coefficient_rows']:,} | {f['coefficient_nonzeros']:,} | {size} | {f.get('MPS_nonzero_integer_entries', '—')} |")
    lines += ['', '## Cold component costs', '',
        'These child-stage sums are included within the fresh clock above; they are',
        'not additional expenses. Rejected proof levels remain included.', '',
        '| Run | Input/orbitals, s | State/upper, s | Magnetic screen, s | Proof preparation, s | Proof solves, s | Exact replay, s | Other, s |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for name in ledger['cold_runs']:
        groups = dict.fromkeys(('input', 'state', 'magnetic', 'prepare', 'solve', 'replay', 'other'), 0.)
        for row in runs:
            if not belongs_to_cold(row['name'], name): continue
            suffix = row['name'][len(name)+1:]
            key = ('input' if suffix.startswith('input_') else 'state' if suffix.startswith('state_screen_')
                else 'magnetic' if '_magnetic_' in '_'+suffix else 'prepare' if suffix.endswith('_prepare')
                else 'solve' if suffix.endswith('_solve') else 'replay' if suffix.endswith('_replay') else 'other')
            groups[key] += row['wall_seconds']
        lines.append('| '+name+' | '+' | '.join(f'{v:.3f}' for v in groups.values())+' |')
    lines += ['', '## Coupling paths', '',
        'Intermediate Hamiltonians are diagnostic models. The endpoint at lambda=1',
        'includes the exact original-model rotation allowance. Each path includes',
        'state continuation and failed proof attempts, but source integral/orbital',
        'preparation is additional and remains in the campaign ledger.', '']
    for path in sorted((OUT/'coupling').glob('*/results.json')):
        d = read(path)
        if lines[-1] != '': lines.append('')
        lines += [f'### {path.parent.name}', '', f"Measured continuation clock: {d['wall_seconds']:.3f} seconds; {link(path)}.", '',
            '| Coupling | Width, mHa | Target met | Selected case |', '|---|---:|---|---|']
        for row in d['outcomes']:
            best = row.get('best')
            if best:
                interval, case = best['interval'], Path(best['case'])
            else:
                case = Path(row.get('case', row.get('base_case', '')))
                ip = case/'original_interval.json'
                if not ip.exists(): ip = case/'exact/interval.json'
                interval = read(ip) if ip.exists() else None
            value = width(interval) if interval else 'No accepted interval'
            target = bool(interval and interval['target_met'])
            label = link(case/'design.json', case.name) if (case/'design.json').exists() else '—'
            lines.append(f"| {row.get('lambda', row.get('coupling'))} | {value} | {target} | {label} |")
    lines += ['', '## H12 post hoc diagnostics', '',
        'These continuations do not replace the failed frozen H12 run. They reuse',
        'its state and selected maps; its complete 1,502.800-second clock remains',
        'additional. Parent-clock sums below exclude intervening work on other',
        'cases and are not a newly timed cold pipeline.', '']
    for filename in ('h12_recovery_result.json', 'h12_main_refinement_result.json'):
        path = OUT/filename
        if not path.exists(): continue
        d = read(path)
        value = width(d['interval'])+' mHa' if d.get('interval') else 'No accepted complete interval'
        seconds = d.get('additional_recovery_parent_seconds', d.get('additional_parent_seconds'))
        clock = f'{seconds:.3f} additional seconds' if seconds is not None else d.get('skipped', 'No measured execution')
        lines.append(f'- {link(path, filename)}: {value}; {clock}.')
        total = d.get('all_H12_diagnostic_parent_seconds', d.get('combined_diagnostic_parent_seconds'))
        if total is not None:
            lines.append(f'  Recorded diagnostic parent-clock sum through this point: {total:.3f} seconds; source cold run additional.')
    lines += ['', '## Disconnected controls', '',
        '| Control | Width, mHa | Local Fock labels over every charge | Global determinants enumerated |', '|---|---:|---:|---:|']
    for path in sorted((OUT/'fragment_controls').glob('*/replay.json')):
        d = read(path)
        lines.append(f"| {link(path, path.parent.name)} | {width(d)} | {d['local_labels_enumerated_including_all_charges']} | {d['global_determinant_labels_enumerated']} |")
    lines += ['', '## External algorithm controls', '',
        'Numerical FCI estimates and exact enumerated certificates are separate',
        'comparisons. FCI runs after proof discovery on the pre-rounding integral',
        'model; it does not supply a rigorous two-sided certificate.', '',
        '| Exact enumerated control | Fresh complete seconds | Accepted exact width, mHa | Literal original input match |',
        '|---|---:|---:|---|']
    for path in sorted((OUT/'models').glob('*_enumerated_control/baseline_result.json')):
        d = read(path)
        value = width(d['interval']) if d.get('interval') else 'No accepted interval'
        lines.append(f"| {link(path, path.parent.name)} | {d['fresh_complete_pipeline_seconds']:.3f} | {value} | {d['matched_original_input']} |")
    lines += ['', 'Numerical FCI timing receipts:', '']
    for row in runs:
        if row['name'].endswith('_numerical_fci_solve'):
            path = OUT/'runs'/(row['name']+'.json')
            lines.append(f"- {link(path, row['name'])}: {row['wall_seconds']:.3f} seconds, {row['status']}.")
    lines += ['', '## Complete campaign ledger', '',
        f"There were {ledger['campaign_attempt_count']} measured child attempts, totaling "
        f"{ledger['campaign_child_stage_seconds']:.3f} seconds, including failed attempts, tests, "
        'and dependency installation. This is a sum of measured compute stages, not one fresh solver run.', '',
        'Editorial work and routine file inspection are excluded from that compute ledger.',
        'Cold orchestration clocks are reported separately and must not be added again.',
        'No new cloud instances or external spending were used.', '',
        f"See {link(OUT/'accounting.json', 'the complete accounting record')} for exact receipt hashes, failed process attempts, "
        'accuracy misses and representation sizes.', '']
    path = ROOT/'research/interacting_scaling_20260915/TABLES.md'
    with path.open('x') as stream: stream.write('\n'.join(lines))
    print(path)


if __name__ == '__main__': build()
