"""Record the changed pricing costs without rewriting the historical control."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]


def run(out):
    previous = ROOT/'results/spin_enrichment_20260913'; manifest_path = previous/'manifest.json'
    old_files = json.loads(manifest_path.read_text())['files']
    changed = [name for name, r in old_files.items() if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != r['sha256']]
    (out/'preservation_after.json').write_text(json.dumps({'checked': len(old_files), 'changed': changed,
        'manifest_sha256': hashlib.sha256(manifest_path.read_bytes()).hexdigest()}, indent=2)+'\n')
    if changed:
        raise ValueError('Historical evidence changed: '+repr(changed))
    campaign = json.loads((out/'campaign/summary.json').read_text())
    summary = json.loads((out/'campaign/continue/summary.json').read_text())
    old_ledger = json.loads((previous/'cost_ledger.json').read_text()); trials = []
    for trial in summary['trials']:
        row = dict(trial)
        if trial['accepted']:
            directory = out/'campaign/continue'/trial['case']
            receipt = json.loads((directory/'receipt.json').read_text()); cert = json.loads((directory/'certificate.json').read_text())
            row.update({k: receipt[k] for k in ('Gram_entries', 'used_map_nonzeros', 'constructed_map_nonzeros',
                'construction_seconds', 'solve_seconds', 'export_seconds', 'accept_seconds', 'wall_seconds',
                'certificate_bytes', 'anti_dimensions', 'status')})
            row['expanded_certificate_bytes'] = receipt['accepted']['expanded_certificate_bytes']
            row['direction_nonzero_coefficients'] = sum(c != 0 for b in cert['anti_blocks'] for v in b['directions'] for c in v)
            row['factor_nonzero_coefficients'] = sum(c != 0 for b in cert['anti_blocks']+cert['base_blocks'] for v in b['factor'] for c in v)
            row['exact_residual_penalty_mHa'] = float(F(receipt['accepted']['retained']['residual_l1'])*1000)
        trials.append(row)
    incremental = campaign['arms'][0]['wall_seconds']; ancestry = old_ledger['prior_causal_descriptor_discovery_seconds']
    result = {'campaign': campaign, 'trials': trials, 'trace_construction': summary['model_construction'],
        'prior_causal_source_discovery_seconds': ancestry, 'new_trace_incremental_discovery_seconds': incremental,
        'new_trace_cumulative_causal_discovery_seconds': ancestry+incremental,
        'frozen_coefficient_incremental_discovery_seconds': old_ledger['incremental_continuation_seconds'],
        'frozen_coefficient_cumulative_causal_discovery_seconds': old_ledger['cumulative_causal_discovery_seconds'],
        'scope': 'Historical, single-run, equal-budget comparison. Source discovery and previous failed controls remain charged in their ledgers; no new full-family SDP or ablation is counted.'}
    if (out/'fresh_replay.json').exists():
        result['fresh_replay_and_comparison_seconds'] = json.loads((out/'fresh_replay.json').read_text())['comparison_wall_seconds']
    if (out/'full_dual_watchdog.json').exists():
        result['full_dual_diagnostic'] = json.loads((out/'full_dual_watchdog.json').read_text())
        result['verification_execution'] = 'Fresh energy replay/comparison and the conditional full-dual diagnostic ran concurrently after discovery. Their individual wall times remain charged but must not be added as calendar elapsed time. Neither overlapped the timed discovery campaign.'
    if (out/'audit.json').exists():
        result['fresh_dual_and_physical_audit_seconds'] = json.loads((out/'audit.json').read_text())['wall_seconds']
    (out/'cost_ledger.json').write_text(json.dumps(result, indent=2)+'\n')
    (out/'environment.json').write_bytes((previous/'environment.json').read_bytes())
    files = {ROOT/p for p in old_files}; files.add(manifest_path)
    files.update(p for p in (ROOT/'research/trace_pricing_20260913').rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    files.update(p for p in out.rglob('*') if p.is_file() and p.name != 'manifest.json')
    manifest = {str(p.resolve().relative_to(ROOT)): {'bytes': p.stat().st_size,
        'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)}
    (out/'manifest.json').write_text(json.dumps({'files': manifest}, indent=2)+'\n'); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); run(args.out.resolve())
