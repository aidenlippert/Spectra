"""Measured incremental costs, causal ancestry and preserved evidence."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]


def run(out):
    previous = ROOT/'results/spin_completion_20260913'
    prior_manifest = previous/'manifest.json'; prior = json.loads(prior_manifest.read_text())['files']
    changed = [p for p, r in prior.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != r['sha256']]
    (out/'preservation_after.json').write_text(json.dumps({'checked': len(prior), 'changed': changed,
        'manifest_sha256': hashlib.sha256(prior_manifest.read_bytes()).hexdigest()}, indent=2)+'\n')
    if changed:
        raise ValueError('Prior evidence changed: '+repr(changed))
    campaign = json.loads((out/'campaign/summary.json').read_text()); old = json.loads((previous/'cost_ledger.json').read_text())
    rows = []
    for arm in campaign['arms']:
        directory = out/'campaign'/arm['arm']; summary = json.loads((directory/'summary.json').read_text())
        for trial in summary['trials']:
            row = {'arm': arm['arm'], **trial}
            if trial['accepted']:
                receipt = json.loads((directory/trial['case']/'receipt.json').read_text())
                cert = json.loads((directory/trial['case']/'certificate.json').read_text())
                row.update({k: receipt[k] for k in ('Gram_entries', 'used_map_nonzeros', 'constructed_map_nonzeros',
                    'construction_seconds', 'solve_seconds', 'export_seconds', 'accept_seconds', 'wall_seconds',
                    'certificate_bytes', 'anti_dimensions', 'status')})
                row['expanded_certificate_bytes'] = receipt['accepted']['expanded_certificate_bytes']
                row['direction_nonzero_coefficients'] = sum(c != 0 for b in cert['anti_blocks'] for v in b['directions'] for c in v)
                row['factor_nonzero_coefficients'] = sum(c != 0 for b in cert['anti_blocks']+cert['base_blocks'] for v in b['factor'] for c in v)
                row['exact_residual_penalty_mHa'] = float(F(receipt['accepted']['retained']['residual_l1'])*1000)
            rows.append(row)
    incremental = next(a['wall_seconds'] for a in campaign['arms'] if a['arm'] == 'continue')
    result = {'campaign': campaign, 'trials': rows,
        'prior_causal_descriptor_discovery_seconds': old['adaptive_seconds_with_counterexample_ancestry'],
        'incremental_continuation_seconds': incremental,
        'cumulative_causal_discovery_seconds': incremental+old['adaptive_seconds_with_counterexample_ancestry'],
        'preceding_whole_pass_seconds_excluding_ancestry_and_validation': old['fresh_lower_discovery_and_controls_seconds'],
        'preceding_validation_seconds': old['fresh_replay_seconds'],
        'scope': 'Current continuation and ablation are additional costs. Previous failed full controls remain in the preceding whole-pass ledger. Single-run measurements; BLAS environment one, Clarabel automatic threads.'}
    if (out/'fresh_replay.json').exists():
        result['fresh_replay_seconds'] = json.loads((out/'fresh_replay.json').read_text())['wall_seconds']
    if (out/'full_dual_watchdog.json').exists():
        result['full_dual_diagnostic'] = json.loads((out/'full_dual_watchdog.json').read_text())
        result['verification_overlap'] = 'The fresh interval replays and the separate dual diagnostic overlapped; their wall costs are reported separately, not summed as elapsed calendar time.'
    if (out/'audit.json').exists():
        result['fresh_dual_and_physical_audit_seconds'] = json.loads((out/'audit.json').read_text())['wall_seconds']
    (out/'cost_ledger.json').write_text(json.dumps(result, indent=2)+'\n')
    (out/'environment.json').write_bytes((previous/'environment.json').read_bytes())
    files = {ROOT/p for p in prior}; files.add(prior_manifest)
    files.update(p for p in (ROOT/'research/spin_enrichment_20260913').rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    files.update(p for p in out.rglob('*') if p.is_file() and p.name != 'manifest.json')
    manifest = {str(p.resolve().relative_to(ROOT)): {'bytes': p.stat().st_size,
        'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)}
    (out/'manifest.json').write_text(json.dumps({'files': manifest}, indent=2)+'\n'); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); run(args.out.resolve())
