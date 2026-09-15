"""Preserve prior evidence and record measured costs of this continuation."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import platform

ROOT = Path(__file__).resolve().parents[2]


def run(out):
    prior_manifest = ROOT/'results/spin_subspace_20260913/manifest.json'
    prior = json.loads(prior_manifest.read_text())['files']
    changed = [name for name, row in prior.items() if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != row['sha256']]
    preservation = {'manifest_sha256': hashlib.sha256(prior_manifest.read_bytes()).hexdigest(),
        'checked': len(prior), 'changed': changed}
    (out/'preservation_after.json').write_text(json.dumps(preservation, indent=2)+'\n')
    if changed:
        raise ValueError('Previous evidence changed: '+repr(changed))
    campaign = json.loads((out/'campaign/summary.json').read_text())
    diagnostic = json.loads((out/'diagnostic_watchdog.json').read_text())
    prior_campaign = json.loads((ROOT/'results/spin_subspace_20260913/campaign/summary.json').read_text())
    source_seconds = next(a['wall_seconds'] for a in prior_campaign['arms'] if a['arm'] == 'full')
    prior_dual_seconds = json.loads((ROOT/'results/spin_subspace_20260913/full_dual_watchdog.json').read_text())['wall_seconds']
    ancestry = source_seconds+prior_dual_seconds
    raw_control = json.loads((out/'raw_full_watchdog.json').read_text()) if (out/'raw_full_watchdog.json').exists() else None
    cases = []
    for arm in campaign['arms']+([raw_control] if raw_control is not None else []):
        directory = out/'raw_full' if arm['arm'] == 'full_raw' else out/'campaign'/arm['arm']
        summary = json.loads((directory/'summary.json').read_text())
        for row in summary['trials']:
            record = {'arm': arm['arm'], **row}
            if row['accepted']:
                receipt = json.loads((directory/row['case']/'receipt.json').read_text())
                cert = json.loads((directory/row['case']/'certificate.json').read_text())
                record.update({key: receipt[key] for key in ('Gram_entries', 'constructed_map_nonzeros',
                    'used_map_nonzeros', 'construction_seconds', 'solve_seconds', 'export_seconds',
                    'accept_seconds', 'wall_seconds', 'certificate_bytes', 'anti_dimensions', 'status')})
                record['expanded_certificate_bytes'] = receipt['accepted']['expanded_certificate_bytes']
                record['direction_nonzero_coefficients'] = sum(c != 0 for b in cert['anti_blocks'] for row in b['directions'] for c in row)
                record['factor_nonzero_coefficients'] = sum(c != 0 for b in cert['base_blocks']+cert['anti_blocks'] for row in b['factor'] for c in row)
                record['exact_residual_penalty_mHa'] = float(F(receipt['accepted']['retained']['residual_l1'])*1000)
            cases.append(record)
    adaptive_seconds = next(a['wall_seconds'] for a in campaign['arms'] if a['arm'] == 'adaptive')
    result = {'diagnostic': diagnostic, 'campaign': campaign, 'raw_full_followup': raw_control, 'trials': cases,
        'prior_counterexample_ancestry': {'full_frame_source_solve_seconds': source_seconds,
            'dual_proposal_and_acceptance_seconds': prior_dual_seconds, 'total_seconds': ancestry},
        'adaptive_seconds_with_new_diagnostic': adaptive_seconds+diagnostic['wall_seconds'],
        'adaptive_seconds_with_counterexample_ancestry': adaptive_seconds+diagnostic['wall_seconds']+ancestry,
        'fresh_lower_discovery_and_controls_seconds': campaign['wall_seconds']+diagnostic['wall_seconds']+(raw_control['wall_seconds'] if raw_control else 0),
        'scope': 'Single-run observations. Shared Hamiltonian construction and old upper discovery are frozen inputs, not newly discovered here. Controls that reuse directions inherit all their discovery cost.'}
    if (out/'fresh_replay.json').exists():
        result['fresh_replay_seconds'] = json.loads((out/'fresh_replay.json').read_text())['wall_seconds']
    if (out/'full_dual_watchdog.json').exists():
        result['new_full_dual'] = json.loads((out/'full_dual_watchdog.json').read_text())
    (out/'cost_ledger.json').write_text(json.dumps(result, indent=2)+'\n')
    versions = {'python': platform.python_version(), 'platform': platform.platform(),
        'BLAS_thread_environment': 1, 'Clarabel_thread_setting': 'unchanged solver default'}
    from importlib.metadata import version
    versions.update({name: version(name) for name in ('numpy', 'scipy', 'cvxpy', 'clarabel')})
    (out/'environment.json').write_text(json.dumps(versions, indent=2)+'\n')
    files = {ROOT/name for name in prior}
    files.add(prior_manifest)
    files.update(p for p in (ROOT/'research/spin_completion_20260913').rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    files.update(p for p in out.rglob('*') if p.is_file() and p.name != 'manifest.json')
    manifest = {str(p.resolve().relative_to(ROOT)): {'bytes': p.stat().st_size,
        'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)}
    (out/'manifest.json').write_text(json.dumps({'files': manifest}, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); run(args.out.resolve())
