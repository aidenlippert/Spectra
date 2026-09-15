"""Bind new selected-CI upper states to existing lower proofs and replay both."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_symbolic import decode
from research.certificate_scaling.cubic_interval_replay import run


def prepare():
    scoreboard = json.loads((ROOT/'results/certificate_scaling/active_space_ladder_latest.json').read_text())
    tasks = []
    for oldrow in scoreboard['rows']:
        n = int(oldrow['system'][1:])
        folder = ROOT/f'results/all_angles_20260913/selected_refinement/h{n}'
        if n == 10:
            folder = ROOT/'results/all_angles_20260913/selected_refinement/h10_extended'
        new = json.loads((folder/'receipt.json').read_text())
        old = json.loads((ROOT/oldrow['interval']).read_text())
        passing = [row for row in new['rows'] if Fraction(row['upper'])-Fraction(old['lower']) <= Fraction(16,10000)]
        if not passing:
            raise ValueError(f'No passing new upper for H{n}')
        chosen = passing[0]
        certificate = ROOT/old['certificate']
        raw = certificate.read_bytes(); cert = json.loads(raw)
        fixture_path = ROOT/f'results/certificate_scaling/active_space_ladder/h{n}/fixture.json'
        fxraw = fixture_path.read_bytes(); fixture = json.loads(fxraw)
        if sha256(raw).hexdigest() != old['certificate_sha256']:
            raise ValueError('Existing lower certificate content changed')
        if sha256(fxraw).hexdigest() != new['fixture_sha256']:
            raise ValueError('New upper fixture mismatch')
        if any(cert[k] != fixture[k] for k in ('modes','particles')):
            raise ValueError('Lower and upper sector mismatch')
        if decode(cert['hamiltonian'],cert['modes'],4) != decode(fixture['hamiltonian'],fixture['modes'],4):
            raise ValueError('Lower and upper Hamiltonian mismatch')
        ref = folder/chosen['witness_file']
        if sha256(ref.read_bytes()).hexdigest() != chosen['witness_sha256']:
            raise ValueError('Upper witness content changed')
        tasks.append({'system': f'H{n}', 'method': old['method'],
                      'certificate': str(certificate.relative_to(ROOT)),
                      'proof': old['proof'], 'reference': str(ref.relative_to(ROOT)),
                      'fixture': str(fixture_path.relative_to(ROOT)),
                      'fixture_sha256': sha256(fxraw).hexdigest(),
                      'selected_basis': chosen['basis_size'],
                      'upper_discovery_through_pass_seconds': chosen['elapsed_seconds'],
                      'all_run_seconds': new['wall_seconds'],
                      'upper_method': new['kind'], 'old_interval': oldrow['interval']})
    return tasks


def replay(task, out):
    result = run(ROOT/task['certificate'], ROOT/task['reference'], out,
                 task['method'], ROOT/task['proof'] if task['proof'] else None)
    # The shared replayer retains an FCI-specific historical scope string.
    # This run instead uses an explicitly sourced selected-CI witness.
    result['scope'] = ('Exact two-sided interval for original frozen finite Hamiltonian. '
                       'Lower certificate reused; new upper discovered by selected CI from H and supplied HF determinant. '
                       'No saved FCI coefficients used by upper discovery. No asymptotic accuracy or physical-model guarantee.')
    result['upper_provenance'] = task
    out.write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--tasks', type=Path, required=True)
    parser.add_argument('--system')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if args.prepare:
        args.tasks.write_text(json.dumps(prepare(), indent=2)+'\n')
    else:
        task = next(x for x in json.loads(args.tasks.read_text()) if x['system'] == args.system)
        replay(task, args.out)
