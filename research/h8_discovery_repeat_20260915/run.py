"""Freeze and execute one bounded repeat of the supplied H8 discovery recipe.

The optimizer is unchanged. The old compact seed, prepared physical maps,
MPS, and nonsinglet proof are inherited dependencies, not fresh discovery.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'research/h8_discovery_repeat_20260915'
OUT = ROOT/'results/h8_discovery_repeat_20260915'
DISCOVERY = OUT/'discovery'
PACKAGE = ROOT/'results/h8_spin_import_20260914/imported/Spectra_H8_0767448mHa'
NUMERIC = ROOT/'.venv-correlated/bin/python'
STDLIB = Path('/opt/homebrew/Caskroom/miniconda/base/bin/python')
MODULE = 'research.h8_discovery_repeat_20260915.run'
PREPARED = ROOT/'results/sector_quotient_20260914/prepared_linear_closure'
INPUTS = {
    'fixture': ROOT/'results/certificate_scaling/active_space_ladder/h8/fixture.json',
    'state': ROOT/'results/correlated_pair_20260913/mps/h8_spatial_warm144/state.json',
    'old_compact_checkpoint': ROOT/'results/sector_quotient_20260914/candidates/optimize_number_linear_closure/checkpoint.npz',
    'nonsinglet': ROOT/'results/collective_completion_20260914/candidates/h8_spin_r1_nonsinglet/round_0/certificate.json',
    'frozen_endpoints': ROOT/'results/sector_quotient_20260914/frozen_inputs.json',
}


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            value.update(chunk)
    return value.hexdigest()


def dump(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def commands():
    script = lambda name: [str(NUMERIC), '-B', str(PACKAGE/'discovery'/name)]
    return [
        ('mps_moments', 90, script('mps_six_moments.py')),
        ('moment_pullback', 90, script('moment_probe.py')),
        ('sparse_build', 60, script('build_spin_sparse.py')),
        ('stage1', 210, script('solve_spin_sparse.py')+['stage1', '--seconds', '130', '--mu', '2']),
        ('stage2', 270, script('solve_spin_sparse.py')+['stage2', '--seconds', '180', '--mu', '.03',
                                                   '--restart', str(DISCOVERY/'stage1/checkpoint.npz')]),
        ('prune', 30, [str(STDLIB), '-B', '-S', str(PACKAGE/'discovery/prune_tiny_factors.py'),
                       str(DISCOVERY/'stage2/export/certificate.json'), str(OUT/'certificate.json')]),
        ('exact_replay', 300, [str(STDLIB), '-B', '-S', '-m', 'research.h8_discovery_repeat_20260915.replay']),
    ]


def configure():
    os.environ['SPECTRA_ROOT'] = str(ROOT)
    os.environ['SPECTRA_OUTPUT'] = str(DISCOVERY)
    os.environ['NUMBA_CACHE_DIR'] = str(OUT/'numba_cache')
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'NUMBA_NUM_THREADS'):
        os.environ[name] = '1'


def freeze():
    OUT.mkdir(exist_ok=False)
    DISCOVERY.mkdir()
    frame = json.loads((PREPARED/'frame.json').read_text())
    if sha(INPUTS['fixture']) != frame['fixture_sha256'] or sha(PREPARED/'seed.npz') != frame['seed_sha256']:
        raise ValueError('Prepared input binding changed')
    fixture = json.loads(INPUTS['fixture'].read_text())
    if (fixture['modes'], fixture['particles']) != (16, 8):
        raise ValueError('This recipe is specialized to the declared H8 frame')
    supplied = json.loads((PACKAGE/'SHA256.json').read_text())
    files = set(INPUTS.values()) | set(PREPARED.glob('*'))
    for base in (ROOT/'research', ROOT/'experiments'):
        files.update(base.rglob('*.py'))
    files.update((PACKAGE/'discovery').glob('*.py'))
    for path in (PACKAGE/'discovery').glob('*.py'):
        if sha(path) != supplied[str(path.relative_to(PACKAGE))]:
            raise ValueError(('Supplied source changed', str(path)))
    versions = subprocess.run([str(NUMERIC), '-B', '-c',
        'import json,platform,numpy,scipy,numba; print(json.dumps(dict(python=platform.python_version(),'
        'platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))'],
        check=True, capture_output=True, text=True, timeout=30)
    procedure = {
        'kind': 'frozen_H8_warm_discovery_repeat_v1', 'frozen_UTC': utc(),
        'question': 'Does a new optimization from the declared old compact inputs produce another complete interval <=1.6 mHa?',
        'target_width_Ha': '1/625', 'target_is_total_width_not_per_electron': True,
        'procedure': 'Unmodified numerical recipe in the supplied README; two stages, one trajectory, no adaptive restart after a miss.',
        'steps': [{'name': n, 'hard_timeout_seconds': t, 'command': c} for n, t, c in commands()],
        'optimizer_stop': 'Each stage stops at its fixed solve-time budget or the supplied numerical width<1.3 mHa and primal_l2<5e-8 check; only exact replay decides acceptance.',
        'selection': 'Supplied stage2 best-score export, followed by the supplied row-pruning rule. No selection from historical winning factors.',
        'all_failed_or_timed_out_steps_retained': True,
        'inherited_inputs': {k: str(v.relative_to(ROOT)) for k, v in INPUTS.items()},
        'inherited_preparation': str(PREPARED.relative_to(ROOT)),
        'fresh_this_run': ['degree-six MPS moment contraction including an empty Numba cache',
                           'physical-functional pullback', 'spin representation and sparse normal construction',
                           'both optimization stages', 'rational export and pruning', 'exact MPS and both-spin replay'],
        'not_reconstructed': ['integrals and rational Hamiltonian', 'physical coefficient maps and seed frame',
                              'old compact checkpoint', 'MPS state discovery', 'nonsinglet proof discovery'],
        'source_and_input_hashes': {str(p.relative_to(ROOT)): {'sha256': sha(p), 'bytes': p.stat().st_size}
                                   for p in sorted(files) if p.is_file()},
        'parent_manifest': 'results/h8_spin_import_20260914/manifest.json',
        'parent_manifest_sha256': sha(ROOT/'results/h8_spin_import_20260914/manifest.json'),
        'runtime': json.loads(versions.stdout), 'thread_count': 1,
        'discovery_is_cold_from_integrals': False,
        'trial_count': 1,
    }
    dump(OUT/'protocol.json', procedure)
    dump(OUT/'protocol_seal.json', {'sha256': sha(OUT/'protocol.json'), 'sealed_UTC': utc()})
    print(json.dumps({'protocol_frozen': True, 'steps': len(commands()),
                      'source_and_input_files': len(procedure['source_and_input_hashes']),
                      'protocol_sha256': sha(OUT/'protocol.json')}), flush=True)


def check_freeze():
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Sealed completed pass')
    seal = json.loads((OUT/'protocol_seal.json').read_text())
    if sha(OUT/'protocol.json') != seal['sha256']:
        raise ValueError('Frozen procedure changed')
    protocol = json.loads((OUT/'protocol.json').read_text())
    for name, entry in protocol['source_and_input_hashes'].items():
        if sha(ROOT/name) != entry['sha256']:
            raise ValueError(('Frozen source or input changed', name))
    return protocol


def step(name):
    from research.correlated_pair_20260913 import budget
    configure()
    protocol = check_freeze()
    selected = next(s for s in protocol['steps'] if s['name'] == name)
    budget.OUT = OUT/'runs'
    return budget.run(name, selected['hard_timeout_seconds'], selected['command'])


def execute():
    configure()
    started = time.monotonic()
    started_utc = utc()
    freeze()
    outcome = 'all_steps_completed'
    for name, _, _ in commands():
        result = subprocess.run([str(STDLIB), '-B', '-S', '-m', MODULE, 'step', name], cwd=ROOT)
        if result.returncode:
            outcome = 'stopped_after_failed_step'
            break
    record = {'started_UTC': started_utc, 'finished_UTC': utc(),
              'elapsed_seconds_including_freeze_and_step_input_checks': time.monotonic()-started,
              'status': outcome, 'last_step': name,
              'cold_from_integrals': False, 'post_run_diagnostics_and_report_excluded': True}
    dump(OUT/'execution.json', record)
    print(json.dumps(record), flush=True)
    return 0 if outcome == 'all_steps_completed' else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('execute', 'step'))
    parser.add_argument('name', nargs='?')
    args = parser.parse_args()
    sys.exit(execute() if args.action == 'execute' else step(args.name))
