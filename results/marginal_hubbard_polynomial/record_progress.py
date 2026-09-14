"""Record final test/replay evidence and refresh only this work's manifest paths."""
from pathlib import Path
import hashlib
import json
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'results/marginal_hubbard_polynomial'
LOG = Path('/tmp/marginal_analytic_hubbard_validation.log')


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


log = LOG.read_text()
match = re.search(r'Ran (\d+) tests in ([\d.]+)s\s+OK\s*$', log)
if not match or int(match.group(1)) != 353:
    raise ValueError('Expected completed353-test full-suite evidence')
tests, seconds = int(match.group(1)), float(match.group(2))
shutil.copyfile(LOG, BASE/'full_validation.log')
checks = []
for row in read(BASE/'scaling.json')['rows']:
    folder = BASE/str(row['sites'])
    receipt = read(folder/'receipt.json')
    assert receipt == read(folder/'independent_replay.json')
    assert receipt['weight_positivity']['residual_l1'] == '0'
    assert receipt['numerator_positivity']['residual_l1'] == '0'
    checks.append(str((folder/'independent_replay.json').relative_to(ROOT)))
assert len(checks) == 7
progress = {
    'status': 'Exact directly generated Hubbard Q-gap certificates accepted through32 sites; general goal remains open.',
    'report': 'research/marginal_hubbard_polynomial.md',
    'tests_passed': tests, 'elapsed_seconds': seconds,
    'stdlib_replays': checks,
    'construction': 'Explicit degree<=4 positive indicators and degree<=2 number multipliers; no optimizer, configuration enumeration or full-population ideal lift.',
    'scaling': 'O(m^2) positive atoms and O(m) number-multiplier terms for this chain formula; all measured exact residuals zero.',
    'accuracy_scope': 'Loose bound U-4t(m-1)-epsilon, not an accurate ground-energy interval or an improvement over elementary uniform row estimates.',
    'remaining': ['Combine controlled construction size with accurate bounds for molecular Hamiltonians.',
                  'Remove explicit finite-sector metric/proof discovery from the optimized H6 path.',
                  'Compress retained response and upper-witness work with checked error.',
                  'General marginal representability, finite-temperature, dynamics and synthesis remain open.']}
save(BASE/'progress.json', progress)
vpath = ROOT/'results/marginal_final_validation.json'
validation = read(vpath)
validation.update(tests_passed=tests, elapsed_seconds=seconds,
    command="PYTHONPYCACHEPREFIX=/tmp/spectra_quantum_pycache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_marginal*.py' -q",
    latest_validation_scope='353 full marginal tests;17 new standard-library replays: four joint/vanishing H6 gaps, four existing-witness energies, two fixed-old-metric SOS obstructions, seven analytic Hubbard-chain Q gaps. Four new energy-ledger entries; no general representability or accuracy/scaling theorem.',
    analytic_hubbard_report='research/marginal_hubbard_polynomial.md',
    analytic_hubbard_progress='results/marginal_hubbard_polynomial/progress.json',
    analytic_hubbard_stdlib_replays=checks)
assert len(validation['stdlib_replays_passed']) == 142
save(vpath, validation)
jpath = ROOT/'results/marginal_h6/polynomial_metric/joint_progress.json'
joint = read(jpath)
joint.update(latest_full_tests_passed=tests, latest_full_elapsed_seconds=seconds,
             analytic_hubbard_followup='results/marginal_hubbard_polynomial/progress.json')
save(jpath, joint)
for name in ['research/marginal_blocker_progress.md', 'research/marginal_polynomial_metric.md']:
    path = ROOT/name
    text = path.read_text().replace('**350 full marginal tests pass.**', '**353 full marginal tests pass.**')
    path.write_text(text)
report = ROOT/'research/marginal_hubbard_polynomial.md'
text = report.read_text()
if 'Full marginal validation:' not in text:
    text += '\nFull marginal validation: **353 tests pass**; seven new chain certificates independently replay with the standard library.\n'
report.write_text(text)
paths = set([
    'experiments/marginal_polynomial_metric.py', 'experiments/marginal_polynomial_sos.py',
    'experiments/marginal_hubbard_polynomial.py', 'tests/test_marginal_polynomial_metric.py',
    'tests/test_marginal_polynomial_sos.py', 'tests/test_marginal_hubbard_polynomial.py',
    'research/marginal_polynomial_metric.md', 'research/marginal_joint_metric.md',
    'research/marginal_projector_degree.md', 'research/marginal_hubbard_polynomial.md',
    'research/marginal_blocker_progress.md', 'results/marginal_final_validation.json'])
for base in [ROOT/'results/marginal_h6/polynomial_metric', BASE]:
    paths.update(str(p.relative_to(ROOT)) for p in base.rglob('*') if p.is_file())
mpath = ROOT/'results/marginal_progress_manifest.json'
manifest = read(mpath)
for path in sorted(paths):
    manifest['sha256'][path] = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
manifest.update(file_count=len(manifest['sha256']), last_refreshed_paths=sorted(paths),
                last_refresh_scope='Joint and vanishing H6 metric proofs, fixed-old-metric quadratic-SOS separators, projector-degree argument, analytic Hubbard-chain certificates, validated source/tests/reports and labeled numerical discovery.')
save(mpath, manifest)
for path in paths:
    assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == manifest['sha256'][path]
print(json.dumps({'tests': tests, 'seconds': seconds, 'new_stdlib_replays':17,
                  'energy_witness_ledger':142, 'manifest_files':manifest['file_count'],
                  'refreshed_paths':len(paths), 'general_goal_complete':False}, indent=2))
