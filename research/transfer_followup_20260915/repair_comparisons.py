"""Repair only NumPy-scalar JSON export, retaining the failed files and jobs."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, CHEM, PREFIX


def typed_dump(target, value):
    native = json.loads(json.dumps(value, default=lambda scalar: scalar.item()))
    # The frozen numerical function prints the same record after writing it.
    value.clear()
    value.update(native)
    dump(target, value)


def run_case(case):
    path = case/'fci_reference.json'
    raw = path.read_bytes()
    try:
        json.loads(raw)
    except json.JSONDecodeError:
        pass
    else:
        raise ValueError('Refusing to replace a complete reference record')
    folder = case/'invalidated_fci_serialization'
    folder.mkdir(exist_ok=False)
    path.rename(folder/'fci_reference.json.partial')
    dump(folder/'repair.json', {'failed_file_sha256': hashlib.sha256(raw).hexdigest(),
        'reason': 'A NumPy bool in the comparison field is not accepted by the standard JSON encoder',
        'change': 'Convert NumPy scalar values to their native Python scalar for serialization only',
        'frozen_benchmark_source_unchanged': True, 'numerical_solver_and_tolerances_unchanged': True,
        'reference_recomputed_after_discovery': True})
    from research.transfer_solver_20260915 import benchmark

    benchmark.dump = typed_dump
    benchmark.run(case)


def run():
    outcomes = []
    for name in ('h4_control', 'h8_cold', 'h6_asymmetric', 'water_asymmetric'):
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', name+'_fci_fixed',
            '--seconds', '180', '--', CHEM, '-B', '-m', 'research.transfer_followup_20260915.repair_comparisons',
            'case', str(OUT/'cases'/name)], cwd=ROOT)
        outcomes.append({'case': name, 'exit_code': result.returncode})
    dump(OUT/'fci_export_repairs.json', outcomes)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'case'))
    p.add_argument('case', type=Path, nargs='?')
    a = p.parse_args()
    run() if a.action == 'run' else run_case(a.case.resolve())
