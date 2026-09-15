"""Freeze changed-input tests together and run every predeclared case."""
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import time
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import CASES, STD


def run():
    names = ['h8_cold', 'h6_asymmetric', 'water_asymmetric', 'h10_size']
    sources = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in sorted((ROOT/'research/transfer_solver_20260915').glob('*.py'))}
    record = {'frozen_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'cases': {name: CASES[name] for name in names}, 'source_hashes': sources,
              'development_case': 'h4_control', 'old_achieved_results_preserved': True,
              'target_total_width_Ha': '1/625', 'all_attempts_reported': True,
              'changed_input_cases_must_use_same_source': True,
              'cold_from_new_integrals_and_empty_state_cache': True,
              'no_existing_coefficient_maps_states_or_optimization_outputs': True,
              'not_a_statistical_reliability_sample': True}
    dump(OUT/'validation_set_protocol.json', record)
    started = time.monotonic()
    outcomes = []
    for name in names:
        for source, digest in sources.items():
            if hashlib.sha256((ROOT/source).read_bytes()).hexdigest() != digest:
                raise ValueError('Frozen validation-set source changed')
        result = subprocess.run([STD, '-B', '-S', '-m', 'research.transfer_solver_20260915.campaign', name], cwd=ROOT)
        outcomes.append({'case': name, 'exit_code': result.returncode})
        (OUT/'validation_set_progress.json').write_text(json.dumps(outcomes, indent=2)+'\n')
        print(json.dumps(outcomes[-1]), flush=True)
    dump(OUT/'validation_set_execution.json', {'seconds': time.monotonic()-started, 'outcomes': outcomes})


if __name__ == '__main__':
    run()
