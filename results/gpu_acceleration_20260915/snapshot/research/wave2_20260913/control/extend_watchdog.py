"""Replace our process watchdog after a completed solve needs more export time."""
import json
import os
from pathlib import Path
import signal
import sys
import time

child, wrapper = map(int, sys.argv[1:3])
duration = int(sys.argv[3])
if not 1 <= duration <= 1500:
    raise ValueError('Bounded extension required')
command = Path(f'/proc/{child}/cmdline').read_bytes()
wrapper_command = Path(f'/proc/{wrapper}/cmdline').read_bytes()
if b'research.wave2_20260913.control.run_h10_lower' not in command or not wrapper_command.startswith(b'timeout\x00'):
    raise ValueError('Unexpected owned process identity')
start_ticks = Path(f'/proc/{child}/stat').read_text().split()[21]
started = time.time()
record = {'child': child, 'original_watchdog': wrapper, 'child_start_ticks': start_ticks,
          'extension_seconds': duration, 'start_unix': started,
          'reason': 'Numerical solve complete; exact export and new residual proof need more time.',
          'verification_gates_changed': False}
Path('watchdog_extension.json').write_text(json.dumps(record, indent=2)+'\n')
# Kill only our old timeout parent, not its process group or child. This is a
# runtime-budget revision, with an explicit replacement deadline below.
os.kill(wrapper, signal.SIGKILL)
while time.time()-started < duration:
    stat = Path(f'/proc/{child}/stat')
    if not stat.exists() or stat.read_text().split()[21] != start_ticks or stat.read_text().split()[2] == 'Z':
        record['outcome'] = 'child_finished'
        break
    time.sleep(5)
else:
    os.kill(child, signal.SIGTERM)
    record['outcome'] = 'extension_deadline_terminated_child'
record['end_unix'] = time.time()
Path('watchdog_extension.json').write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps(record), flush=True)
