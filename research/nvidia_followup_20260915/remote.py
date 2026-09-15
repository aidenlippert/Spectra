"""Run bounded commands on this experiment's owned instance and retain receipts."""
import argparse
import datetime
import json
from pathlib import Path
import shlex
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/nvidia_followup_20260915'
KEY = Path('/Users/aidenlippert/.ssh/id_ed25519')


def run(name, script, seconds):
    if not name.replace('_', '').isalnum() or (OUT/f'{name}.json').exists():
        raise ValueError('A new alphanumeric receipt name is required')
    if (OUT/'termination.json').exists():
        raise ValueError('Instance termination already requested')
    owned = json.loads((OUT/'instance.json').read_text())
    ready = json.loads((OUT/'ready.json').read_text())
    if ready['id'] != owned['id']:
        raise ValueError('Owned instance mismatch')
    command = ['ssh', '-i', str(KEY), '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15',
        '-o', 'StrictHostKeyChecking=yes', '-o', f'UserKnownHostsFile={OUT}/known_hosts',
        f"ubuntu@{ready['ip']}", 'timeout', '--signal=TERM', '--kill-after=15s',
        str(seconds), 'bash', '-se']
    started = time.monotonic()
    record = {'name': name, 'instance_id': owned['id'], 'started_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'timeout_seconds': seconds, 'script': script, 'status': 'running'}
    (OUT/f'{name}.json').write_text(json.dumps(record, indent=2)+'\n')
    with (OUT/f'{name}.log').open('x') as log:
        try:
            proc = subprocess.run(command, input=script, text=True, stdout=log, stderr=subprocess.STDOUT,
                timeout=seconds+45)
            record.update(exit_code=proc.returncode, status='passed' if proc.returncode == 0 else 'failed')
        except subprocess.TimeoutExpired:
            record.update(exit_code=None, status='transport_timeout_remote_limit_still_applies')
    record.update(wall_seconds=time.monotonic()-started,
        finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (OUT/f'{name}.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record), flush=True)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('name')
    parser.add_argument('script', type=Path)
    parser.add_argument('--seconds', type=int, default=900)
    args = parser.parse_args()
    run(args.name, args.script.read_text(), args.seconds)
