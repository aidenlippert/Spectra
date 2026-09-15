"""Run outstanding accepting checks as separately timed, bounded processes."""
import json
from pathlib import Path
import subprocess
import sys
from research.reconstruction_compression_20260914.inputs import OUT

def main():
    frozen=json.loads((OUT/'frozen_inputs.json').read_text())
    for cert in sorted((OUT/'candidates').glob('*/certificate.json')):
        output=cert.parent/'interval.json'
        if output.exists():continue
        disc=json.loads((cert.parent/'discovery.json').read_text());case=disc['case']
        if case not in ('h6_actual','h8'):continue
        inputs=frozen['h6' if case=='h6_actual' else 'h8']
        receipts=[json.loads(p.read_text()) for p in (OUT/'runs').glob('*.json')]
        used=sum(r.get('wall_seconds',0) for r in receipts)
        if used+90>2400:raise RuntimeError('Aggregate resource cap would be exceeded')
        name='accept_'+cert.parent.name
        command=[sys.executable,'-m','research.reconstruction_compression_20260914.budget','--name',name,'--seconds','90','--',
                 sys.executable,'-m','research.reconstruction_compression_20260914.replay',inputs['fixture'],str(cert),str(output),
                 '--upper-receipt',inputs['upper_receipt']]
        print(json.dumps({'next_exact_replay':cert.parent.name,'completed_process_seconds':used}),flush=True)
        subprocess.run(command,check=True)

if __name__=='__main__':main()
