set -eu
cd /home/ubuntu/spectra-nvidia
.venv/bin/python - <<'PY'
import hashlib,json,tarfile
from pathlib import Path
files=[Path('installed.txt'),Path('gpu.txt')]
for root in [Path('results/transfer_solver_20260915/adaptive/h10_correlated_guide'),Path('results/transfer_solver_20260915/cases/h8_cold'),Path('results/nvidia_followup_20260915')]:
    files.extend(p for p in root.rglob('*') if p.is_file() and not set(p.parts)&{'prepared','__pycache__','numba_cache','mps'})
files=sorted(set(files))
manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
Path('final_download_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
with tarfile.open('final_download.tar.gz','w:gz') as stream:
    for p in files+[Path('final_download_manifest.json')]:stream.add(p,arcname=str(p),recursive=False)
Path('final_download.sha256').write_text(hashlib.sha256(Path('final_download.tar.gz').read_bytes()).hexdigest()+'\n')
print(json.dumps({'files':len(files),'bytes':Path('final_download.tar.gz').stat().st_size}))
PY
