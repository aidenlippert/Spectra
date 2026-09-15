set -eu
cd /home/ubuntu/spectra-nvidia
.venv/bin/python - <<'PY'
import hashlib,json,tarfile
from pathlib import Path
root=Path('results/transfer_solver_20260915/adaptive/h10_correlated_guide/sparse_cudss_adaptive')
assert (root/'discovery.json').exists()
files=[p for p in root.rglob('*') if p.is_file() and 'exact_compiled' not in p.parts]
manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
Path('proposal_download_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
with tarfile.open('proposal_download.tar.gz','w:gz') as stream:
    for p in files+[Path('proposal_download_manifest.json')]:stream.add(p,arcname=str(p),recursive=False)
Path('proposal_download.sha256').write_text(hashlib.sha256(Path('proposal_download.tar.gz').read_bytes()).hexdigest()+'\n')
print(json.dumps({'files':len(files),'bytes':Path('proposal_download.tar.gz').stat().st_size}))
PY
