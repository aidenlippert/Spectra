set -eu
cd /home/ubuntu/spectra-nvidia
.venv/bin/python - <<'PY'
import hashlib,json,tarfile
from pathlib import Path
for name,key in [('source','sha256'),('cases','archive_sha256')]:
    data=json.loads(Path(name+'_upload_manifest.json').read_text())
    archive=Path(name+'.tar.gz')
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==data[key]
    with tarfile.open(archive) as stream:
        stream.extractall('.',filter='data')
    if isinstance(data.get('files'),dict):
        for path, expected in data['files'].items():
            assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected,path
print('Uploaded archive and case hashes verified')
import cupy as cp, numpy as np, nvmath, inspect, flint
from cuquantum.tensornet import Network
from nvmath.sparse.advanced import DirectSolver, DirectSolverOptions, DirectSolverMatrixType, DirectSolverMatrixViewType
from sksparse.cholmod import cholesky
print('CuPy',cp.__version__,'nvmath',nvmath.__version__,'flint',flint.__version__)
print('GPU',cp.cuda.runtime.getDeviceProperties(0)['name'])
print(inspect.signature(DirectSolverOptions))
print(list(DirectSolverMatrixType),list(DirectSolverMatrixViewType))
print('nvmath reset',inspect.signature(DirectSolver.reset_operands))
print('tensor network',inspect.signature(Network))
print('FP64 eigh',cp.asnumpy(cp.linalg.eigh(cp.eye(4,dtype=cp.float64))[0]).tolist())
PY
