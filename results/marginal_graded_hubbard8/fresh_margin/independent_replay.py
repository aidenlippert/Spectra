"""Independent exact replay of the fresh margin certificate."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[3]))
from experiments.marginal_polynomial_metric import replay

HERE = Path(__file__).parent
certificate = json.loads((HERE / 'proof' / 'certificate.json').read_text())
receipt = replay(certificate)
result = {'replayed': True, 'receipt': receipt,
          'scope': 'Independent exact replay of this certificate; no numerical solver result is treated as proof.'}
(HERE / 'independent_replay.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
