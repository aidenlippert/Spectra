"""Independent command-line acceptance receipt for the generated certificate."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.marginal_joint_charge_dp import replay

base = Path(__file__).resolve().parent
receipt = replay(json.loads((base / 'certificate.json').read_text()))
(base / 'independent_replay.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))
