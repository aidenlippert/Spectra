"""Higher-metric search seeded by this chain's newly accepted D-metric atoms."""
from pathlib import Path
from unittest.mock import patch
import json
from experiments.marginal_joint_spinflip import flip_label
from experiments.marginal_moment_pricing import MomentDictionary, construct

root = Path('results/marginal_graded_hubbard8')
data = json.loads((root/'hamiltonian.json').read_text())
seed_path = root/'sharp_degree4/proof/certificate.json'
proof = json.loads(seed_path.read_text()); assert proof['hamiltonian'] == data['hamiltonian']
seeds = set()
for part in ('weight_proof', 'numerator_proof'):
    for family, key in [('positive', 'positive_indicators'), ('charge', 'charge_indicators')]:
        for atom in proof[part][key]:
            label = (family, atom['required'], atom['occupied'])
            seeds.add(min(label, tuple(flip_label(label, data['modes']//2))))
original = MomentDictionary.labels

def seeded_labels(self, degree=6):
    return sorted(set(original(self, degree)) | seeds)

out = root/'seeded_quadratic'
with patch.object(MomentDictionary, 'labels', seeded_labels), \
     patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No states')), \
     patch('experiments.marginal_spin_constructor.spin_states', side_effect=AssertionError('No state generation')), \
     patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError('No full-population lift')):
    result = construct(data, '-16', out, feature_degree=2, proof_degree=6, time_limit=120)
result['seed_source'] = str(seed_path)
result['seed_directions'] = len(seeds)
result['scope'] = 'Higher-metric search at gamma=-16 with atom directions seeded from this same chain’s accepted gamma=-18.01 D-metric certificate. Metric and positive weights are solved afresh; this run is not bare unseeded discovery. Full rational export remains acceptance.'
(out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result), flush=True)
