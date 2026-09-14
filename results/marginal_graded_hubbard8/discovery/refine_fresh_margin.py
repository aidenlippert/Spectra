"""Fresh bare-H margin refinement on the crossover atom-direction checkpoint.

This driver imports only the 4,148 canonical atom directions.  Metric
coefficients and proof weights are rebuilt, then exact export is the sole
acceptance gate.
"""
from pathlib import Path
import json
import sys
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).parents[3]))
from joint_reynolds_search import construct

ROOT = Path(__file__).parents[1]
CHECKPOINT = ROOT / "joint_reynolds_crossover" / "checkpoint.json"
OUT = ROOT / "fresh_margin"

def run():
    data = json.loads((ROOT / "hamiltonian.json").read_text())
    with patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No states')), \
         patch('experiments.marginal_spin_constructor.spin_states', side_effect=AssertionError('No state generation')), \
         patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError('No full-population lift')), \
         patch('experiments.marginal_joint_coefficient_constructor.prepare', side_effect=AssertionError('No full atom dictionary')):
        return construct(data, '-14', OUT, time_limit=180, batch=128,
                         hard_normalization=True, method='ipm', audit=True,
                         coordinates='coefficient', weight_degree=6,
                         resume=CHECKPOINT, refine=True, crossover=True,
                         weight_floor=.001, numerator_floor=.001)

if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
