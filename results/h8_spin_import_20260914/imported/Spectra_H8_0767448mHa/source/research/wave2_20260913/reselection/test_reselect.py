import json
from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from reselect import run


def test_reselection_rejects_wrong_sector_size(tmp_path):
    root = Path(__file__).resolve().parents[3]
    fixture = root / 'results/certificate_scaling/active_space_ladder/h8/fixture.json'
    source = root / 'results/wave2_20260913/upper_states/h8_pt2_hf/step_05_upper.json'
    data = json.loads(source.read_text()); data['independent_upper']['states'] = data['independent_upper']['states'][:-1]
    bad = tmp_path / 'bad.json'; bad.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='exactly 1024'):
        run(fixture, bad, tmp_path / 'out')
