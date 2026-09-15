"""A preceding family cap cannot be attached to a three-spectator energy."""
from pathlib import Path
import json
import subprocess
import sys


def test_two_spectator_family_driver_refuses_three_spectator_energy(tmp_path):
    root=Path(__file__).resolve().parents[1]
    p=root/'results/marginal_graded_hubbard8/two_spectator/W_zero/final'
    (tmp_path/'diagonal_family_limit_proposal.json').write_bytes((p/'diagonal_family_limit_proposal.json').read_bytes())
    c=json.loads((p/'profile_joint_r1_2_certificate.json').read_text())
    c['kind']='hubbard_projector_extension_v16';c['three_spectator_hopping']={'0,1,2,3,4,1,1,2':'1/1000'}
    (tmp_path/'profile_joint_r1_2_certificate.json').write_text(json.dumps(c))
    result=subprocess.run([sys.executable,str(root/'results/marginal_graded_hubbard8/discovery/range_two_family_limit_replay.py'),str(tmp_path)],capture_output=True,text=True,timeout=90)
    assert result.returncode!=0
    assert 'Three-spectator energy requires a matching enlarged family cap' in result.stderr
    assert not (tmp_path/'range_two_family_limit_replay.json').exists()
