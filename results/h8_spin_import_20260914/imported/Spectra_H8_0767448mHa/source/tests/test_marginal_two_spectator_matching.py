"""Old family caps cannot be attached to the new energy operator class."""
from pathlib import Path
import subprocess
import sys


def test_pair_family_driver_refuses_two_spectator_energy(tmp_path):
    root=Path(__file__).resolve().parents[1]
    p=root/'results/marginal_graded_hubbard8/two_spectator/W_zero'
    (tmp_path/'diagonal_family_limit_proposal.json').write_bytes(
        (p/'previous_family/diagonal_family_limit_proposal.json').read_bytes())
    (tmp_path/'profile_joint_r1_2_certificate.json').write_bytes(
        (p/'profile_joint_r1_2_certificate.json').read_bytes())
    result=subprocess.run([sys.executable,str(root/'results/marginal_graded_hubbard8/discovery/range_two_family_limit_replay.py'),str(tmp_path)],
                          capture_output=True,text=True,timeout=60)
    assert result.returncode!=0
    assert 'Two-spectator energy requires a matching enlarged family cap' in result.stderr
    assert not (tmp_path/'range_two_family_limit_replay.json').exists()
