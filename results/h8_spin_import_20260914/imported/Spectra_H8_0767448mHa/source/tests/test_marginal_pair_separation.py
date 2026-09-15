"""The separation claim must keep the original fixed certificate family."""
import copy
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
path = ROOT/'results/marginal_graded_hubbard8/discovery/pair_transfer_separation.py'
spec = importlib.util.spec_from_file_location('pair_separation',path)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


def inputs():
    p = ROOT/'results/marginal_graded_hubbard8/pair_transfer/W_zero'
    return [json.loads((p/name).read_text()) for name in
            ['profile_joint_r1_2_certificate.json','previous_family/profile_joint_r1_2_certificate.json',
             'previous_family/range_two_family_limit_certificate.json']]


def test_actual_fixed_family_matches():
    module.match(*inputs())


@pytest.mark.parametrize('field', ['target','chain_sites','vector','windows','projector_sum_ceiling','joint','telescoping_diagonal'])
def test_changed_comparison_scope_refused(field):
    new,old,family = inputs(); new = copy.deepcopy(new)
    if field == 'target': new[field]['W'] = '1'
    elif field == 'chain_sites': new[field] += 2
    elif field == 'vector': new[field] = {'63':1}
    elif field == 'windows': new[field] += 1
    elif field == 'projector_sum_ceiling': new[field] = '2'
    elif field == 'joint': new[field]['ratio'] = '1'
    else: new[field]['0'] = '1/1000'
    with pytest.raises(ValueError): module.match(new,old,family)


def test_old_family_matching_driver_refuses_pair_energy(tmp_path):
    p = ROOT/'results/marginal_graded_hubbard8/pair_transfer/W_zero'
    (tmp_path/'diagonal_family_limit_proposal.json').write_bytes(
        (p/'previous_family/diagonal_family_limit_proposal.json').read_bytes())
    (tmp_path/'profile_joint_r1_2_certificate.json').write_bytes(
        (p/'profile_joint_r1_2_certificate.json').read_bytes())
    result = subprocess.run([sys.executable,str(ROOT/'results/marginal_graded_hubbard8/discovery/range_two_family_limit_replay.py'),str(tmp_path)],
                            capture_output=True,text=True,timeout=60)
    assert result.returncode != 0
    assert 'Pair-transfer energy requires a matching enlarged family cap' in result.stderr
    assert not (tmp_path/'range_two_family_limit_replay.json').exists()
