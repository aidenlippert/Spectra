import pytest

from research.rsi_discovery_20260916.v2.correlation_search import validate_design, select
from research.rsi_discovery_20260916.v2.workflow import failure_state


def test_nonlocal_support_selection_does_not_require_a_global_map():
    data = {"modes": 12, "particles": 6, "hamiltonian": [
        {"word": [[1, 0], [0, 0]], "coefficient": "1/3"}]}
    source = '''def propose(payload):
    n = payload["spatial"]
    return {"clusters": [[0, 1, 2]], "collective_clusters": [[0, n-2, n-1]]}
'''
    assert select(data, source) == {"clusters": [[0, 1, 2]], "collective_clusters": [[0, 4, 5]]}


def test_representation_cannot_edit_target_or_exceed_the_declared_support_budget():
    with pytest.raises(ValueError, match="only local and collective"):
        validate_design({"clusters": [[0, 1]], "collective_clusters": [], "particles": 2}, 4)
    with pytest.raises(ValueError, match="original orbital domain"):
        validate_design({"clusters": [[0, 4]], "collective_clusters": []}, 4)
    with pytest.raises(ValueError, match="budget"):
        validate_design({"clusters": [[0, 1], [0, 2], [1, 2]], "collective_clusters": []}, 4)


def test_physical_failures_keep_resource_numerical_and_semantic_meanings_separate():
    assert failure_state("solve_0", True, "") == "resource_obstruction"
    assert failure_state("prepare_0", False, "Gram allocation envelope exceeded") == "resource_obstruction"
    assert failure_state("solve_0", False, "No numerical iterate produced") == "numerical_failure"
    assert failure_state("solve_0", False, "ModuleNotFoundError") == "implementation_failure"
    assert failure_state("prepare_0", False, "Proposed map changed an actual solver coefficient") == "counterexample_to_encoded_claim"
