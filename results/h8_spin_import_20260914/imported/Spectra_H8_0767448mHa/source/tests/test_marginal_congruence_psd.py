import pytest

from experiments.marginal_congruence_psd import CongruenceProofs, matrix_digest


def test_congruence_proves_spd_without_original_diagonal_dominance():
    a = [[2, 3], [3, 5]]
    proof = CongruenceProofs({matrix_digest(a): [[1, -3], [0, 2]]})
    result = proof.check(a)
    assert result['rank'] == 2
    assert result['minimum_exact_margin'] == '2'
    assert proof.finish() == {'supplied': 1, 'used': 1}


@pytest.mark.parametrize('a', [[[1, 2], [2, 1]], [[1, 0], [0, 0]], [[-1]]])
def test_witness_cannot_accept_indefinite_or_singular_matrix(a):
    identity = [[int(i == j) for j in range(len(a))] for i in range(len(a))]
    with pytest.raises(ValueError, match='dominance'):
        CongruenceProofs({matrix_digest(a): identity}).check(a)


@pytest.mark.parametrize('a', [[[0, 1], [1, 0]], [[1, 2], [2, 1]]])
def test_fallback_preserves_indefinite_refusal(a):
    with pytest.raises(ValueError):
        CongruenceProofs({}).check(a)


def test_singular_fallback_and_exact_divisibility_gate():
    result = CongruenceProofs({}).check([[1, 1], [1, 1]])
    assert result['rank'] == 1 and result['nullity'] == 1
    with pytest.raises(ValueError, match='Nonexact'):
        CongruenceProofs({}).check([[1, 0], [0, 1]], initial_divisor=2)


@pytest.mark.parametrize('factor', [[[1, 0], [1, 1]], [[1, 0], [0, 0]],
                                    [[1.0]], [[True]], [[10**31]], [[1, 2]], []])
def test_malformed_factors_refused(factor):
    with pytest.raises(ValueError):
        CongruenceProofs({matrix_digest([[1]]): factor})


def test_stale_witnesses_and_wrong_dimension_refused():
    proof = CongruenceProofs({matrix_digest([[1]]): [[1]]})
    proof.check([[2]])
    with pytest.raises(ValueError, match='Unused'):
        proof.finish()
    with pytest.raises(ValueError, match='dimension'):
        CongruenceProofs({matrix_digest([[1]]): [[1, 0], [0, 1]]}).check([[1]])


@pytest.mark.parametrize('a', [[[1, 2], [0, 1]], [[1.0]], [[True]], [[1, 0]], []])
def test_invalid_matrices_refused(a):
    with pytest.raises(ValueError):
        CongruenceProofs({}).check(a)


def test_positive_initial_divisor_gate_also_applies_to_witness():
    with pytest.raises(ValueError, match='divisor'):
        CongruenceProofs({matrix_digest([[1]]): [[1]]}).check([[1]], 0)


def test_extended_energy_replay_mixes_witness_and_singular_fallback():
    from experiments.marginal_projector_extendibility import replay
    certificate = {
        'kind': 'hubbard_projector_extension_v2', 'chain_sites': 8,
        'target': {'U': '0', 't': '0', 'V': '0'},
        'local_window': {'kind': 'local_hubbard_block_v1', 'sites': 4, 'U': '0', 't': '0', 'V': '0'},
        'vector': {85: 1}, 'windows': 2, 'projector_sum_ceiling': '2',
        'penalty': '0', 'penalized_lower': '-1',
    }
    witnesses = {matrix_digest([[1]]): [[1]]}
    exact = replay(certificate)
    accelerated = replay(certificate, psd_witnesses=witnesses)
    assert accelerated['open_lower_density'] == exact['open_lower_density']
    assert accelerated['local_sum_dimensions'] == 256
    assert accelerated['congruence_witnesses']['used'] == 1
    assert any(b['psd']['nullity'] for b in accelerated['overlap']['sectors'])
    assert any('method' in b['psd'] for b in accelerated['local_sectors'])
    with pytest.raises(ValueError, match='Unused'):
        replay(certificate, psd_witnesses={**witnesses, matrix_digest([[7]]): [[1]]})
    with pytest.raises(ValueError):
        replay({**certificate, 'penalized_lower': '1'}, psd_witnesses=witnesses)
