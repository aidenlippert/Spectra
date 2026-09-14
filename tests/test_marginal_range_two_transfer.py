from fractions import Fraction as F
import pytest

from experiments.marginal_range_two_transfer import compile_state, contract, enclose
from experiments.marginal_boundary_transfer import compile_state as compile_nearest, contract as nearest_contract
from experiments.marginal_transfer_verify import apply_word


SOURCE = {15: 1, 51: 2, 204: -1, 240: 1}


def hopping_image(vector, bonds):
    result = {}
    for s, a in vector.items():
        for left, right, coefficient in bonds:
            for spin in (0, 1):
                for i, j in ((left, right), (right, left)):
                    image = apply_word(((1, 2*i+spin), (0, 2*j+spin)), s)
                    if image:
                        target, sign = image
                        result[target] = result.get(target, F(0))+a*coefficient*sign
    return {s: a for s, a in result.items() if a}


def physical_state(blocks, a, b):
    vector = {0: F(1)}
    for block in range(blocks):
        vector = {s+(r << (8*block)): x*y for s, x in vector.items() for r, y in SOURCE.items()}
    for cut in range(1, blocks):
        first = hopping_image(vector, [(4*cut-1, 4*cut, -1)])
        second = hopping_image(first, [(4*cut-1, 4*cut, -1)])
        vector = {s: vector.get(s, 0)-a*first.get(s, 0)+b*second.get(s, 0)
                  for s in set(vector) | set(first) | set(second)}
        vector = {s: x for s, x in vector.items() if x}
    return vector


@pytest.mark.parametrize('blocks,W', [(1, F(1, 7)), (2, F(1, 7)), (3, F(-1, 5))])
def test_range_two_transfer_matches_expanded_filtered_physical_state(blocks, W):
    U, t, V, a, b = F(3), F(2, 3), F(-1, 4), F(1, 10), F(1, 50)
    compiled = compile_state(SOURCE, 4, U, t, V, W)
    result = contract(compiled, a, b, blocks)
    vector = physical_state(blocks, a, b)
    sites = 4*blocks
    norm = sum(x*x for x in vector.values())
    numerator = sum(x*vector.get(s, 0) for s, x in hopping_image(vector, [(i, i+1, -t) for i in range(sites-1)]).items())
    for s, x in vector.items():
        occupations = [(s >> (2*i)) & 3 for i in range(sites)]
        charges = [v.bit_count()-1 for v in occupations]
        diagonal = U*sum(v == 3 for v in occupations)
        diagonal += V*sum(charges[i]*charges[i+1] for i in range(sites-1))
        diagonal += W*sum(charges[i]*charges[i+2] for i in range(sites-2))
        numerator += x*x*diagonal
    assert F(result['norm']) == norm
    assert F(result['energy']) == numerator/norm
    interval = enclose(compiled, a, b, blocks, 128)
    assert F(interval['energy_lower']) <= numerator/norm <= F(interval['energy_upper'])


def test_zero_range_two_agrees_with_original_transfer():
    original = nearest_contract(compile_nearest(SOURCE, 4, 4, 1, F(1, 2)), F(1, 10), F(1, 50), 3)
    extended = contract(compile_state(SOURCE, 4, 4, 1, F(1, 2), 0), F(1, 10), F(1, 50), 3)
    assert extended['norm'] == original['norm']
    assert extended['energy'] == original['energy']
    assert extended['target']['W'] == '0'


@pytest.mark.parametrize('W', [0.1, True, '1/1000001', '1000001'])
def test_invalid_range_two_coefficients_refused(W):
    with pytest.raises(ValueError):
        compile_state(SOURCE, 4, W=W)
