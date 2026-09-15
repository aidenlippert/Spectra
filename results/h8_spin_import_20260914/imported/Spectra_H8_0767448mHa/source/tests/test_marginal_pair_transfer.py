"""Independent pair bit swaps, full-space projection and energy refusal gates."""
from fractions import Fraction as F
import pytest
from experiments.marginal_pair_transfer import LABELS, coefficients, actions
from experiments.marginal_hopping_telescope import projected_matrix
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_charged_projectors import _particlehole, _spinflip
from test_marginal_hopping_telescope import toy


def direct(state, terms):
    image = {}
    for key,value in terms.items():
        i,j = LABELS[key]
        for a,b,sign in [(i,j,1),(4-j,4-i,-1),(i+1,j+1,-1),(5-j,5-i,1)]:
            occupations = ((state>>(2*a))&3, (state>>(2*b))&3)
            if occupations in ((3,0),(0,3)):
                target = state ^ (3<<(2*a)) ^ (3<<(2*b))
                image[target] = image.get(target,0) + value*sign
    return {s:a for s,a in image.items() if a}


def test_full_fock_pair_signs_and_discrete_symmetries():
    for terms in [{key:1} for key in LABELS] + [{key:F(i+1,17) for i,key in enumerate(LABELS)}]:
        data = actions(terms)
        for s in range(4096):
            assert data[s] == direct(s,terms)
            for transform in (_particlehole,_spinflip):
                us,phase = transform(s)
                right = {transform(t)[0]: transform(t)[1]*a for t,a in data[s].items()}
                assert right == {t:phase*a for t,a in data[us].items()}


def test_projection_reconstructs_entire_physical_image():
    terms = {key:F(i+1,17) for i,key in enumerate(LABELS)}
    data = actions(terms)
    for _,_,cols in sector_matrices(6,0,0).values():
        matrix = projected_matrix(cols,data)
        for j,col in enumerate(cols):
            image = {}
            for s,a in col.items():
                for t,b in direct(s,terms).items(): image[t] = image.get(t,0)+a*b
            rebuilt = {s:F(matrix[i][j]*a,sum(v*v for v in basis.values()))
                       for i,basis in enumerate(cols) for s,a in basis.items()}
            assert {s:a for s,a in image.items() if a} == {s:a for s,a in rebuilt.items() if a}


def test_translates_cancel_for_arbitrary_periodic_length():
    for sites in (8,9,10,17,101):
        for i,j in LABELS.values():
            totals = {}
            for offset in range(sites):
                for a,b,sign in [(i,j,1),(4-j,4-i,-1),(i+1,j+1,-1),(5-j,5-i,1)]:
                    key = tuple(sorted(((a+offset)%sites,(b+offset)%sites)))
                    totals[key] = totals.get(key,0)+sign
            assert not any(totals.values())


def pair_toy():
    c = toy()
    c.update(kind='hubbard_projector_extension_v14', spin_telescope={'0,1':'1/1000'},
             spectator_hopping={'0,1,2':'1/1000'}, pair_transfer={'0,1':'1/1000'})
    return c


def test_energy_full_fock_acceptance_and_negative_psd_refusal():
    c = pair_toy(); r = replay(c)
    assert r['accepted'] and r['pair_transfer'] == c['pair_transfer']
    assert r['local_sum_dimensions'] == 4096 and len(r['local_sectors']) == 94
    assert r['local_maximum_psd_dimension'] == 200
    assert r['target'] == c['target']
    assert r['open_lower_density'] == r['periodic_lower_density'] == '-9/5'
    # A large actual pair correction must enter the PSD matrix and be refused.
    c['pair_transfer'] = {'0,1':'100'}
    with pytest.raises(ValueError): replay(c)


@pytest.mark.parametrize('version', [1,2,4,5,6,7,8,9,10,11,12,13])
def test_old_version_cannot_silently_ignore_pair_transfer(version):
    c = pair_toy(); c['kind'] = f'hubbard_projector_extension_v{version}'
    with pytest.raises(ValueError,match='requires v14'): replay(c)


@pytest.mark.parametrize('source', [None,{}, {'0,0':1},{'3,4':1},{'0,1':0},
                                  {'0,1':True},{'0,1':.1},{'0,1':'1/1000001'},
                                  {'0,1':'1000001'}])
def test_invalid_exact_pair_coefficients_refused(source):
    with pytest.raises(ValueError): coefficients(source)


def test_existing_family_cannot_claim_pair_constraints():
    from test_marginal_spectator_hopping import family
    from experiments.marginal_range_two_family_limit import replay as family_replay
    c = family(); c['pair_transfer'] = {'0,1':'1/1000'}
    with pytest.raises(ValueError,match='constrain moments'): family_replay(c)
