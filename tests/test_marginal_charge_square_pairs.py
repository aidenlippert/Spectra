from fractions import Fraction as F
import pytest

from experiments.marginal_charge_square_pairs import coefficients, five_site_value, local_value, LABELS
from experiments.marginal_projector_extendibility import replay


def direct(state, terms):
    q = [((state // 4**i) % 4 == 3)-((state // 4**i) % 4 == 0) for i in range(6)]
    value = F(0)
    for label, a in terms.items():
        i,j=map(int,label.split(','))
        value += a*(q[i]**2*q[j]**2-q[4-j]**2*q[4-i]**2-q[i+1]**2*q[j+1]**2+q[5-j]**2*q[5-i]**2)
    return value


def test_all_four_pairs_match_full_fock_expansion_and_psd_minimum():
    terms = {label: F(i+1,7) for i,label in enumerate(LABELS)}
    assert len(terms) == 4
    values=[]
    for s in range(4096):
        assert local_value(s,terms)==direct(s,terms)
        q=[((s >> (2*i)) & 3).bit_count()-1 for i in range(6)]
        quadratic=F(q[0]*q[3]-2*q[1]*q[4]+q[2]*q[5],7)
        sparse=int(s%1024==102)-int(s%1024==612)-int(s//4==102)+int(s//4==612)
        values.append(direct(s,terms)+quadratic+sparse)
    c={'kind':'hubbard_projector_extension_v8','chain_sites':10,
       'target':{'U':'0','t':'0','V':'0','W':'0'},
       'local_window':{'kind':'local_hubbard_range2_block_v1','sites':6,'U':'0','t':'0','V':'0'},
       'vector':{0x999:1,0x666:1},'windows':2,'projector_sum_ceiling':'2','penalty':'0','penalized_lower':str(min(values)),
       'joint':{'vector':{62:1,3008:1},'windows':2,'ratio':'1/2','projector_sum_ceiling':'2','penalty':'0'},
       'telescoping_diagonal':{102:1,612:-1},'quadratic_charge_telescope':{'0,3':'1/7'},
       'charge_square_pair_telescope':{key:str(value) for key,value in terms.items()}}
    result=replay(c)
    assert result['local_sum_dimensions']==4096
    assert result['local_maximum_psd_dimension']==200
    assert result['charge_square_pair_components']==4
    c['penalized_lower']=str(min(values)+F(1,7))
    with pytest.raises(ValueError):replay(c)
    c['kind']='hubbard_projector_extension_v7'
    with pytest.raises(ValueError,match='require v8'):replay(c)


def test_square_pair_symmetry_and_independent_periodic_cancellation():
    terms=coefficients({label:i+1 for i,label in enumerate(LABELS)})
    for s in range(1024):
        reflected=sum(((s>>(2*i))&3)<<(2*(4-i)) for i in range(5))
        assert five_site_value(reflected,terms)==-five_site_value(s,terms)
        assert five_site_value(s^1023,terms)==five_site_value(s,terms)
    for bits in range(256):
        p=[(bits>>i)&1 for i in range(8)]
        for i,j in LABELS.values():
            assert sum(p[(k+i)%8]*p[(k+j)%8]-p[(k+4-j)%8]*p[(k+4-i)%8]-p[(k+i+1)%8]*p[(k+j+1)%8]+p[(k+5-j)%8]*p[(k+5-i)%8] for k in range(8))==0


@pytest.mark.parametrize('terms',[{}, {'0,0':1},{'1,1':1},{'0,4':1},{'3,0':1},{'0,1':0},{'0,1':0.1},{'0,1':True},{'0,1':'1/1000001'}])
def test_noncanonical_square_pairs_refused(terms):
    with pytest.raises(ValueError):coefficients(terms)
