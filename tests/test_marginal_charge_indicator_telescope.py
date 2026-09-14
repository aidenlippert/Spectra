from fractions import Fraction as F
import pytest

from experiments.marginal_charge_indicator_telescope import LABELS, ALL_LABELS, coefficients, five_site_value, local_value
from experiments.marginal_projector_extendibility import replay


def direct(s,terms):
    p=[int((s//4**i)%4 in (0,3)) for i in range(6)]
    result=F(0)
    for key,a in terms.items():
        sites=list(map(int,key.split(',')))
        def term(offset,reflected):
            value=1
            for i in sites:value*=p[offset+(4-i if reflected else i)]
            return value
        result+=a*(term(0,False)-term(0,True)-term(1,False)+term(1,True))
    return result


def test_all_higher_indicators_against_full_fock_psd_and_refusal():
    assert len(LABELS)==6
    terms={key:F(i+1,11) for i,key in enumerate(LABELS)}
    values=[]
    for s in range(4096):
        assert local_value(s,terms)==direct(s,terms)
        q=[((s>>(2*i))&3).bit_count()-1 for i in range(6)]
        p=[v*v for v in q]
        old=F(q[0]*q[3]-2*q[1]*q[4]+q[2]*q[5],7)+F(p[0]*p[1]-p[3]*p[4]-p[1]*p[2]+p[4]*p[5],13)
        sparse=int(s%1024==102)-int(s%1024==612)-int(s//4==102)+int(s//4==612)
        values.append(direct(s,terms)+old+sparse)
    c={'kind':'hubbard_projector_extension_v9','chain_sites':10,
       'target':{'U':'0','t':'0','V':'0','W':'0'},
       'local_window':{'kind':'local_hubbard_range2_block_v1','sites':6,'U':'0','t':'0','V':'0'},
       'vector':{0x999:1,0x666:1},'windows':2,'projector_sum_ceiling':'2','penalty':'0','penalized_lower':str(min(values)),
       'joint':{'vector':{62:1,3008:1},'windows':2,'ratio':'1/2','projector_sum_ceiling':'2','penalty':'0'},
       'telescoping_diagonal':{102:1,612:-1},'quadratic_charge_telescope':{'0,3':'1/7'},
       'charge_square_pair_telescope':{'0,1':'1/13'},
       'higher_charge_indicator_telescope':{key:str(value) for key,value in terms.items()}}
    result=replay(c)
    assert result['local_sum_dimensions']==4096
    assert result['local_maximum_psd_dimension']==200
    assert result['higher_indicator_components']==6
    c['penalized_lower']=str(min(values)+F(1,1001))
    with pytest.raises(ValueError):replay(c)
    c['kind']='hubbard_projector_extension_v8'
    with pytest.raises(ValueError,match='require v9'):replay(c)


def test_complete_binary_basis_symmetries_and_periodic_cancellation():
    assert len(ALL_LABELS)==12
    for key in ALL_LABELS:
        for s in range(1024):
            reflected=sum(((s>>(2*i))&3)<<(2*(4-i)) for i in range(5))
            assert five_site_value(reflected,{key:F(1)})==-five_site_value(s,{key:F(1)})
            assert five_site_value(s^1023,{key:F(1)})==five_site_value(s,{key:F(1)})
        for mask in range(256):
            # Choose empty (indicator 1) or singly occupied (indicator 0).
            cells=[0 if mask>>i&1 else 1 for i in range(8)]
            assert sum(direct(sum(cells[(start+i)%8]*4**i for i in range(6)),{key:F(1)}) for start in range(8))==0


@pytest.mark.parametrize('terms',[{}, {'0,1':1},{'0,1,4':0},{'4,3,2':1},{'0,0,1':1},{'0,1,2,3,4':1},{'0,1,2':0.1},{'0,1,2':True},{'0,1,2':'1/1000001'}])
def test_malformed_indicator_terms_refused(terms):
    with pytest.raises(ValueError):coefficients(terms)
