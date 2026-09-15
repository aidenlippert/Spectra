from fractions import Fraction as F
from itertools import product
import pytest
from experiments.marginal_signed_charge_telescope import PATTERNS,ORBIT,coefficients,five_site_value,local_value
from experiments.marginal_projector_extendibility import replay


def direct_five(s,terms):
    q=tuple(((s//4**i)%4==3)-((s//4**i)%4==0) for i in range(5))
    value=F(0)
    for label,a in terms.items():
        pattern=tuple(map(int,label.split(',')));negative=tuple(-x for x in pattern)
        value+=a*(int(q in (pattern,negative))-int(q in (pattern[::-1],negative[::-1])))
    return value


def test_disjoint_complete_pattern_basis_symmetry():
    assert len(PATTERNS)==52 and len(ORBIT)==208
    assert len(set(ORBIT))==4*len(PATTERNS)
    for q in product((-1,0,1),repeat=5):
        if q in ORBIT:
            key,sign=ORBIT[q]
            assert ORBIT[tuple(-v for v in q)]==(key,sign)
            assert ORBIT[q[::-1]]==(key,-sign)
        else:assert q==q[::-1] or q==tuple(-v for v in q[::-1])


def test_full_fock_operator_and_psd_against_independent_pattern_sum():
    terms={key:F(i+1,53) for i,key in enumerate(PATTERNS)}
    five=[direct_five(s,terms) for s in range(1024)]
    assert all(five_site_value(s,terms)==five[s] for s in range(1024))
    values=[]
    for s in range(4096):
        assert local_value(s,terms)==five[s%1024]-five[s//4]
        q=[((s>>(2*i))&3).bit_count()-1 for i in range(6)];p=[v*v for v in q]
        old=F(q[0]*q[3]-2*q[1]*q[4]+q[2]*q[5],7)+F(p[0]*p[1]-p[3]*p[4]-p[1]*p[2]+p[4]*p[5],13)+F(p[0]*p[1]*p[2]-p[2]*p[3]*p[4]-p[1]*p[2]*p[3]+p[3]*p[4]*p[5],11)
        sparse=int(s%1024==102)-int(s%1024==612)-int(s//4==102)+int(s//4==612)
        values.append(five[s%1024]-five[s//4]+old+sparse)
    c={'kind':'hubbard_projector_extension_v10','chain_sites':10,'target':{'U':'0','t':'0','V':'0','W':'0'},'local_window':{'kind':'local_hubbard_range2_block_v1','sites':6,'U':'0','t':'0','V':'0'},'vector':{0x999:1,0x666:1},'windows':2,'projector_sum_ceiling':'2','penalty':'0','penalized_lower':str(min(values)),'joint':{'vector':{62:1,3008:1},'windows':2,'ratio':'1/2','projector_sum_ceiling':'2','penalty':'0'},'telescoping_diagonal':{102:1,612:-1},'quadratic_charge_telescope':{'0,3':'1/7'},'charge_square_pair_telescope':{'0,1':'1/13'},'higher_charge_indicator_telescope':{'0,1,2':'1/11'},'signed_charge_telescope':{key:str(v) for key,v in terms.items()}}
    result=replay(c)
    assert result['local_sum_dimensions']==4096 and result['local_maximum_psd_dimension']==200
    assert result['signed_charge_components']==52
    c['penalized_lower']=str(min(values)+F(1,53053))
    with pytest.raises(ValueError):replay(c)
    c['kind']='hubbard_projector_extension_v9'
    with pytest.raises(ValueError,match='require v10'):replay(c)


@pytest.mark.parametrize('value',[{}, {'0,0,0,0,0':1},{'1,0,0,0,0':1},{next(iter(PATTERNS)):0},{next(iter(PATTERNS)):0.1},{next(iter(PATTERNS)):True},{next(iter(PATTERNS)):'1/1000001'}])
def test_noncanonical_or_inexact_pattern_coefficients_refused(value):
    with pytest.raises(ValueError):coefficients(value)
