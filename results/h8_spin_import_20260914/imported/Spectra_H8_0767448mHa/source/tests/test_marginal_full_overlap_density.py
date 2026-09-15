"""Independent tensor contractions and exact positive-projector diagnostic checks."""
from pathlib import Path
from fractions import Fraction as F
import sys
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'results/marginal_graded_hubbard8/discovery'))
from full_overlap_density import partial_upper,reconstruct,best_sparse_projector,quadratic,entry


def test_both_contiguous_partial_traces_match_dense_tensor_contractions():
    vector={s:(s%7)-3 for s in range(64) if s.bit_count()==2 and (s%7)-3}
    dense=np.zeros(64,dtype=np.int64)
    for s,a in vector.items():dense[s]=a
    left=dense.reshape(4,16);right=dense.reshape(16,4)
    for actual,expected in [(partial_upper(vector,3,True),left.T@left),
                            (partial_upper(vector,3,False),right@right.T)]:
        assert actual=={(i,j):int(expected[i,j]) for i in range(16) for j in range(i,16) if expected[i,j]}


def test_stationary_vacuum_full_mixture_has_equal_unit_trace_marginals():
    denominator,left,right,difference,_=reconstruct([
        {'weight':'1/3','vector':{'0':1}}, {'weight':'2/3','vector':{'4095':1}}])
    assert not difference and left==right
    assert F(left[0,0],denominator)==F(1,2)==F(left[1023,1023],denominator)


def test_asymmetric_physical_source_gives_valid_projector_probabilities():
    denominator,left,right,difference,_=reconstruct([{'weight':'1','vector':{'1':1}}])
    vector,numerator,norm=best_sparse_projector(difference)
    p=F(quadratic(left,vector),denominator*norm);q=F(quadratic(right,vector),denominator*norm)
    assert 0<=p<=1 and 0<=q<=1 and p!=q
    assert p-q==F(numerator,denominator*norm)


def test_marginals_are_invariant_under_integer_source_rescaling():
    results=[reconstruct([{'weight':'1','vector':{'1':2*k,'4':3*k}}]) for k in (1,7)]
    d1,l1,r1,_,_=results[0];d2,l2,r2,_,_=results[1]
    for a,b in [(l1,l2),(r1,r2)]:
        assert {k:F(v,d1) for k,v in a.items()}=={k:F(v,d2) for k,v in b.items()}


def test_sparse_projector_uses_coherence_when_diagonal_expectations_cancel():
    difference={(0,0):2,(1,1):-2,(0,1):3}
    vector,numerator,norm=best_sparse_projector(difference)
    assert abs(F(numerator,norm))==3
    assert quadratic(difference,vector)==numerator and norm==2


@pytest.mark.parametrize('mixture',[
    [{'weight':'-1','vector':{'1':1}}],
    [{'weight':'1/2','vector':{'1':1}}],
    [{'weight':'1','vector':{'0':1,'1':1}}],
    [{'weight':'1','vector':{'4096':1}}],
])
def test_invalid_mixture_sources_are_refused(mixture):
    with pytest.raises(ValueError):reconstruct(mixture)
