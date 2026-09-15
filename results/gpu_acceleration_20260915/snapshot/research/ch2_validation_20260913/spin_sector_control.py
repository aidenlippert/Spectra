"""Exact two-electron spin decomposition; no chemical-energy claim."""
from fractions import Fraction as F
from research.ch2_validation_20260913.certify import matrix, spin_square
from research.certificate_scaling.commutator_dual_witness import psd
if __name__ == "__main__":
    raw,den,_=matrix(spin_square(4),[3,6,9,12])
    S=[[F(x,den) for x in row] for row in raw]
    assert psd(S)["rank"]==1
    assert sum(S[i][i] for i in range(4))==2
    assert [[sum(S[i][k]*S[k][j] for k in range(4)) for j in range(4)] for i in range(4)]==[[2*x for x in row] for row in S]
    print({"Ms":0,"singlet_dimension":3,"triplet_dimension":1,"exact":True})
