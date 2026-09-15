import unittest
from fractions import Fraction as F
import itertools
import math
import numpy as np
from experiments.v7_certificate import Generator, Piece, derive_certificate, check_certificate
I=np.eye(2,dtype=complex); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]],complex); Z=np.diag([1,-1]).astype(complex)
PM={"I":I,"X":X,"Y":Y,"Z":Z}
def kron_label(s):
    a=np.array([[1]],complex)
    for c in s: a=np.kron(a,PM[c])
    return a
def dense_generator(terms,n,gamma):
    labels=[''.join(x) for x in itertools.product('IXYZ',repeat=n)]; mats=[kron_label(p) for p in labels]; d=2**n; G=np.zeros((d*d,d*d),complex)
    H=sum((float(c)*kron_label(p) for p,c in terms),np.zeros((d,d),complex))
    for j,(p,q) in enumerate(zip(labels,mats)):
        v=1j*(H@q-q@H)-gamma*sum(c!='I' for c in p)*q
        G[:,j]=[np.trace(a.conj().T@v)/d for a in mats]
    return labels,mats,G
def expm_taylor(A):
    """Independent scaling/squaring exponential (no scipy/eigenvector inverse)."""
    s=max(0,int(np.ceil(np.log2(max(1.0,np.linalg.norm(A,np.inf))))))
    B=A/(2**s); E=np.eye(A.shape[0],dtype=complex); term=E.copy()
    for k in range(1,120):
        term=term@B/k; E=E+term
        if np.linalg.norm(term,np.inf)<1e-16: break
    for _ in range(s): E=E@E
    return E
def evolve(O,labels,mats,G,t):
    d=O.shape[0]; c=np.array([np.trace(q.conj().T@O)/d for q in mats])
    return sum((x*q for x,q in zip(expm_taylor(G*t)@c,mats)),np.zeros_like(O))
class DenseCertificateTests(unittest.TestCase):
    def test_two_and_three_qubit_constant_candidates(self):
        for terms,target in [([('ZI',F(1,2)),('XX',F(1,3))],'YI'),([('ZII',F(1,2)),('XXI',F(1,3)),('IZZ',F(1,5))],'YII')]:
            n=len(target)
            for gamma in (F(0),F(1,5),F(2)):
                for t in (F(1,10),F(1)):
                    g=Generator(dict(terms),gamma,n); coeff=[{target:F(1)}]; current=coeff[0]
                    for k in range(1,13):
                        current=g.apply(current); coeff.append({p:c/F(math.factorial(k)) for p,c in current.items()})
                    pieces=[Piece(t,tuple(coeff))]; cert=derive_certificate(Generator(dict(terms),gamma,n),{target:F(1)},pieces,integration_basis='bernstein')
                    got=check_certificate(Generator(dict(terms),gamma,n),{target:F(1)},pieces,cert,F(1),expected_time=t)
                    self.assertEqual(got['status'],'certified')
                    labels,mats,G=dense_generator(terms,n,float(gamma)); exact=evolve(kron_label(target),labels,mats,G,float(t))
                    candidate=sum((float(c)*(float(t)**k)*kron_label(p) for k,op in enumerate(coeff) for p,c in op.items()),np.zeros_like(exact))
                    # At gamma=0 independently confirm the generator exponential
                    # against direct unitary conjugation.
                    if gamma == 0:
                        H=sum((float(c)*kron_label(p) for p,c in terms),np.zeros_like(exact)); U=expm_taylor(-1j*H*float(t))
                        self.assertLess(np.linalg.norm(exact-U.conj().T@kron_label(target)@U,ord=2),1e-10)
                    # The degree-12 Taylor candidate has a genuinely small bound.
                    err=np.linalg.norm(exact-candidate,ord=2)
                    self.assertLessEqual(err,float(F(got['bound']))+1e-8,(n,gamma,t,err,got['bound']))
if __name__=='__main__': unittest.main()
