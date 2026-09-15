import json
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from experiments.marginal_symbolic import add, adj, canonical, mono, product, scale, encode
from research.certificate_scaling.projected_product_parent import recognize, edge_factor


def hamiltonian(q,g,c=0):
    # Independently use physical n_i*n_j and hopping expansion, not B squares.
    terms=[mono((),F(c))]
    for i,w in enumerate(g):
        ni=mono(((1,i),(0,i)));nj=mono(((1,i+1),(0,i+1)))
        terms.extend([scale(ni,w*q[i+1]**2),scale(nj,w*q[i]**2),
            scale(product(ni,nj),-w*(q[i]**2+q[i+1]**2)),
            mono(((1,i),(0,i+1)),-w*q[i]*q[i+1]),
            mono(((1,i+1),(0,i)),-w*q[i]*q[i+1])])
    return canonical(add(*terms))


def apply(poly,vector):
    # Independent occupation-bit CAR action, used only in small tests.
    result={}
    for bits,amplitude in vector.items():
        for word,coefficient in poly.items():
            state=bits;value=amplitude*coefficient
            for create,i in reversed(word):
                if ((state>>i)&1)==create:value=0;break
                if (state&((1<<i)-1)).bit_count()%2:value=-value
                state^=1<<i
            if value:result[state]=result.get(state,F(0))+value
    return {s:v for s,v in result.items() if v}


class ParentTests(unittest.TestCase):
    def test_inferred_permuted_signed_path_and_transported_fock_state(self):
        q=tuple(map(F,(1,2,1,3)));g=tuple(map(F,(2,5,7)))
        order=(2,0,3,1);signs=(1,-1,-1,1);source=hamiltonian(q,g,4)
        transformed={}
        for w,c in source.items():
            for _,i in w:c*=signs[i]
            transformed[tuple((create,order[i]) for create,i in w)]=c
        h=canonical(transformed);r=recognize(h,4,2)
        self.assertTrue(r['accepted']);self.assertEqual(r['lower'],'4')
        self.assertTrue(r['signed_permutation_round_trip_checked'])
        # Construct transported amplitudes independently, including wedge sign.
        v={}
        for i,j in combinations(range(4),2):
            v[(1<<order[i])|(1<<order[j])]=q[i]*q[j]*signs[i]*signs[j]*(1 if order[i]<order[j] else -1)
        self.assertEqual(apply(h,v),{s:4*x for s,x in v.items()})

    def test_treewidth_one_star_has_no_common_null_state(self):
        # Adversarial boundary: the unchanged ansatz fails on a branching tree.
        squares=[]
        for j in (1,2,3):
            B=add(product(mono(((0,0),)),add(mono(()),mono(((1,j),(0,j)),-1))),
                  scale(product(mono(((0,j),)),add(mono(()),mono(((1,0),(0,0)),-1))),-1))
            squares.append(product(canonical(adj(B)),B))
        h=add(*squares);states=[sum(1<<i for i in s) for s in combinations(range(4),2)]
        cols=[apply(h,{s:F(1)}) for s in states]
        matrix=[[col.get(s,F(0)) for col in cols] for s in states]
        work=[row[:] for row in matrix];det=F(1)
        for i in range(6):
            pivot=next(k for k in range(i,6) if work[k][i])
            if pivot!=i:work[i],work[pivot]=work[pivot],work[i];det=-det
            value=work[i][i];det*=value
            for k in range(i+1,6):
                ratio=work[k][i]/value
                work[k]=[x-ratio*y for x,y in zip(work[k],work[i])]
        self.assertEqual(det,4)
        shifted=add(h,mono((),-2));cube=product(product(shifted,shifted),shifted)
        identity=add(cube,scale(shifted,-3))
        for s in states:self.assertEqual(apply(identity,{s:F(1)}),{})
        self.assertFalse(recognize(h,4,2)['accepted'])

    def test_nonuniform_fixed_two_particle_state_independent_fock(self):
        q=tuple(map(F,(1,2,1,3)));g=tuple(map(F,(2,5,7)));h=hamiltonian(q,g,4)
        r=recognize(h,4,2);self.assertTrue(r['accepted']);self.assertEqual(r['lower'],'4')
        self.assertEqual(r['q'],['1','2','1','3'])
        v={sum(1<<i for i in inds):q[inds[0]]*q[inds[1]] for inds in combinations(range(4),2)}
        self.assertEqual(F(r['normalization']),sum(x*x for x in v.values()))
        self.assertEqual(apply(h,v),{s:4*x for s,x in v.items()})
        for i in range(3):self.assertEqual(apply(edge_factor(i,q[i],q[i+1]),v),{})
        self.assertEqual(r['dp_updates'],8)
        # A fixed-N pure Gaussian state is a Slater determinant. Its two-form
        # obeys the Pluecker relation; this parent state's two-form does not.
        p=lambda i,j:v[(1<<i)|(1<<j)]
        self.assertNotEqual(p(0,1)*p(2,3)-p(0,2)*p(1,3)+p(0,3)*p(1,2),0)

    def test_all_number_sectors_and_rational_weights(self):
        q=(F(1),F(2,3),F(5,7));g=(F(3,5),F(7,11));h=hamiltonian(q,g,-2)
        for n in range(4):
            r=recognize(h,3,n);self.assertTrue(r['accepted']);self.assertEqual(r['upper'],'-2')
        r=recognize(hamiltonian(tuple(map(F,(1,2,3))),tuple(map(F,(5,7)))),3,2)
        self.assertEqual(r['normalization'],'49')

    def test_wrong_density_sign_and_mutations_refuse(self):
        h=hamiltonian(tuple(map(F,(1,2,1,3))),tuple(map(F,(2,5,7))))
        for w in [((1,0),(1,1),(0,0),(0,1)),((1,0),(0,1)),((1,0),(0,0)),
                  ((1,0),(0,3)),((1,3),(0,3)),((1,0),(1,1))]:
            bad=dict(h);bad[w]=bad.get(w,F(0))+1
            self.assertFalse(recognize(bad,4,2)['accepted'])
        bad={w:(-v if len(w)==4 else v) for w,v in h.items()}
        self.assertFalse(recognize(bad,4,2)['accepted'])
        for m,n in [(4,-1),(4,5),(True,1),(1,1)]:
            with self.assertRaises(ValueError):recognize(h,m,n)
        for bad in [{():1.0},{():True},{((1,8),):F(1)}]:
            with self.assertRaises(ValueError):recognize(bad,4,2)

    def test_cli_python_without_site_packages(self):
        h=hamiltonian(tuple(map(F,(1,2,3))),tuple(map(F,(5,7))))
        with tempfile.TemporaryDirectory() as d:
            fixture=Path(d)/'fixture.json';out=Path(d)/'out.json'
            fixture.write_text(json.dumps({'modes':3,'particles':2,'hamiltonian':encode(h)}))
            subprocess.run([sys.executable,'-S','-m','research.certificate_scaling.projected_product_parent',
                '--fixture',str(fixture),'--out',str(out)],check=True,capture_output=True)
            self.assertTrue(json.loads(out.read_text())['accepted'])

if __name__=='__main__':unittest.main()
