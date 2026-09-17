from fractions import Fraction as F
import json
from pathlib import Path
import unittest
from research.local_response_20260916.projector_correction import grouped_overlap,norm2,verify_correction,profile_information_ceiling
from research.constructive_response_20260916.composition_exact import is_psd


class ProjectorTests(unittest.TestCase):
    def test_group_overlap_independent_small_tensor(self):
        # Embed four logical qubits into fixed-charge two-site rungs.
        logical=[[1,2,-1,3],[2,1,0,1],[1,1,2,-1]]
        vectors=[]
        for v in logical:
            w=[0]*256
            for a in range(2):
                for b in range(2):w[(1,4)[a]+16*(1,4)[b]]=v[a+2*b]
            vectors.append(w)
        p=grouped_overlap(vectors)
        psi=[logical[0][x&3]*logical[2][x>>2] for x in range(16)]
        norm=norm2(psi);midnorm=norm2(logical[1]);actual=F(0)
        for x in range(16):
            for y in range(16):
                if x&9 == y&9:
                    actual+=F(psi[x]*psi[y]*logical[1][(x>>1)&3]*logical[1][(y>>1)&3],norm*midnorm)
        self.assertEqual(p,actual)
        # Independently instantiate the full 16-dimensional projection sum.
        g=[F(1,2),F(4,5),F(3,4)];gamma=F(1,1000)
        self.assertGreaterEqual((min(g[0],g[2])-gamma)*(g[1]-gamma),min(g[0],g[2])*g[1]*p)
        matrix=[[F(0) for _ in range(16)] for _ in range(16)]
        for x in range(16):
            for y in range(16):
                value=(sum(g)-gamma) if x==y else 0
                for i in range(3):
                    mask=3<<i
                    if x&~mask == y&~mask:
                        value-=g[i]*F(logical[i][(x>>i)&3]*logical[i][(y>>i)&3],norm2(logical[i]))
                matrix[x][y]=value
        self.assertTrue(is_psd(matrix))

    def test_actual_certificate_and_invalid_profile(self):
        root=Path(__file__).resolve().parents[2]
        cert=json.loads((root/'results/local_response_20260916/boundary/certificate.json').read_text())
        correction=json.loads((root/'results/local_response_20260916/projector/correction.json').read_text())
        rec=verify_correction(cert,correction)
        self.assertGreater(F(rec['physical_energy_minus_old_family_ceiling_at_least_t']),F(1,4))
        correction['local_profiles'][0]['gap_over_t']='10'
        with self.assertRaises(ValueError):verify_correction(cert,correction)

    def test_profile_counterexample_by_sparse_physical_state(self):
        root=Path(__file__).resolve().parents[2]
        correction=json.loads((root/'results/local_response_20260916/projector/correction.json').read_text())
        vectors=[entry['integer_vector'] for entry in correction['local_profiles']]
        gaps=[F(entry['gap_over_t']) for entry in correction['local_profiles']]
        info=profile_information_ceiling(vectors,gaps)
        a,d=info['external_local_labels']
        psi={a+(center<<4)+(d<<12):x for center,x in enumerate(vectors[1]) if x}
        norm=sum(x*x for x in psi.values())
        energy=F(0)
        for bi,(v,g) in enumerate(zip(vectors,gaps)):
            shift=4*bi;mask=255<<shift
            overlap=sum(F(x*y*v[(sx>>shift)&255]*v[(sy>>shift)&255],norm*norm2(v))
                        for sx,x in psi.items() for sy,y in psi.items() if sx&~mask==sy&~mask)
            energy+=g*(1-overlap)
        self.assertEqual(energy,F(info['profile_model_upper_over_t']))
        self.assertLessEqual(energy,gaps[0]+gaps[2])
        self.assertFalse(info['physical_H_upper_claimed'])

if __name__=='__main__':unittest.main()
