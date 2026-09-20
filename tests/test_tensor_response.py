from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json,subprocess,sys,tempfile,unittest
import numpy as np
from spectra_tensor import tensor as t,physics as p,exact,io
from spectra_tensor.solver import screen
from spectra_tensor.variational import sweep_solve,enrich

def model(n=4):
    return dict(kind='spin_independent_hubbard_v1',sites=n,particles=n,U=['7/3']*n,edges=[[i,i+1,'2/3'] for i in range(n-1)]+([[0,n-1,'-1/5']] if n>2 else []),mode_order='up_then_down',units='t')
def basis(n):return sorted(sum(1<<(2*i) for i in up)+sum(1<<(2*i+1) for i in down) for up in combinations(range(n),n//2) for down in combinations(range(n),n//2))
def physical(state,n):return [(state>>(2*i))&3 for i in range(n)]
def oracle(spec):
    n=spec['sites'];states=basis(n);lookup={s:i for i,s in enumerate(states)};H=np.zeros((len(states),len(states)))
    for col,state in enumerate(states):
        H[col,col]=sum(float(F(u)) for i,u in enumerate(spec['U']) if ((state>>(2*i))&3)==3)
        for a,b,v in spec['edges']:
            for spin in (0,1):
                for create,destroy in ((2*a+spin,2*b+spin),(2*b+spin,2*a+spin)):
                    if not (state>>destroy)&1:continue
                    sign=(-1)**((state&((1<<destroy)-1)).bit_count());mid=state^(1<<destroy)
                    if (mid>>create)&1:continue
                    sign*=(-1)**((mid&((1<<create)-1)).bit_count());H[lookup[mid^(1<<create)],col]-=float(F(v))*sign
    return H,states

def vector(x,states):return np.array([t.amplitude(x,physical(s,len(x.cores))) for s in states])
def proposal(spec):
    b=t.product(p.source(spec));hb=t.apply(p.compile_hubbard(spec),b)
    return t.compress(t.add(b,hb,alpha=.13-.2j,beta=.04+.015j),32,0.)[0]

class TensorResponseTests(unittest.TestCase):
    def test_independent_fermion_action_every_column(self):
        spec=model();H,states=oracle(spec);mpo=p.compile_hubbard(spec);imp=exact.integer_hamiltonian(spec)
        for j,state in enumerate(states):
            b=t.product(physical(state,4));np.testing.assert_allclose(vector(t.apply(mpo,b),states),H[:,j],atol=1e-14)
            action=exact.apply_integer(imp,io.export_program(b,8));D=action['local_denominator']**4*imp['den'];vals=[exact.source_amplitude(action,physical(s,4))[0]/D for s in states]
            np.testing.assert_allclose(vals,H[:,j],atol=1e-14)
    def test_streaming_and_canonicalization(self):
        spec=model();x=proposal(spec);mpo=p.compile_hubbard(spec);_,states=oracle(spec)
        y,_=t.apply_compressed(mpo,x,64,0.);np.testing.assert_allclose(vector(y,states),vector(t.apply(mpo,x),states),atol=2e-13)
        np.testing.assert_allclose(vector(t.right_canonical(x),states),vector(x,states),atol=2e-13)
        rounded,_=t.compress(y,3);rounded.validate();exact.validate_program(io.export_program(rounded),4)
    def test_exact_python_gmp_and_dense_containment(self):
        spec=model();z=2+1.5j;x=proposal(spec);c=io.candidate(spec,[z],[x]);r=io.request(spec,[z]);a=exact.verify(c,r,'python');b=exact.verify(c,r,'gmp');self.assertEqual(a['queries'],b['queries'])
        H,states=oracle(spec);source=vector(t.product(p.source(spec)),states);truth=source@np.linalg.solve(z*np.eye(len(H))-H,source);row=b['queries'][0];center=complex(float(F(row['center_real'])),float(F(row['center_imag'])))
        self.assertLessEqual(abs(truth-center),float(F(row['radius']))+1e-13);self.assertAlmostEqual(screen(spec,x,z)['radius'],float(F(row['radius'])),places=11)
    def test_lazy_action(self):
        spec=model();x=io.export_program(proposal(spec));h=exact.integer_hamiltonian(spec);a=exact.apply_integer(h,x);b=exact.apply_integer(h,x,lazy=False)
        self.assertEqual(list(a['cores']),b['cores']);self.assertEqual(exact.exact_overlap(a,a),exact.exact_overlap(b,b))
    def test_tampering_rejected(self):
        spec=model();z=2+1.5j;c=io.candidate(spec,[z],[proposal(spec)]);r=io.request(spec,[z]);pairs=[]
        bad=deepcopy(r);bad['model']['U'][0]='3';pairs.append((c,bad))
        bad=deepcopy(r);bad['source'][0]=2;pairs.append((c,bad))
        bad=deepcopy(r);bad['frequencies'][0]['eta']='0';pairs.append((c,bad))
        bad=deepcopy(r);bad['frequencies'][0]['omega']='7';pairs.append((c,bad))
        bad=deepcopy(c);bad['queries'][0]['program']['cores'][0].append(bad['queries'][0]['program']['cores'][0][0]);pairs.append((bad,r))
        bad=deepcopy(c);bad['queries'][0]['program']['cores'][0][0][1]=3;pairs.append((bad,r))
        bad=deepcopy(r);bad['model']['U'][0]=2.;pairs.append((c,bad))
        for a,b in pairs:
            with self.assertRaises((ValueError,AssertionError)):exact.verify(a,b,'python')
    def test_zero_proposal_wide_bound(self):
        spec=model();z=2+1.5j;x=t.product(p.source(spec)).scaled(0);r=exact.verify(io.candidate(spec,[z],[x]),io.request(spec,[z]),'python');self.assertEqual(F(r['queries'][0]['radius']),F(2,3))
    def test_local_refinement(self):
        spec=model();z=2+1.5j;x,_=enrich(spec,z,proposal(spec),32);x,_=sweep_solve(spec,z,x,sweeps=2,target=1e-9,verbose=False)
        r=exact.verify(io.candidate(spec,[z],[x]),io.request(spec,[z]),'python');self.assertLess(F(r['queries'][0]['radius']),F(1,10**7))
    def test_zero_hamiltonian(self):
        from spectra_tensor.solver import cocg
        spec=model();spec['U']=['0']*4;spec['edges']=[];z=2+1.5j;x,_=cocg(spec,z,verbose=False)
        r=exact.verify(io.candidate(spec,[z],[x]),io.request(spec,[z]),'python');self.assertLess(F(r['queries'][0]['radius']),F(1,10**12))
    def test_public_api_separate_check_and_fraction(self):
        from spectra_tensor.api import solve
        spec=model();r=io.request(spec,[complex(1/3,1.5)]);r['frequencies'][0]['omega']='1/3'
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory)/'run';result=solve(r,d,bonds=(16,32),seed_steps=4,sweeps=2,verbose=False);self.assertEqual(result['status'],'target_met')
            receipt=json.loads((d/result['queries'][0]['receipt']).read_text());self.assertEqual(receipt['queries'][0]['omega'],'1/3');self.assertEqual(receipt['numerical_packages_loaded'],[])
            with self.assertRaises(FileExistsError):solve(r,d,verbose=False)
    def test_reference_matches_independent_oracle(self):
        from spectra_tensor.reference import solve_reference
        spec=model();z=2+1.5j;r=solve_reference(spec,z,target=1e-9);H,states=oracle(spec);b=vector(t.product(p.source(spec)),states);truth=b@np.linalg.solve(z*np.eye(len(H))-H,b)
        self.assertLess(abs(complex(r['center_real'],r['center_imag'])-truth),r['radius']+1e-13)
        limit=solve_reference(model(20),z)
        self.assertEqual(limit['status'],'resource_limit_preflight');self.assertEqual(limit['single_complex_vector_bytes'],546156472576)
    def test_insufficient_rank_is_not_accepted(self):
        from spectra_tensor.api import solve
        with tempfile.TemporaryDirectory() as directory:
            r=solve(io.request(model(),[2+1.5j]),Path(directory)/'run',bonds=(1,),seed_steps=2,sweeps=1,verbose=False)
            self.assertEqual(r['status'],'target_not_met');self.assertGreater(F(r['max_radius']),F(1,1000))
if __name__=='__main__':unittest.main()
