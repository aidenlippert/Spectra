"""Sparse response lift of a closed-shell determinant into QH|HF>."""
from fractions import Fraction as F
from itertools import combinations
import json, time
from experiments.marginal_symbolic import decode, hermitian
from research.molecular_collective_20260913.core import digest
from research.compact_response_20260913 import program

def apply(poly, vec, accounting=None):
    out={}
    for state,amp in vec.items():
        for word,c in poly.items():
            target=state; sign=1
            for creation,mode in reversed(word):
                occupied=(target>>mode)&1
                if occupied==creation: break
                if (target & ((1<<mode)-1)).bit_count()%2: sign=-sign
                target ^= 1<<mode
            else: out[target]=out.get(target,F(0))+sign*c*amp
    filtered={s:a for s,a in out.items() if a}
    if accounting is not None:
        accounting.update(raw_output_labels=len(out),returned_labels=len(filtered),source_term_checks=len(poly)*len(vec))
    return filtered

def hf_state(data):
    # interleaved spin orbitals: occupy the first particles/2 spatial orbitals
    n=data['particles']//2; return sum(3 << (2*i) for i in range(n))

def lift(data, alpha=None):
    m,n=data['modes'],data['particles']
    if type(m) is not int or type(n) is not int or m<4 or m%2 or n<2 or n%2 or n>m-2:
        raise ValueError('Even closed-shell N with an empty highest orbital is required')
    start=time.monotonic(); H=decode(data['hamiltonian'],data['modes'],4); hf=hf_state(data)
    if not hermitian(H) or any(sum(2*c-1 for c,_ in w) for w in H):
        raise ValueError('A Hermitian number-conserving Hamiltonian is required')
    qmask=3 << (data['modes']-2)
    if hf&qmask: raise ValueError('The HF seed must lie entirely in P')
    first_work={};second_work={}
    h0=apply(H,{hf:F(1)},first_work); e0=h0.get(hf,F(0))
    chi={s:a for s,a in h0.items() if s!=hf and (s&qmask)==qmask}
    # QH|HF> is a finite excitation list; no fixed-N basis is built.
    Hchi=apply(H,chi,second_work); m1=sum(a*a for a in chi.values()); m2=sum(a*Hchi.get(s,F(0)) for s,a in chi.items())
    if alpha is None:
        # Fixed rational grid, including zero; no optimum outside the grid
        # or even a stationary-point claim is inferred.
        candidates=[F(i,1000) for i in range(-10000,10001)]
        alpha=min(candidates,key=lambda a:(e0-2*a*m1+a*a*m2)/(1+a*a*m1))
    if type(alpha) is not F or abs(alpha)>10 or alpha.denominator>10**8:
        raise ValueError('Bounded rational lift parameter required')
    num=e0-2*alpha*m1+alpha*alpha*m2; den=1+alpha*alpha*m1
    return {'kind':'sparse_response_lift_upper_v1','fixture_sha256':digest(data),'modes':data['modes'],'particles':data['particles'],
      'hf_state':hf,'q_mask':qmask,'alpha':str(alpha),'hf_energy_Ha':str(e0),'chi_terms':len(chi),
      'chi_support':sorted(chi),'chi_coefficients':[str(chi[s]) for s in sorted(chi)],'chi_norm_squared_Ha2':str(m1),
      'chi_H_chi_Ha3':str(m2),'numerator_Ha':str(num),'denominator':str(den),'upper_Ha':str(num/den),
      'construction_seconds':time.monotonic()-start,'h0_labels':len(h0),'hchi_labels':len(Hchi),
      'union_materialized_labels':len(set(h0)|set(Hchi)),'source_term_checks':len(H)*(1+len(chi)),
      'peak_stored_labels':max(1+first_work['raw_output_labels']+len(h0),len(h0)+len(chi)+second_work['raw_output_labels']+len(Hchi)),
      'peak_label_metric':'Maximum simultaneously stored action-dictionary label entries, including raw and filtered outputs; excludes scalar/recipe serialization.',
      'states_constructed':len(set(h0)|set(Hchi)),
      'scope':'Sparse excitation support only; no fixed-N determinant enumeration.'}

def check(data, cert):
    if cert.get('kind')!='sparse_response_lift_upper_v1' or cert.get('fixture_sha256')!=digest(data):
        raise ValueError('Lift fixture binding failed')
    if type(cert.get('alpha')) is not str: raise ValueError('Exact rational alpha string required')
    fresh=lift(data,F(cert['alpha']))
    if set(cert)!=set(fresh): raise ValueError('Unknown or missing lift fields')
    for key in ('modes','particles','hf_state','q_mask','hf_energy_Ha','chi_terms','chi_support','chi_coefficients','chi_norm_squared_Ha2','chi_H_chi_Ha3','numerator_Ha','denominator','upper_Ha'):
        if fresh[key]!=cert[key]: raise ValueError('Lift exact field mismatch: '+key)
    return {'upper_Ha':cert['upper_Ha'],'chi_terms':cert['chi_terms'],'h0_labels':fresh['h0_labels'],
      'hchi_labels':fresh['hchi_labels'],'union_materialized_labels':fresh['union_materialized_labels'],
      'peak_stored_labels':fresh['peak_stored_labels'],'source_term_checks':fresh['source_term_checks'],
      'many_body_basis_enumerated':0,'replay_seconds':fresh['construction_seconds']}

def run():
    out=program.ROOT/'results/global_response_20260913/lifted_upper';out.mkdir(parents=True,exist_ok=True); rows=[]
    for name in ('h6','fresh_h6_1p6','h8','fresh_h6_1p73'):
        if name.startswith('fresh_'):
            src=program.ROOT/'results/response_consistency_20260913'/name
            if not (src/'fixture.json').exists(): src=program.ROOT/'results/global_response_20260913'/name
            data=json.loads((src/'fixture.json').read_text())
        else: data,tail,ref=program.load_case(name)
        cert=lift(data); path=out/f'{name}.json';path.write_text(json.dumps(cert,indent=2)+'\n');rows.append({'certificate':cert,'replay':check(data,cert)})
    (out/'discovery.json').write_text(json.dumps({'cases':rows},indent=2)+'\n')
    print(json.dumps([(r['certificate']['modes'],r['certificate']['upper_Ha'],r['certificate']['chi_terms']) for r in rows]))
if __name__=='__main__':run()
