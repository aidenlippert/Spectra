"""Exact signed-orbit dimension of the real five-site stationary test space."""
from collections import defaultdict
from pathlib import Path
import hashlib,json,sys
from full_overlap_telescope import ROOT,hole5,spin5
from experiments.marginal_local_hubbard_block import _sector,_reflection

def main():
    sectors=defaultdict(list)
    for s in range(1024):sectors[_sector(s,5)].append(s)
    pairs={(s,t) for group in sectors.values() for i,s in enumerate(group) for t in group[i:]}
    total=len(pairs);unused=set(pairs);basis=[];forbidden=0;orbit_sizes=defaultdict(int)
    generators=[(hole5,1),(spin5,1),(lambda s:_reflection(s,5),-1)]
    while unused:
        seed=min(unused);orbit={seed:1};queue=[seed];consistent=True
        while queue:
            pair=queue.pop();coefficient=orbit[pair]
            for action,character in generators:
                s,a=action(pair[0]);t,b=action(pair[1]);target=(min(s,t),max(s,t));value=coefficient*a*b*character
                if target not in pairs:raise ValueError('Signed orbit leaves Hermitian spin-number space')
                if target in orbit:
                    if orbit[target]!=value:consistent=False
                else:orbit[target]=value;queue.append(target)
        if not set(orbit)<=unused or len(orbit)>8:raise ValueError('Nonpartitioning or unbounded symmetry orbit')
        unused.difference_update(orbit);orbit_sizes[len(orbit)]+=1
        if consistent:
            for action,character in generators:
                transformed={}
                for (r,c),value in orbit.items():
                    s,a=action(r);t,b=action(c);transformed[min(s,t),max(s,t)]=a*b*value
                if transformed!={key:character*value for key,value in orbit.items()}:
                    raise ValueError('Signed orbit does not have the target symmetry')
            basis.append([[r,c,a] for (r,c),a in sorted(orbit.items())])
        else:forbidden+=1
    diagonal=sum(all(r==c for r,c,a in orbit) for orbit in basis)
    digest=hashlib.sha256(json.dumps(basis,separators=(',',':')).encode()).hexdigest()
    files={Path(__file__).resolve(),Path(sys.modules['full_overlap_telescope'].__file__).resolve(),Path(sys.modules['full_overlap_density'].__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name and str(Path(name).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(name).resolve())
    result={'accepted':True,'real_spin_number_preserving_dimension':total,
            'particle_hole_even_spin_flip_even_reflection_odd_dimension':len(basis),
            'diagonal_dimension':diagonal,'offdiagonal_dimension':len(basis)-diagonal,
            'forbidden_signed_orbits':forbidden,'orbit_size_counts':dict(orbit_sizes),
            'basis_sha256':digest,
            'dimension_proof':'Real Hermitian spin-number-preserving matrix units partition into signed permutation orbits under H,S,R. A consistent character orbit contributes exactly one independent vector, with disjoint support from all other orbits. A sign-conflicted orbit contributes zero. Every matrix unit is consumed exactly once and every retained orbit satisfies H+,S+,R- exactly.',
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
            'scope':'Exact dimension census of real five-site operators under the stated symmetries, not a positivity, representability or runtime-scalability result. No claim that all these constraints are integrated in the current energy/family implementation.'}
    out=ROOT/'results/marginal_graded_hubbard8/spin_word/overlap_space.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ('accepted','real_spin_number_preserving_dimension','particle_hole_even_spin_flip_even_reflection_odd_dimension','diagonal_dimension','offdiagonal_dimension')}))

if __name__=='__main__':main()
