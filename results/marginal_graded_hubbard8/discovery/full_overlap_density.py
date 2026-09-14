"""Exact full five-site overlap diagnostic and sparse positive-projector witness."""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
from collections import Counter, defaultdict
import argparse
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import _reflection, _sector
from experiments.marginal_charged_projectors import _particlehole, _spinflip


def transformed(vector,action):
    return {action(s)[0]:a*action(s)[1] for s,a in vector.items() if a}


def canonical(vector):
    items=sorted(vector.items())
    sign=1 if items[0][1]>0 else -1
    return tuple((s,sign*a) for s,a in items)


def orbit(vector):
    """Eight projective symmetry images; verify closure on each input source."""
    generators=(_particlehole,_spinflip,lambda s:_reflection(s,6))
    vectors=[vector]
    for action in generators:vectors += [transformed(v,action) for v in list(vectors)]
    images=Counter(canonical(v) for v in vectors)
    for action in generators:
        actual=Counter()
        for image,multiplicity in images.items():actual[canonical(transformed(dict(image),action))]+=multiplicity
        if actual!=images:raise ValueError('Symmetry images are not a closed density-matrix orbit')
    return images


def partial_upper(vector,sites,left):
    """Unnormalized contiguous partial trace, upper triangle, exact integers."""
    if type(sites) is not int or not 2<=sites<=6:raise ValueError('Two through six sites required')
    if not vector or any(type(s) is not int or not 0<=s<4**sites or type(a) is not int for s,a in vector.items()):
        raise ValueError('Bounded integer Fock vector required')
    groups=defaultdict(dict);mask=4**(sites-1)-1
    for s,a in vector.items():
        retained,environment=(s&mask,s>>(2*(sites-1))) if left else (s>>2,s&3)
        groups[environment][retained]=a
    result=defaultdict(int)
    for values in groups.values():
        items=sorted(values.items())
        for i,(s,a) in enumerate(items):
            for t,b in items[i:]:result[s,t]+=a*b
    return {key:value for key,value in result.items() if value}


def reconstruct(mixture):
    if not 1<=len(mixture)<=137:raise ValueError('At most137 physical mixture sources required')
    sources=[]
    for item in mixture:
        weight=F(item['weight']);v={int(s):a for s,a in item['vector'].items() if a}
        if weight<=0 or not v or len(v)>4096 or any(type(a) is not int or abs(a)>10**9 for a in v.values()):
            raise ValueError('Positive bounded integer physical source required')
        if any(not 0<=s<4096 for s in v) or len({_sector(s,6) for s in v})!=1:
            raise ValueError('One physical spin-number sector per source required')
        norm=sum(a*a for a in v.values());sources.append((weight/(8*norm),v))
    if sum(F(item['weight']) for item in mixture)!=1:raise ValueError('Normalized mixture required')
    denominator=lcm(*(w.denominator for w,v in sources))
    if denominator.bit_length()>12000:raise ValueError('Common denominator exceeds12000-bit replay budget')
    matrices=[defaultdict(int),defaultdict(int)];orbit_images=0
    for weight,vector in sources:
        integer=weight*denominator
        if integer.denominator!=1:raise ValueError('Noninteger common-denominator coefficient')
        for image,multiplicity in orbit(vector).items():
            orbit_images+=1;v=dict(image);coefficient=int(integer)*multiplicity
            for matrix,left in zip(matrices,(True,False)):
                for key,value in partial_upper(v,6,left).items():matrix[key]+=coefficient*value
    left,right=({key:value for key,value in matrix.items() if value} for matrix in matrices)
    sectors={s:_sector(s,5) for s in range(1024)}
    for matrix in (left,right):
        if sum(value for (s,t),value in matrix.items() if s==t)!=denominator:
            raise ValueError('Partial trace lost unit trace')
        if any(sectors[s]!=sectors[t] for s,t in matrix):raise ValueError('Reduced spin-sector support changed')
    difference={key:left.get(key,0)-right.get(key,0) for key in left.keys()|right.keys()}
    difference={key:value for key,value in difference.items() if value}
    return denominator,left,right,difference,orbit_images


def entry(matrix,s,t):return matrix.get((min(s,t),max(s,t)),0)


def quadratic(matrix,vector):
    return sum(a*b*entry(matrix,s,t) for s,a in vector.items() for t,b in vector.items())


def best_sparse_projector(difference):
    best=None
    for (s,t),value in sorted(difference.items()):
        if s==t:
            nominees=[({s:1},value,1)]
        else:
            diagonal=entry(difference,s,s)+entry(difference,t,t)
            nominees=[({s:1,t:sign},diagonal+2*sign*value,2) for sign in (1,-1)]
        for vector,numerator,norm in nominees:
            if best is None or abs(numerator)*best[2]>abs(best[1])*norm:
                best=(vector,numerator,norm)
    if best is None or not best[1]:raise ValueError('No nonzero sparse projector witness')
    return best


def matrix_digest(matrix):
    digest=hashlib.sha256()
    for (s,t),value in sorted(matrix.items()):
        # Signed decimal integers with unambiguous separators; fixed upper triangle.
        digest.update(f'{s},{t}:{value}\n'.encode())
    return digest.hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    cp=args.directory/'range_two_family_limit_certificate.json';rp=args.directory/'range_two_family_limit_replay.json'
    c=json.loads(cp.read_text());r=json.loads(rp.read_text())
    if not r.get('accepted') or not r.get('three_spectator_hopping'):raise ValueError('Accepted three-spectator family required')
    for name,h in r['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=h:raise ValueError('Stale family receipt source')
    denominator,left,right,difference,images=reconstruct(c['mixture'])
    if not difference:
        files={Path(__file__).resolve(),cp.resolve(),rp.resolve()}
        for module in tuple(sys.modules.values()):
            name=getattr(module,'__file__',None)
            if name and str(Path(name).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(name).resolve())
        result={'accepted':True,'stationary_extension_refuted':False,'five_site_overlap_agrees':True,
                'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
                'scope':'Exact full five-site overlap agreement is necessary, but does not prove a stationary quantum extension.'}
        (args.output/'full_overlap_replay.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({'accepted':True,'five_site_overlap_agrees':True}),flush=True)
        return
    vector,numerator,norm=best_sparse_projector(difference)
    left_probability=F(quadratic(left,vector),denominator*norm)
    right_probability=F(quadratic(right,vector),denominator*norm)
    moment=left_probability-right_probability
    if moment!=F(numerator,denominator*norm) or not 0<=left_probability<=1 or not 0<=right_probability<=1:
        raise ValueError('Positive-projector probability check failed')
    if not moment:raise ValueError('Stationary obstruction is zero')
    diagonal_count=sum(s==t for s,t in difference)
    frobenius_square=F(sum((1 if s==t else 2)*v*v for (s,t),v in difference.items()),denominator**2)
    files={Path(__file__).resolve(),cp.resolve(),rp.resolve()}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name and str(Path(name).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(name).resolve())
    result={'accepted':True,'stationary_extension_refuted':True,'symmetrization':'Equal eight-image projective orbit under six-site particle-hole, spin flip and reflection.',
            'sites':6,'overlap_sites':5,'mixture_sources':len(c['mixture']),'unique_orbit_images_summed':images,
            'common_denominator':str(denominator),'common_denominator_bits':denominator.bit_length(),
            'left_upper_nonzeros':len(left),'right_upper_nonzeros':len(right),
            'difference_upper_nonzeros':len(difference),'difference_diagonal_nonzeros':diagonal_count,
            'difference_frobenius_square':str(frobenius_square),'difference_frobenius_square_float':float(frobenius_square),
            'left_integer_matrix_sha256':matrix_digest(left),'right_integer_matrix_sha256':matrix_digest(right),
            'difference_integer_matrix_sha256':matrix_digest(difference),
            'witness':{'five_site_vector':{str(s):a for s,a in vector.items()},'norm_square':norm,
                       'spin_sector':list(_sector(next(iter(vector)),5)),
                       'left_probability':str(left_probability),'right_probability':str(right_probability),
                       'exact_difference':str(moment),'difference_float':float(moment),
                       'six_site_projector_difference_norm_bound':1},
            'seconds':time.monotonic()-started,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
            'scope':'Complete contiguous five-site partial traces of an exact positive normalized symmetry-averaged six-site mixture. A nonzero difference of expectations of the same normalized positive projector refutes a stationary quantum extension of this particular mixture. Both embedded projectors lie between0 andI, so their difference has norm at most1. No new energy certificate, family optimum or general representability conclusion.'}
    (args.output/'full_overlap_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    (args.output/'projector_witness.json').write_text(json.dumps({'accepted':False,'vector':result['witness']['five_site_vector'],
        'source_receipt':str((args.output/'full_overlap_replay.json').resolve().relative_to(ROOT)),
        'scope':'Source for a possible new stationary projector telescope; not an energy certificate.'},indent=2)+'\n')
    print(json.dumps({'accepted':True,'nonzero_upper_entries':len(difference),'nonzero_diagonal_entries':diagonal_count,
        'witness':result['witness']['five_site_vector'],'moment':float(moment),'seconds':result['seconds']}),flush=True)


if __name__=='__main__':main()
