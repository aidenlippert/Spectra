"""Identify exact singlet-annihilating channels in two retained dictionaries."""
from fractions import Fraction
import hashlib
import json
import numpy as np
from experiments.marginal_symbolic import add,canonical,mono,product
from research.acceptance_channels_20260915.campaign import OUT,dump
from research.acceptance_channels_20260915.channel_membership import coordinates


def channels(spatial,creation):
    raising=add(*(mono(((1,2*p),(0,2*p+1))) for p in range(spatial)))
    result=[]
    for p in range(spatial):
        operator=mono(((1,2*p),)) if creation else mono(((0,2*p+1),))
        left,right=product(raising,operator),product(operator,raising)
        if canonical(left)!=canonical(right):raise ValueError('Required exact spin commutator failed')
        result.append(canonical(left))
    return result


def run():
    prepared=OUT/'h12_cached/prepared'
    frame=json.loads((prepared/'frame.json').read_text())
    bases=np.load(prepared/'bases.npz')
    groups=[]
    for index,creation in ((22,False),(24,True)):
        block=frame['blocks'][index]
        if block['kind']!='spin_three_half_highest':raise ValueError('Unexpected block kind')
        group=frame['groups'][block['physical_group']]
        words=group['words']
        if not np.array_equal(bases[f'V_{index}'],np.eye(len(words))):
            raise ValueError('An exact identity retained basis is required')
        vectors=[]
        for polynomial in channels(12,creation):
            v=coordinates(polynomial,group)
            if v is None:raise ValueError('A required null channel is absent')
            recovered=add(*(mono(tuple(map(tuple,w)),c) for w,c in zip(words,v) if c))
            if canonical(recovered)!=polynomial:raise ValueError('Coordinates fail exact reconstruction')
            vectors.append(v)
        for i,v in enumerate(vectors):
            for j,w in enumerate(vectors):
                if sum(a*b for a,b in zip(v,w))!=Fraction(11 if i==j else 0):
                    raise ValueError('Exact independence check failed')
        groups.append({'block':index,'dimension':len(words),'independent_null_vectors':12,
            'operator':'S_plus a_p_alpha_dagger' if creation else 'S_plus a_p_beta',
            'euclidean_Gram':'11 I_12 exactly',
            'coordinates':[{str(i):str(v) for i,v in enumerate(row) if v} for row in vectors]})
    result={'scope':'Exact singlet-annihilating operators in existing H12 dictionaries; not a repaired dual or new energy certificate.',
        'argument':'S_plus annihilates every singlet. The exact CAR calculation verifies that it commutes with a_p_beta and a_p_alpha_dagger. Thus each displayed cubic channel annihilates every singlet. Its vector is in the kernel of every physical singlet Gram matrix for the corresponding dictionary.',
        'groups':groups,'exact_commutators_checked':True,'exact_dictionary_coordinates_checked':True,
        'independent_null_vectors_total':24,'candidate_dual_repaired':False,'family_obstruction_proved':False,
        'remaining_trace_kernel_blocks':'The two 588-dimensional blocks and other trace kernels are not repaired by this receipt.',
        'source_sha256':{name:hashlib.sha256((prepared/name).read_bytes()).hexdigest() for name in ('frame.json','bases.npz')}}
    dump(OUT/'singlet_null_channels.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('groups','source_sha256')}),flush=True)


if __name__=='__main__':run()
