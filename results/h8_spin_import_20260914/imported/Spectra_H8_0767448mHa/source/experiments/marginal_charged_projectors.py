"""Exact bounded Gram certificates for four charged six-site projectors.

Sources occupy distinct spin-number sectors, so their normalized rank-one
projectors are orthogonal on a single window. Every exterior configuration
is retained when several windows overlap.
"""
from fractions import Fraction as F
from experiments.marginal_local_hubbard_block import _exact, _reflection, _sector
from experiments.marginal_projector_extendibility import _psd, _projector_vector


def _spinflip(s):
    modes = [m ^ 1 for m in range(12) if s >> m & 1]
    sign = (-1)**sum(a > b for i, a in enumerate(modes) for b in modes[i+1:])
    return sum(1 << m for m in modes), sign


def _particlehole(s):
    return 4095 ^ s, (-1)**sum(m+m//2 for m in range(12) if s >> m & 1)


def _transform(vector, action):
    out = {}
    for s, a in vector.items():
        target, sign = action(s)
        if target in out:
            raise ValueError('Source transformation must be bijective')
        out[target] = sign*a
    return out


def charged_vectors(source):
    if type(source) is not dict or not 1 <= len(source) <= 300:
        raise ValueError('Bounded charged source required')
    vector = {}
    for label, a in source.items():
        if type(label) not in (str, int):
            raise ValueError('Integer determinant label required')
        s = int(label)
        if s in vector or not 0 <= s < 4096 or _sector(s, 6) != (2, 3):
            raise ValueError('Charged source must have spin numbers (2,3)')
        if type(a) is not int or abs(a) > 10**9:
            raise ValueError('Bounded integer charged amplitudes required')
        vector[s] = a
    vector = {s: a for s, a in vector.items() if a}
    norm = sum(a*a for a in vector.values())
    if not norm:
        raise ValueError('Nonzero charged source required')
    vectors = [vector, _transform(vector, _spinflip), _transform(vector, _particlehole)]
    vectors.append(_transform(vectors[1], _particlehole))
    for v, sector in zip(vectors, [(2,3), (3,2), (4,3), (3,4)]):
        if {_sector(s, 6) for s in v} != {sector} or sum(a*a for a in v.values()) != norm:
            raise ValueError('Charged transformation lost spin or norm')
        reflected = _transform(v, lambda s: _reflection(s, 6))
        if not any(reflected == {s: p*a for s, a in v.items()} for p in (-1, 1)):
            raise ValueError('Charged projectors must preserve reflection sectors')
        if _transform(_transform(v, _particlehole), _particlehole) != v:
            raise ValueError('Particle-hole transformation is not an involution')
    return vectors, norm


def _family_grams(vectors, windows):
    if type(windows) is not int or windows not in (2, 3, 4):
        raise ValueError('Charged window count must be 2, 3 or 4')
    norms = [sum(a*a for a in v.values()) for v in vectors]
    groups = {}
    for offset in range(windows):
        for index, vector in enumerate(vectors):
            for environment in range(4**(windows-1)):
                left = environment & ((1 << (2*offset))-1)
                right = environment >> (2*offset)
                column = {left | (s << (2*offset)) | (right << (2*(offset+6))): a
                          for s, a in vector.items()}
                sectors = {_sector(s, windows+5) for s in column}
                if len(sectors) != 1 or sum(a*a for a in column.values()) != norms[index]:
                    raise ValueError('Invalid charged isometry embedding')
                groups.setdefault(next(iter(sectors)), []).append((column,index))
    grams = {}
    for key, columns in sorted(groups.items()):
        grams[key] = {'gram': [[sum(a*right.get(s,0) for s,a in left.items())
                               for right,_ in columns] for left,_ in columns],
                      'norms': [norms[index] for _,index in columns],
                      'sources': [index for _,index in columns]}
    if sum(len(g['gram']) for g in grams.values()) != len(vectors)*windows*4**(windows-1):
        raise ValueError('Incomplete charged Gram coverage')
    return grams


def charged_overlap_grams(source, windows):
    vectors, norm = charged_vectors(source)
    groups = _family_grams(vectors, windows)
    grams = {key: data['gram'] for key,data in groups.items()}
    return norm, grams


def joint_overlap_grams(half_source, charged_source, windows):
    half,_ = _projector_vector(half_source,6)
    half = {s:a for s,a in half.items() if a}
    if {_sector(s,6) for s in half} != {(3,3)}:
        raise ValueError('Joint half-filled source must have spin numbers (3,3)')
    charged,_ = charged_vectors(charged_source)
    return _family_grams([half]+charged,windows)


def joint_projector_bound(half_source, charged_source, windows, ratio, ceiling,proof=None):
    ratio,ceiling = _exact(ratio,10**9),_exact(ceiling,10**9)
    if type(windows) is not int or windows not in (2,3,4) or ratio<=0 or not max(1,ratio)<=ceiling<=windows*max(1,ratio):
        raise ValueError('Invalid joint ratio, ceiling or window count')
    groups=joint_overlap_grams(half_source,charged_source,windows);sectors=[]
    for key,data in groups.items():
        # C D C^T <= B I iff B D^{-1} - C^T C >= 0, for D>0.
        diagonal=[ceiling*n/(1 if index==0 else ratio) for n,index in zip(data['norms'],data['sources'])]
        matrix=[[(diagonal[i] if i==j else F(0))-x for j,x in enumerate(row)]
                for i,row in enumerate(data['gram'])]
        sectors.append({'sector':list(key),'psd':_psd(matrix,proof)})
    return {'accepted':True,'windows':windows,'support_sites':windows+5,
            'charged_ratio':str(ratio),'ceiling':str(ceiling),'average_ceiling':str(ceiling/windows),
            'gram_dimension':sum(len(g['gram']) for g in groups.values()),
            'maximum_psd_dimension':max(len(g['gram']) for g in groups.values()),'sectors':sectors,
            'scope':'All-state joint overlapping bound on Phalf + ratio*Pcharged. Unequal source norms and weights are retained by the exact diagonal congruence. This bound retains cross-family overlaps.'}


def charged_projector_bound(source, windows, ceiling):
    ceiling = _exact(ceiling, 10**9)
    if type(windows) is not int or windows not in (2, 3, 4) or not 1 <= ceiling <= windows:
        raise ValueError('Invalid charged projector ceiling or window count')
    norm, grams = charged_overlap_grams(source, windows)
    sectors = []
    for key, gram in grams.items():
        shifted = [[(ceiling*norm if i == j else F(0))-x for j, x in enumerate(row)]
                   for i, row in enumerate(gram)]
        sectors.append({'sector': list(key), 'psd': _psd(shifted)})
    return {'accepted': True, 'windows': windows, 'support_sites': windows+5,
            'ceiling': str(ceiling), 'average_fidelity_ceiling': str(ceiling/windows),
            'source_count': 4, 'vector_norm': str(norm),
            'gram_dimension': sum(len(g) for g in grams.values()),
            'maximum_psd_dimension': max(len(g) for g in grams.values()), 'sectors': sectors,
            'scope': 'All-state sum of charged local projectors on consecutive six-site windows. Four orthogonal sources are derived by signed spin exchange and particle-hole transforms. Odd source embedding phases only conjugate the Gram by diagonal signs.'}
