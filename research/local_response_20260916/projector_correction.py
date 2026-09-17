"""Exact correction for incompatible neighboring local ground subspaces.

Accepting functions use rational arithmetic and local blocks only. Numerical
proposal imports are confined to construct(). This is not a general closure of
recursive response bounds.
"""
from fractions import Fraction as F
import json
from math import floor, sqrt
from pathlib import Path
from research.local_response_20260916.local_exact import (
    local_blocks, require_psd, spin_counts, verify, verify_dual)


def norm2(v):
    return sum(x*x for x in v)


def rung_reduction(v, keep):
    if keep == 0:
        return [[sum(v[a+16*c]*v[b+16*c] for c in range(16))
                 for b in range(16)] for a in range(16)]
    return [[sum(v[c+16*a]*v[c+16*b] for c in range(16))
             for b in range(16)] for a in range(16)]


def grouped_overlap(vectors):
    left = rung_reduction(vectors[0], 1)
    right = rung_reduction(vectors[2], 0)
    middle = [(i&15, i>>4, x) for i,x in enumerate(vectors[1]) if x]
    value = sum(x*y*left[a][c]*right[b][d]
                for a,b,x in middle for c,d,y in middle)
    p = F(value, norm2(vectors[0])*norm2(vectors[1])*norm2(vectors[2]))
    if not 0 <= p <= 1:
        raise AssertionError('Projector overlap probability outside [0,1]')
    return p


def profile_information_ceiling(vectors, gaps):
    """Variational upper for the coarsened profile model, NOT physical H.

    |a>_A tensor u1_BC tensor |d>_D is in the middle projector's range.
    Restrict a,d to give the actual global Nup=Ndown=4. This is a finite
    local-product witness that limits what the three profile inequalities
    alone can imply about the correction above the baseline.
    """
    n0,n1,n2 = map(norm2,vectors)
    u0,u1,u2 = vectors
    center = {spin_counts(i,4) for i,x in enumerate(u1) if x}
    if len(center) != 1:
        raise ValueError('Definite middle charge required')
    cu,cd = next(iter(center))
    left = [F(sum(sum(u0[a+16*b]*u1[b+16*c] for b in range(16))**2
                  for c in range(16)),n0*n1) for a in range(16)]
    right = [F(sum(sum(u1[b+16*c]*u2[c+16*d] for c in range(16))**2
                   for b in range(16)),n1*n2) for d in range(16)]
    choices = []
    for a in range(16):
        au,ad = spin_counts(a,2)
        for d in range(16):
            du,dd = spin_counts(d,2)
            if au+du+cu == 4 and ad+dd+cd == 4:
                value = gaps[0]*(1-left[a])+gaps[2]*(1-right[d])
                choices.append((value,a,d))
    if not choices:
        return None
    value,a,d = min(choices)
    return {'profile_model_upper_over_t':str(value),
            'external_local_labels':[a,d],
            'physical_H_upper_claimed':False,
            'meaning':'No universal correction larger than this follows from these three fixed spectral-profile inequalities alone; the full local spectra contain additional information.'}


def profile_check(certificate, correction):
    blocks = local_blocks(certificate['messages'])
    shifts = list(map(F, certificate['lower_shifts']))
    entries = correction['local_profiles']
    if not isinstance(entries, list) or len(entries) != 3:
        raise ValueError('Exactly three local profiles required')
    vectors=[];gaps=[]
    for bi, entry in enumerate(entries):
        v = entry['integer_vector']
        if (not isinstance(v,list) or len(v)!=256 or
                any(type(x) is not int for x in v) or norm2(v)==0):
            raise ValueError('Nonzero 256-entry integer local vector required')
        support = {spin_counts(i,4) for i,x in enumerate(v) if x}
        if len(support)!=1:
            raise ValueError('Local vector must have a definite spin population')
        if type(entry['gap_over_t']) not in (int,str,F):
            raise ValueError('Exact local gap required')
        gap=F(entry['gap_over_t'])
        if gap<=0:
            raise ValueError('Positive local profile gap required')
        n=norm2(v)
        for labels,ham in blocks[bi].values():
            matrix=[[ham[i][j]+gap*F(v[x]*v[y],n)
                     -(shifts[bi]+gap if i==j else 0)
                     for j,y in enumerate(labels)] for i,x in enumerate(labels)]
            require_psd(matrix)
        vectors.append(v);gaps.append(gap)
    return vectors,gaps


def verify_correction(certificate, correction):
    primal=verify(certificate)
    if correction.get('kind')!='grouped_local_projector_correction_v1':
        raise ValueError('Wrong correction kind')
    vectors,gaps=profile_check(certificate,correction)
    p=grouped_overlap(vectors)
    gamma=F(correction['correction_over_t'])
    g=min(gaps[0],gaps[2]);h=gaps[1]
    if not 0<=gamma<=min(g,h) or (g-gamma)*(h-gamma)<g*h*p:
        raise ValueError('Unproved grouped-projector correction')
    lower=F(primal['energy_lower_over_t'])+gamma
    result={'status':'accepted_exact_nonlocal_projector_correction',
            'energy_lower_over_t':str(lower),
            'base_lower_over_t':primal['energy_lower_over_t'],
            'correction_over_t':str(gamma),'overlap_squared':str(p),
            'local_profile_gaps_over_t':list(map(str,gaps)),
            'max_PSD_dimension':36,'global_determinants_enumerated':0,
            'scope':'Fully coupled open2x4 Hubbard; local spectral profiles plus one grouped-projector inequality; no general recursion/scaling theorem'}
    result['coarsened_profile_information_limit'] = profile_information_ceiling(vectors,gaps)
    if 'dual_marginals' in certificate:
        dual=verify_dual(certificate)
        ceiling=F(dual['family_lower_ceiling_over_t'])
        result['old_family_lower_ceiling_over_t']=str(ceiling)
        result['physical_energy_minus_old_family_ceiling_at_least_t']=str(lower-ceiling)
    return result


def construct(certificate_path, out):
    import numpy as np
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    certificate=json.loads(certificate_path.read_text())
    verify(certificate)
    blocks=local_blocks(certificate['messages'])
    shifts=list(map(F,certificate['lower_shifts']))
    vectors=[];gaps=[];diagnostics=[]
    for bi,local in enumerate(blocks):
        spectra=[]
        for key,(labels,ham) in local.items():
            vals,vecs=np.linalg.eigh(np.array(ham,dtype=float)-float(shifts[bi])*np.eye(len(labels)))
            for j,e in enumerate(vals):spectra.append((float(e),key,vecs[:,j]))
        spectra.sort(key=lambda x:x[0])
        value,key,vec=spectra[0]
        integer=[0]*256
        for index,x in zip(local[key][0],vec):integer[index]=int(round(x*10**10))
        vectors.append(integer)
        gaps.append(F(floor((spectra[1][0]-1e-7)*10**10),10**10))
        diagnostics.append({'numerical_lowest':value,'numerical_second':spectra[1][0],
                            'lowest_charge':list(key)})
    correction={'kind':'grouped_local_projector_correction_v1'}
    failures=[]
    for factor in (F(1),F(999,1000),F(99,100),F(9,10),F(1,2)):
        correction['local_profiles']=[{'integer_vector':v,'gap_over_t':str(g*factor)}
                                      for v,g in zip(vectors,gaps)]
        try:
            verified_vectors,verified_gaps=profile_check(certificate,correction)
            break
        except ValueError as error:failures.append({'factor':str(factor),'reason':str(error)})
    else:raise RuntimeError('No positive local spectral profile accepted')
    p=grouped_overlap(verified_vectors)
    g=float(min(verified_gaps[0],verified_gaps[2]));h=float(verified_gaps[1])
    proposed=(g+h-sqrt((g-h)**2+4*g*h*float(p)))/2
    correction['correction_over_t']=str(F(max(0,floor((proposed-1e-9)*10**10)),10**10))
    receipt=verify_correction(certificate,correction)
    receipt['local_numerical_proposals']=diagnostics
    receipt['profile_repair_failures']=failures
    (out/'correction.json').write_text(json.dumps(correction,separators=(',',':'))+'\n')
    (out/'construction.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    construct(a.certificate,a.out)
