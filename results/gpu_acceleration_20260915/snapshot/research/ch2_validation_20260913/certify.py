"""Exact spin-pure upper bounds and spin-penalized full-sector lower bounds."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import argparse
import json
import sys
import time

from experiments.marginal_symbolic import mono, add, scale, product, decode
from research.molecular_collective_20260913.core import digest
from research.certificate_scaling.streaming_reference_upper import compile_term, upper
from research.compact_response_20260913.closure import check_factor

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/ch2_validation_20260913'


def spin_square(modes):
    plus = add(*(mono(((1,2*i),(0,2*i+1))) for i in range(modes//2)))
    minus = add(*(mono(((1,2*i+1),(0,2*i))) for i in range(modes//2)))
    z = add(*(mono(((1,i),(0,i)), F(1 if i%2 == 0 else -1, 2)) for i in range(modes)))
    return add(product(z,z),scale(add(product(plus,minus),product(minus,plus)),F(1,2)))


def matrix(poly, states):
    den = lcm(*(c.denominator for c in poly.values())); where = {s:i for i,s in enumerate(states)}
    terms = [t for w,c in poly.items() if (t := compile_term(w,int(c*den))) is not None]
    M = [[0]*len(states) for _ in states]
    for j,s in enumerate(states):
        for required,occupied,flip,parity,c in terms:
            if s & required == occupied:
                dest = s^flip
                if dest not in where:
                    raise ValueError('Operator leaves the complete spin-projection sector')
                M[where[dest]][j] += c*(-1 if (s&parity).bit_count()%2 else 1)
    if any(M[i][j] != M[j][i] for i in range(len(states)) for j in range(i)):
        raise ValueError('Non-Hermitian exact sector matrix')
    return M,den,len(states)*len(terms)


def prepare(data, spin):
    start = time.monotonic()
    if spin not in (0,1) or type(spin) is not int:
        raise ValueError('Only the declared singlet and triplet are certified')
    if (data['modes'],data['particles'],data['physical_electrons'],data['total_spatial_orbitals'],data['frozen_core_orbitals']) != (12,6,8,7,1):
        raise ValueError('Neutral CH2/core/active electron accounting mismatch')
    alpha = list(combinations(range(6),3+spin)); beta = list(combinations(range(6),3-spin))
    states = sorted(sum(1<<(2*i) for i in a)+sum(1<<(2*i+1) for i in b) for a in alpha for b in beta)
    H,den,checks = matrix(decode(data['hamiltonian'],12,4),states)
    S,sd,schecks = matrix(spin_square(12),states)
    if any(x%sd for row in S for x in row):
        raise ValueError('This even-electron S2 matrix was not integral')
    S = [[x//sd for x in row] for row in S]; sr = [[(j,x) for j,x in enumerate(row) if x] for row in S]
    for i in range(len(states)):
        for j in range(len(states)):
            if sum(H[i][k]*v for k,v in sr[j]) != sum(v*H[k][j] for k,v in sr[i]):
                raise ValueError('Rational Hamiltonian does not commute with total spin')
    return states,H,S,den,{'determinants_enumerated':len(states),'Hamiltonian_word_state_checks':checks,
        'spin_word_state_checks':schecks,'H_and_S2_matrix_entries':2*len(states)**2,
        'preparation_seconds':time.monotonic()-start}


def spin_action(S,v):
    return [sum(x*a for x,a in zip(row,v) if x) for row in S]


def project_pure(S,v,spin):
    # N=6 in six active spatial orbitals: S runs up to 3. M_S=spin
    # excludes smaller spins. Integer polynomials suffice; normalization cancels.
    for other in range(spin+1,4):
        Sv=spin_action(S,v); v=[other*(other+1)*x-y for x,y in zip(v,Sv)]
    if not any(v):
        raise ValueError('Spin projection annihilated the numerical trial')
    common=gcd(*v); return [x//common for x in v]


def penalized(H,S,den,spin,penalty):
    return [[h+penalty*den*(S[i][j]-(spin*(spin+1) if i==j else 0)) for j,h in enumerate(row)] for i,row in enumerate(H)]


def propose(data,spin):
    import numpy as np
    start=time.monotonic(); states,H,S,den,cost=prepare(data,spin); penalty=20
    T=penalized(H,S,den,spin,penalty); tf=np.array(T,dtype=float)/float(den)
    e,V=np.linalg.eigh(tf); L=F(round(float(e[0])*10**10),10**10)-F(1,100000)
    v=project_pure(S,[int(round(float(x)*10**10)) for x in V[:,0]],spin)
    if spin_action(S,v) != [spin*(spin+1)*x for x in v]:
        raise ValueError('Exact trial does not have the declared total spin')
    witness={'states':[s for s,a in zip(states,v) if a],'amplitudes':[a for a in v if a]}
    U,upper_cost=upper(data,witness)
    shifted_den=lcm(den,L.denominator); multiplier=shifted_den//den
    K=[[x*multiplier-(int(L*shifted_den) if i==j else 0) for j,x in enumerate(row)] for i,row in enumerate(T)]
    margin=F(1,1000000); fd=10**12
    C=np.linalg.cholesky(tf-(float(L)+float(margin))*np.eye(len(T)))
    factor=[[int(round(float(C[i,j])*fd)) for j in range(i+1)] for i in range(len(T))]
    exact_margin=check_factor(K,shifted_den,F(0),factor,fd)
    cert={'kind':'ch2_total_spin_energy_interval_v1','fixture_sha256':digest(data),'total_spin':spin,
        'spin_projection':spin,'penalty_Ha':penalty,'lower_Ha':str(L),'upper_Ha':str(U),
        'independent_upper':witness,'factor_denominator':fd,'lower_factor':factor,'exact_factor_margin_Ha':str(exact_margin)}
    return cert,cost|{'upper_replay':upper_cost,'construction_seconds':time.monotonic()-start,
        'numerical_penalized_eigenvalue_Ha':float(e[0]),'width_Ha':str(U-L),'factor_entries':sum(map(len,factor))}


def check(data,cert):
    start=time.monotonic()
    if cert.get('kind')!='ch2_total_spin_energy_interval_v1' or cert['fixture_sha256']!=digest(data):
        raise ValueError('CH2 certificate binding failed')
    spin=cert['total_spin']
    if cert['spin_projection']!=spin:
        raise ValueError('Wrong spin projection for the declared proof')
    if type(cert['factor_denominator']) is not int or not 0<cert['factor_denominator']<=10**14:
        raise ValueError('Bounded positive exact factor denominator required')
    states,H,S,den,cost=prepare(data,spin); penalty=cert['penalty_Ha']
    if type(penalty) is not int or not 0<penalty<=100:
        raise ValueError('Invalid spin penalty')
    L=F(cert['lower_Ha']); T=penalized(H,S,den,spin,penalty)
    kd=lcm(den,L.denominator); mult=kd//den
    K=[[x*mult-(int(L*kd) if i==j else 0) for j,x in enumerate(row)] for i,row in enumerate(T)]
    if str(check_factor(K,kd,F(0),cert['lower_factor'],cert['factor_denominator']))!=cert['exact_factor_margin_Ha']:
        raise ValueError('Spin lower factor margin failed')
    witness=cert['independent_upper']; U,uc=upper(data,witness)
    if str(U)!=cert['upper_Ha'] or U<L:
        raise ValueError('CH2 upper did not reproduce')
    mapping=dict(zip(witness['states'],witness['amplitudes']))
    if any(s not in states for s in mapping):
        raise ValueError('Upper witness has wrong spin projection')
    v=[mapping.get(s,0) for s in states]
    if spin_action(S,v)!=[spin*(spin+1)*x for x in v]:
        raise ValueError('Upper witness is not exactly spin pure')
    return {'total_spin':spin,'spin_projection':spin,'lower_Ha':str(L),'upper_Ha':str(U),
        'width_Ha':str(U-L),'width_mHa':float(1000*(U-L)),'exact_spin_purity':True,
        'exact_H_commutes_S2':True,'cost':cost,'upper_replay':uc,'replay_seconds':time.monotonic()-start}


def run(replay=False):
    start=time.monotonic(); data=json.loads((OUT/'fixture.json').read_text()); rows=[]
    for spin in (0,1):
        path=OUT/f'spin_{spin}_certificate.json'; discovery=None
        if replay:
            cert=json.loads(path.read_text())
        else:
            cert,discovery=propose(data,spin);path.write_text(json.dumps(cert,separators=(',',':'))+'\n')
        receipt=check(data,cert);rows.append({'spin':spin,'discovery':discovery,'receipt':receipt,'certificate_bytes':path.stat().st_size})
    singlet,triplet=[x['receipt'] for x in rows]
    low=F(triplet['lower_Ha'])-F(singlet['upper_Ha']); high=F(triplet['upper_Ha'])-F(singlet['lower_Ha'])
    result={'states':rows,'gap_convention':'E_T-E_S','gap_lower_Ha':str(low),'gap_upper_Ha':str(high),
        'gap_width_mHa':float(1000*(high-low)),'gap_lower_approx_kcal_mol':float(low)*627.509474,
        'gap_upper_approx_kcal_mol':float(high)*627.509474,'wall_seconds':time.monotonic()-start,
        'numerical_packages_loaded':[x for x in ('numpy','scipy','pyscf','cvxpy') if x in sys.modules],
        'scope':'Exact gap for the declared frozen-core fixed-geometry rational Hamiltonian. No physical-model error or zero-point correction is certified.'}
    if replay and result['numerical_packages_loaded']:
        raise AssertionError('Numerical import during CH2 exact replay')
    (OUT/('exact_replay.json' if replay else 'certificate_discovery.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='states'}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--replay',action='store_true');run(parser.parse_args().replay)
