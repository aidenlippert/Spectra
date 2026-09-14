"""Bounded exact overlapping-projector inequalities and Hubbard energy proofs.

Sparse isometry columns replace a full many-body projector matrix. All Gram
columns and symmetry sectors are reconstructed; supplied numerical spectra are
never trusted. Two through five consecutive four- or six-site windows are admitted.
"""
from fractions import Fraction as F
from math import lcm

from experiments.marginal_local_hubbard_block import _exact, _reflection, _sector, sector_matrices, _profiles, _density_profile
from experiments.marginal_polynomial_sos import integer_psd
from experiments.marginal_window_family_bound import _vector, local_certificate


def _psd(matrix,proof=None):
    scale = lcm(*(x.denominator for row in matrix for x in row))
    return (proof.check if proof is not None else integer_psd)([[int(x*scale) for x in row] for row in matrix])


def _rank_one_psd(base,weights,penalty,norm,proof=None):
    """Retain the rank-one denominator instead of powering it into minors."""
    scale=lcm(*(x.denominator for row in base for x in row))
    coefficient=F(penalty*scale,norm)
    r,s=coefficient.numerator,coefficient.denominator
    matrix=[[s*int(x*scale)+r*weights[i]*weights[j] for j,x in enumerate(row)]
            for i,row in enumerate(base)]
    return (proof.check if proof is not None else integer_psd)(matrix,initial_divisor=s)


def _projector_vector(vector,block_sites):
    if type(block_sites) is not int or block_sites not in (4,6):raise ValueError('Projector block must have4 or6 sites')
    if block_sites==4:return _vector(vector)
    if type(vector) is not dict or not 1<=len(vector)<=4096:raise ValueError('Bounded physical projector vector required')
    result={}
    for key,a in vector.items():
        if type(key) not in (int,str):raise ValueError('Integer determinant label required')
        s=int(key)
        if s in result or not 0<=s<4096 or s.bit_count()!=6 or type(a) is not int or abs(a)>10**12:
            raise ValueError('Bounded half-filled integer projector required')
        result[s]=a
    norm=sum(a*a for a in result.values())
    if not norm:raise ValueError('Nonzero projector vector required')
    return result,norm


def overlap_grams(vector,windows,block_sites=4):
    """Construct every Gram sector from physical contiguous embeddings."""
    if type(windows) is not int or windows not in (2,3,4,5):
        raise ValueError('Projector window count must be2,3,4,5')
    amplitudes,norm = _projector_vector(vector,block_sites)
    amplitudes = {s:a for s,a in amplitudes.items() if a}
    spins = {_sector(s,block_sites) for s in amplitudes}
    fixed_spin = len(spins) == 1
    groups = {}
    for offset in range(windows):
        for environment in range(4**(windows-1)):
            left = environment & ((1<<(2*offset))-1)
            right = environment >> (2*offset)
            column = {left | (s<<(2*offset)) | (right<<(2*(offset+block_sites))):a
                      for s,a in amplitudes.items()}
            # Every column has fixed total particle number; use the finer spin
            # split only when the actual input vector has that symmetry.
            keys = {_sector(s,windows+block_sites-1) if fixed_spin else (s.bit_count(),)
                    for s in column}
            if len(keys) != 1 or sum(a*a for a in column.values()) != norm:
                raise ValueError('Invalid embedded isometry column')
            groups.setdefault(next(iter(keys)),[]).append(column)
    grams={}
    for sector,columns in sorted(groups.items()):
        grams[sector]=[[sum(a*right.get(s,0) for s,a in left.items()) for right in columns] for left in columns]
    if sum(map(len,grams.values()))!=windows*4**(windows-1):raise ValueError('Incomplete projector Gram coverage')
    return norm,grams,fixed_spin


def projector_bound(vector, windows, ceiling,block_sites=4,proof=None):
    """Prove sum of m adjacent projectors <= ceiling I, using their Gram."""
    ceiling = _exact(ceiling,10**9)
    norm,grams,fixed_spin=overlap_grams(vector,windows,block_sites)
    if not 1<=ceiling<=windows:raise ValueError('Projector-sum ceiling must lie between1 and window count')
    receipts=[]
    for sector,gram in grams.items():
        matrix=[[(ceiling*norm if i==j else F(0))-x for j,x in enumerate(row)] for i,row in enumerate(gram)]
        receipts.append({'sector':list(sector),'psd':_psd(matrix,proof)})
    dimension = sum(r['psd']['dimension'] for r in receipts)
    if dimension != windows*4**(windows-1):
        raise ValueError('Incomplete projector Gram coverage')
    return {'accepted':True,'windows':windows,'block_sites':block_sites,'support_sites':windows+block_sites-1,
            'ceiling':str(ceiling),'average_fidelity_ceiling':str(ceiling/windows),
            'vector_norm':str(norm),'gram_dimension':dimension,
            'maximum_psd_dimension':max(r['psd']['dimension'] for r in receipts),
            'spin_sector_split':fixed_spin,'sectors':receipts,
            'scope':'All-state inequality for consecutive contiguous even fixed-length projectors. Cyclic copies are defined by conjugation with the fermionic mode-translation unitary.'}


def replay(certificate,*,psd_witnesses=None):
    if type(certificate) is dict and certificate.get('kind') in ('hubbard_projector_extension_v2','hubbard_projector_extension_v4','hubbard_projector_extension_v5','hubbard_projector_extension_v6'):
        from experiments.marginal_congruence_psd import CongruenceProofs
        proof=CongruenceProofs(psd_witnesses) if psd_witnesses is not None else None
        result=_replay_extended(certificate,proof)
        if proof is not None:result['congruence_witnesses']=proof.finish()
        return result
    if psd_witnesses is not None:
        raise ValueError('Congruence witnesses require an extended certificate')
    if type(certificate) is not dict or certificate.get('kind') != 'hubbard_projector_extension_v1':
        raise ValueError('Unsupported projector extension certificate')
    sites = certificate.get('chain_sites')
    windows = certificate.get('windows')
    if type(sites) is not int or sites%2 or not 8 <= sites <= 10**9:
        raise ValueError('Even chain length between8 and10^9 required')
    vector = certificate.get('vector')
    amplitudes,norm = _vector(vector)
    amplitudes = {s:a for s,a in amplitudes.items() if a}
    if len({_sector(s,4) for s in amplitudes}) != 1:
        raise ValueError('Energy projector must preserve spin-number sectors')
    reflected = {}
    for s,a in amplitudes.items():
        r,phase = _reflection(s,4)
        reflected[r] = phase*a
    if not any(reflected == {s:parity*a for s,a in amplitudes.items()} for parity in (1,-1)):
        raise ValueError('Energy projector must preserve reflection sectors')
    penalty = _exact(certificate.get('penalty'),10**9)
    lower = _exact(certificate.get('penalized_lower'),10**9)
    if penalty < 0:
        raise ValueError('Projector penalty must be nonnegative')
    overlap = projector_bound(vector,windows,certificate.get('projector_sum_ceiling'))
    if sites <= overlap['support_sites']:
        raise ValueError('Chain must exceed projector inequality support')
    local = local_certificate(certificate.get('a'),certificate.get('b'),lower)
    data = sector_matrices(4,3,1,local['onsite_profile'],local['hopping_profile'])
    sectors = []
    for key,(_,matrix,columns) in data.items():
        weights = [sum(a*amplitudes.get(s,0) for s,a in c.items()) for c in columns]
        norms = [sum(a*a for a in c.values()) for c in columns]
        shifted = [[x+penalty*F(weights[i]*weights[j],norm)-(lower*norms[i] if i==j else 0)
                    for j,x in enumerate(row)] for i,row in enumerate(matrix)]
        sectors.append({'sector':list(key),'psd':_psd(shifted)})
    theta = F(overlap['average_fidelity_ceiling'])
    periodic = (lower-penalty*theta)/3
    opened = periodic-F(2,sites)
    return {'accepted':True,'chain_sites':sites,'U':'4','t':'1',
            'penalty':str(penalty),'penalized_lower':str(lower),
            'overlap':overlap,'local_sectors':sectors,
            'local_sum_dimensions':sum(r['psd']['dimension'] for r in sectors),
            'local_maximum_psd_dimension':max(r['psd']['dimension'] for r in sectors),
            'periodic_lower_density':str(periodic),'open_lower_density':str(opened),
            'open_lower_energy':str(opened*sites),
            'scope':'Exact lower bound for the half-filled uniform Hubbard chain at U=4,t=1. Local positivity holds on the full Fock space; the centered bulk shift vanishes at half filling. This is an extendibility constraint, not a general representability oracle.'}


def _telescoping_diagonal(source):
    """A reflection-odd five-site diagonal produces an even six-site difference."""
    if type(source) is not dict or not 1<=len(source)<=64:
        raise ValueError('One through64 five-site diagonal entries required')
    result={}
    for label,value in source.items():
        if type(label) not in (str,int):raise ValueError('Integer five-site determinant required')
        state=int(label)
        if state in result or not 0<=state<1024:raise ValueError('Distinct five-site determinants required')
        result[state]=_exact(value)
    result={s:v for s,v in result.items() if v}
    if {_reflection(s,5)[0]:v for s,v in result.items()}!={s:-v for s,v in result.items()}:
        raise ValueError('Five-site diagonal must be reflection-odd')
    return result


def _replay_extended(c,proof=None):
    sites=c.get('chain_sites');local=c.get('local_window');target=c.get('target')
    range_two=c.get('kind')=='hubbard_projector_extension_v6'
    expected_local='local_hubbard_range2_block_v1' if range_two else 'local_hubbard_block_v1'
    if type(sites) is not int or sites%2 or not 8<=sites<=10**9:raise ValueError('Bounded even chain required')
    if type(local) is not dict or local.get('kind')!=expected_local or type(target) is not dict:
        raise ValueError('Explicit local window and target required')
    L=local.get('sites')
    if type(L) is not int or L not in (4,6):raise ValueError('Projector local window must have4 or6 sites')
    U,t,V=(_exact(target.get(k)) for k in ('U','t','V'))
    from experiments.marginal_range_two_density import local_profile, diagonal_value
    W=F(0);range_profile=[]
    if range_two:
        W=_exact(target.get('W'));range_profile=local_profile(L,W,local.get('range_two_density_profile'))
    elif 'W' in target or 'range_two_density_profile' in local:
        raise ValueError('Range-two density requires v6 certificate')
    u,hop,v=(_exact(local.get(k,0)) for k in ('U','t','V'))
    if U<0 or t<0 or L*u!=(L-1)*U or hop!=t or v!=V:
        raise ValueError('Local window does not reproduce the target')
    onsite,hopping=_profiles(L,u,hop,local.get('onsite_profile'),local.get('hopping_profile'))
    _,density=_density_profile(L,V,local.get('density_profile'))
    vector,norm=_projector_vector(c.get('vector'),L);vector={s:a for s,a in vector.items() if a}
    if len({_sector(s,L) for s in vector})!=1:raise ValueError('Energy projector must preserve spin-number sectors')
    reflected={}
    for s,a in vector.items():
        r,phase=_reflection(s,L);reflected[r]=phase*a
    if not any(reflected=={s:p*a for s,a in vector.items()} for p in (-1,1)):
        raise ValueError('Energy projector must preserve reflection sectors')
    penalty=_exact(c.get('penalty'),10**9);lower=_exact(c.get('penalized_lower'),10**9)
    if penalty<0:raise ValueError('Nonnegative projector penalty required')
    joint = c.get('kind') in ('hubbard_projector_extension_v4','hubbard_projector_extension_v5','hubbard_projector_extension_v6')
    telescope={}
    if c.get('kind') in ('hubbard_projector_extension_v5','hubbard_projector_extension_v6'):
        telescope=_telescoping_diagonal(c.get('telescoping_diagonal'))
    elif 'telescoping_diagonal' in c:
        raise ValueError('Telescoping correction requires v5 certificate')
    beta=F(0);ratio=F(0);charged=[];charged_norm=1
    if joint:
        from experiments.marginal_charged_projectors import charged_vectors, joint_projector_bound
        spec=c.get('joint')
        if L!=6 or {_sector(s,L) for s in vector}!={(3,3)} or type(spec) is not dict:
            raise ValueError('Joint energy proof requires six-site half and charged sources')
        beta=_exact(spec.get('penalty'),10**9);ratio=_exact(spec.get('ratio'),10**9)
        if beta<0 or ratio<=0:raise ValueError('Nonnegative joint penalty and positive ratio required')
        charged,charged_norm=charged_vectors(spec.get('vector'))
        joint_overlap=joint_projector_bound(c.get('vector'),spec.get('vector'),spec.get('windows'),
                                            ratio,spec.get('projector_sum_ceiling'),proof=proof)
        if sites<=joint_overlap['support_sites']:
            raise ValueError('Chain must exceed joint projector inequality support')
    overlap=projector_bound(c.get('vector'),c.get('windows'),c.get('projector_sum_ceiling'),L,proof=proof)
    if sites<=overlap['support_sites']:raise ValueError('Chain must exceed projector inequality support')
    data=sector_matrices(L,u,hop,onsite,hopping,V,density);receipts=[]
    for key,(_,matrix,columns) in data.items():
        # Disjoint spin sectors permit at most one rank-one update per block.
        source,source_norm,coefficient=vector,norm,penalty+beta
        if joint:
            matches=[v for v in charged if _sector(next(iter(v)),L)==key[:2]]
            if len(matches)>1:raise ValueError('Charged sources must occupy disjoint spin sectors')
            if matches:source,source_norm,coefficient=matches[0],charged_norm,beta*ratio
        weights=[sum(a*source.get(s,0) for s,a in col.items()) for col in columns]
        norms=[sum(a*a for a in col.values()) for col in columns]
        corrections=[]
        for col in columns:
            values={telescope.get(s&1023,F(0))-telescope.get(s>>2,F(0)) + (diagonal_value(s,range_profile) if range_two else 0) for s in col}
            if len(values)!=1:raise ValueError('Telescoping term must preserve reflection blocks')
            corrections.append(values.pop())
        base=[[x+((corrections[i]-lower)*norms[i] if i==j else 0) for j,x in enumerate(row)] for i,row in enumerate(matrix)]
        receipts.append({'sector':list(key),'psd':_rank_one_psd(base,weights,coefficient,source_norm,proof)})
    theta=F(overlap['average_fidelity_ceiling'])
    joint_cost=beta*F(joint_overlap['average_ceiling']) if joint else F(0)
    periodic=(lower-penalty*theta-joint_cost)/(L-1);opened=periodic-(2*t+abs(V)+2*abs(W))/sites
    result={'accepted':True,'target':{'U':str(U),'t':str(t),'V':str(V)},'chain_sites':sites,
            'local_sites':L,'penalty':str(penalty),'penalized_lower':str(lower),'overlap':overlap,
            'local_sectors':receipts,'local_sum_dimensions':sum(r['psd']['dimension'] for r in receipts),
            'local_maximum_psd_dimension':max(r['psd']['dimension'] for r in receipts),
            'periodic_lower_density':str(periodic),'open_lower_density':str(opened),
            'open_lower_energy':str(opened*sites),
            'scope':'Exact half-filled open U,t,V lower bound from full-Fock local positivity with a rank-one projector penalty and a separately certified all-state overlapping-projector ceiling. This adds quantitative extendibility information beyond local profile or telescoping-boundary corrections.'}
    if joint:
        result.update(joint_overlap=joint_overlap,joint_penalty=str(beta),
                      half_local_penalty=str(penalty+beta),charged_local_penalty=str(beta*ratio),
                      scope='Exact half-filled open U,t,V lower bound. Full-Fock six-site positivity uses disjoint half/charged rank-one updates. Independently translated half and joint weighted overlapping-projector inequalities are combined on the periodic chain before removing the closing Hamiltonian bond. No general representability or accuracy-versus-cost guarantee is asserted.')
    if range_two:
        result['target']['W']=str(W)
        result.update(range_two_density_profile=list(map(str,range_profile)),closing_interaction_norm=str(2*t+abs(V)+2*abs(W)),
                      scope='Exact half-filled open U,t,V,W chain lower with next-nearest density W*q_i*q_(i+2). Fresh all-Fock local positivity and fixed overlapping-projector ceilings are combined. Local range-two profiles sum to five W for six-site windows. Opening removes one nearest bond and two range-two bonds; no general molecular representability claim.')
    if c.get('kind') in ('hubbard_projector_extension_v5','hubbard_projector_extension_v6'):
        result.update(telescoping_diagonal={str(s):str(v) for s,v in sorted(telescope.items())},
                      telescoping_nonzero_entries=len(telescope),
                      telescoping_scope='The six-site diagonal is Y(first five sites)-Y(last five sites). Every translated five-site term cancels exactly on the periodic chain, including cyclic copies defined by fermionic translation. Reflection-odd Y makes the local correction reflection-even; all-Fock positivity includes this correction before opening the Hamiltonian chain.')
    return result
