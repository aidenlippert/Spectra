"""Exact fixed-width norm and energy transfer for disjoint Hubbard cut filters.

Local tensors are derived from an explicit bounded physical block state. No
supplied marginal or transfer matrix is accepted as a physical certificate.
"""
from fractions import Fraction as F
from experiments.marginal_local_hubbard_block import _exact,_terms,_actions,replay as replay_local
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_boundary_unitary import _physical_source


def _move(operator,state):
    out={}
    for s,a in state.items():
        for t,b in operator[s].items():out[t]=out.get(t,0)+a*b
    return {s:a for s,a in out.items() if a}


def _compose(a,b):return {s:_move(a,v) for s,v in b.items()}


def _add(a,b):
    return {s:{t:a[s].get(t,0)+b[s].get(t,0) for t in a[s].keys()|b[s].keys()
               if a[s].get(t,0)+b[s].get(t,0)} for s in a}


def _act(terms,state):
    out={}
    for s,a in state.items():
        for word,b in terms:
            image=apply_word(word,s)
            if image:
                t,sign=image;out[t]=out.get(t,0)+sign*a*b
    return {s:a for s,a in out.items() if a}


def _partial(v,w,selected,sites):
    mask=sum(3<<(2*i) for i in selected)
    def groups(state):
        result={}
        for s,a in state.items():
            index=sum(((s>>(2*i))&3)<<(2*j) for j,i in enumerate(selected))
            result.setdefault(s&~mask,{})[index]=a
        return result
    a,b=groups(v),groups(w);out={}
    for middle,rows in a.items():
        for i,x in rows.items():
            for j,y in b.get(middle,{}).items():out[i,j]=out.get((i,j),0)+x*y
    return {ij:x for ij,x in out.items() if x}


def _edge_tensor(v,w,sites):
    rho=_partial(v,w,(0,sites-1),sites)
    return [[rho.get((l+4*r,lp+4*rp),0) for rp in range(4) for r in range(4)]
            for lp in range(4) for l in range(4)]


def _cut_operators(U=4,t=1,V=0):
    h=_actions(2,0,1)
    identity={s:{s:1} for s in range(16)}
    powers=[identity,{s:{t:-a for t,a in v.items()} for s,v in h.items()},_compose(h,h)]
    embedded=[{s:{(s&195)+(t<<2):a for t,a in op[(s>>2)&15].items()}
               for s in range(256)} for op in powers]
    terms=_terms(4,U,t,[0,U,U,0],[t,t,t],[V,V,V])
    patch={s:_act(terms,{s:1}) for s in range(256)}
    result=[]
    for m in range(3):
        for n in range(m,3):
            norm=_compose(powers[m],powers[n])
            energy=_compose(embedded[m],_compose(patch,embedded[n]))
            if m!=n:
                norm=_add(norm,_compose(powers[n],powers[m]))
                energy=_add(energy,_compose(embedded[n],_compose(patch,embedded[m])))
            if any(energy[t].get(s,0)!=a for s,v in energy.items() for t,a in v.items()):
                raise ValueError('Dressed patch must be exactly Hermitian')
            B=[[norm[r+4*l].get(rp+4*lp,0) for lp in range(4) for l in range(4)]
               for rp in range(4) for r in range(4)]
            result.append((m,n,B,energy))
    return result


def compile_state(state,sites,U=4,t=1,V=0):
    """Compile an actual fixed-spin integer state; support bounded by H8."""
    if type(sites) is not int or sites not in (4,8):raise ValueError('Block size must be4 or8')
    U,t,V=_exact(U),_exact(t),_exact(V)
    if U<0 or t<0:raise ValueError('Nonnegative onsite repulsion and hopping required')
    if type(state) is not dict or not 1<=len(state)<=4900:raise ValueError('Bounded physical state required')
    if any(type(s) is not int or not 0<=s<4**sites or s.bit_count()!=sites
           or sum((s>>(2*i))&1 for i in range(sites))!=sites//2
           or type(a) is not int or abs(a)>10**512 for s,a in state.items()):
        raise ValueError('Exact bounded fixed-spin physical amplitudes required')
    state={s:a for s,a in state.items() if a}
    norm=sum(a*a for a in state.values())
    if not norm:raise ValueError('Nonzero physical state required')
    G=_edge_tensor(state,state,sites)
    onsite=[0]+[U]*(sites-2)+[0]
    hopping=[0]+[t]*(sites-3)+[0]
    mid=_act(_terms(sites,U,t,onsite,hopping,[0]+[V]*(sites-3)+[0]),state)
    left=_act(_terms(sites,U,t,[U]+[0]*(sites-1),[t]+[0]*(sites-2),[V]+[0]*(sites-2)),state)
    right=_act(_terms(sites,U,t,[0]*(sites-1)+[U],[0]*(sites-2)+[t],[0]*(sites-2)+[V]),state)
    J,L,R=(_edge_tensor(v,state,sites) for v in (mid,left,right))
    rhoA=_partial(state,state,(0,sites-2,sites-1),sites)
    rhoB=_partial(state,state,(0,1,sites-1),sites)
    # Group each three-site partial trace by its two inner sites. Its
    # remaining paired single-site indices are the transfer's16 coordinates.
    A={};B={}
    for (s,bra),x in rhoA.items():A.setdefault((s>>2,bra>>2),[]).append((s%4+4*(bra%4),x))
    for (s,bra),x in rhoB.items():B.setdefault((s%16,bra%16),[]).append((s//16+4*(bra//16),x))
    kernels=[]
    for m,n,norm_operator,energy in _cut_operators(U,t,V):
        W=[[0]*16 for _ in range(16)]
        for s,image in energy.items():
            for bra,z in image.items():
                for i,x in A.get((s%16,bra%16),[]):
                    for j,y in B.get((s//16,bra//16),[]):W[i][j]+=x*y*z
        kernels.append((m,n,norm_operator,W))
    return {'sites':sites,'block_norm':norm,'G':G,'J':J,'L':L,'R':R,'kernels':kernels,
            'physical_support':len(state),'target':{'U':str(U),'t':str(t),'V':str(V)}}


def compile_block(c,h,U=4,t=1,V=0):
    oracle,state,checked=_physical_source(c,h)
    full={s:a*phase for r,a in state.items() for s,phase in oracle.orbit(r)[3].items()}
    result=compile_state(full,8,U,t,V);result['source_upper_replay']=checked
    return result


def _vm(v,M):return [sum(a*M[i][j] for i,a in enumerate(v) if a) for j in range(16)]
def _mm(a,b):return [_vm(row,b) for row in a]
def _sum(*vectors):return [sum(items) for items in zip(*vectors)]


def matrices(compiled,a,b):
    a,b=_exact(a),_exact(b);coefficients=[1,a,b]
    B=[[0]*16 for _ in range(16)];W=[[0]*16 for _ in range(16)]
    for m,n,norm,energy in compiled['kernels']:
        coefficient=coefficients[m]*coefficients[n]
        if coefficient:
            for i in range(16):
                for j in range(16):
                    B[i][j]+=coefficient*norm[i][j]
                    W[i][j]+=coefficient*energy[i][j]
    return B,W


def contract(compiled,a,b,blocks):
    """Exact finite recurrence; integer bit growth is bounded by blocks<=64."""
    if type(blocks) is not int or not 1<=blocks<=64:raise ValueError('Exact chain requires1..64 blocks')
    B,W=matrices(compiled,a,b)
    G,J,L,R=(compiled[k] for k in ('G','J','L','R'))
    e=[int(i%4==i//4) for i in range(16)]
    trace=lambda v:sum(x*y for x,y in zip(v,e))
    current=_vm(e,G)
    energy=_sum(_vm(e,J),_vm(e,L))
    if blocks==1:
        numerator=trace(_sum(energy,_vm(e,R)));norm=trace(current)
    else:
        T,BJ,BW=_mm(B,G),_mm(B,J),_mm(B,W)
        previous=current
        energy=_sum(_vm(energy,T),_vm(current,BJ),_vm(e,W))
        current=_vm(current,T)
        for _ in range(3,blocks+1):
            energy=_sum(_vm(energy,T),_vm(current,BJ),_vm(previous,BW))
            previous,current=current,_vm(current,T)
        numerator=trace(_sum(energy,_vm(_vm(previous,B),R)));norm=trace(current)
    if norm<=0:raise ValueError('Filtered physical state has zero norm')
    return {'accepted':True,'sites':blocks*compiled['sites'],'blocks':blocks,'a':str(F(a)),'b':str(F(b)),
            'norm':str(norm),'energy_numerator':str(numerator),'energy':str(F(numerator,norm)),
            'upper_per_site':str(F(numerator,norm)/(blocks*compiled['sites'])),'target':compiled['target'],
            'scope':'Exact physical norm and specified nearest-neighbor U,t,V energy from fixed-width boundary tensors. Exact recurrence capped at64 blocks; no large-chain precision or complexity guarantee.'}


def replay(c,h,a,b,blocks,U=4,t=1,V=0):return contract(compile_block(c,h,U,t,V),a,b,blocks)


def _dyadic(matrix,bits):
    unit=1<<bits
    lo=[];hi=[]
    for row in matrix:
        values=[F(x)*unit for x in row]
        lo.append([x.numerator//x.denominator for x in values])
        hi.append([-((-x.numerator)//x.denominator) for x in values])
    return lo,hi,0


def _interval_product(a,b,bits):
    """Outward integer-grid multiplication with one common binary scale."""
    al,ah,ae=a;bl,bh,be=b
    rows,inner,cols=len(al),len(bl),len(bl[0]);unit=1<<bits
    low=[[0]*cols for _ in range(rows)];high=[[0]*cols for _ in range(rows)]
    for i in range(rows):
        for k in range(inner):
            x,y=al[i][k],ah[i][k]
            if not x and not y:continue
            for j in range(cols):
                u,v=bl[k][j],bh[k][j]
                if not u and not v:continue
                products=(x*u,x*v,y*u,y*v)
                low[i][j]+=min(products);high[i][j]+=max(products)
    low=[[x//unit for x in row] for row in low]
    high=[[-((-x)//unit) for x in row] for row in high]
    largest=max(abs(x) for matrix in (low,high) for row in matrix for x in row)
    shift=largest.bit_length()-bits if largest else 0
    if shift>0:
        factor=1<<shift
        low=[[x//factor for x in row] for row in low]
        high=[[-((-x)//factor) for x in row] for row in high]
    elif shift<0:
        low=[[x<<(-shift) for x in row] for row in low]
        high=[[x<<(-shift) for x in row] for row in high]
    return low,high,ae+be+shift


def enclose(compiled,a,b,blocks,bits=160):
    """Enclose the physical energy using outward-rounded transfer powers.

    A positive common scale is carried symbolically and cancels in E/norm.
    A failed positivity enclosure is a refusal, never an accepted energy.
    """
    if type(blocks) is not int or not 1<=blocks<=125000000:
        raise ValueError('Enclosed chain requires1..125000000 blocks')
    if type(bits) is not int or not 64<=bits<=512:raise ValueError('Precision must be64..512 bits')
    B,W=matrices(compiled,a,b)
    if blocks==1:
        r=contract(compiled,a,b,1);energy=F(r['energy'])
        return {'accepted':True,'sites':r['sites'],'blocks':1,'a':r['a'],'b':r['b'],
                'energy_lower':str(energy),'energy_upper':str(energy),
                'upper_per_site':str(energy/r['sites']),'width_per_site':'0',
                'precision_bits':bits,'matrix_squarings':0,'transfer_dimension':48,'target':compiled['target']}
    norm=compiled['block_norm']
    G,J,L,R=([[F(x,norm) for x in row] for row in compiled[k]] for k in ('G','J','L','R'))
    W=[[F(x,norm*norm) for x in row] for row in W]
    T,BJ,BW=_mm(B,G),_mm(B,J),_mm(B,W)
    e=[int(i%4==i//4) for i in range(16)]
    first=_vm(e,G)
    energy=_sum(_vm(_sum(_vm(e,J),_vm(e,L)),T),_vm(first,BJ),_vm(e,W))
    seed=_vm(first,T)+first+energy
    matrix=[[0]*48 for _ in range(48)]
    for i in range(16):
        matrix[i][16+i]=1
        for j in range(16):
            matrix[i][j]=T[i][j]
            matrix[i][32+j]=BJ[i][j]
            matrix[16+i][32+j]=BW[i][j]
            matrix[32+i][32+j]=T[i][j]
    row=_dyadic([seed],bits);power=_dyadic(matrix,bits)
    exponent=blocks-2;squarings=0;products=0
    while exponent:
        if exponent&1:row=_interval_product(row,power,bits);products+=1
        exponent>>=1
        if exponent:power=_interval_product(power,power,bits);squarings+=1
    BR=_mm(B,R)
    tail=[sum(x*y for x,y in zip(v,e)) for v in BR]
    closing=[[e[i] if i<16 else 0,
              tail[i-16] if 16<=i<32 else (e[i-32] if i>=32 else 0)] for i in range(48)]
    result=_interval_product(row,_dyadic(closing,bits),bits)
    lo,hi,scale=result
    nl,nh=lo[0][0],hi[0][0];el,eh=lo[0][1],hi[0][1]
    if nl<=0:raise ValueError('Rounded contraction could not certify strictly positive norm')
    quotients=[F(x,y) for x in (el,eh) for y in (nl,nh)]
    lower,upper=min(quotients),max(quotients);sites=blocks*compiled['sites']
    return {'accepted':True,'sites':sites,'blocks':blocks,'a':str(F(a)),'b':str(F(b)),
            'energy_lower':str(lower),'energy_upper':str(upper),'upper_per_site':str(upper/sites),
            'width_per_site':str((upper-lower)/sites),'precision_bits':bits,
            'matrix_squarings':squarings,'row_power_products':products,'transfer_dimension':48,
            'norm_scaled_interval':[str(nl),str(nh)],'numerator_scaled_interval':[str(el),str(eh)],
            'shared_binary_exponent':scale-bits,'target':compiled['target'],
            'scope':'Outward-rounded enclosure of a specified physical filtered-state energy, not a ground-energy lower bound. The upper endpoint is variational. Fixed48-dimensional transfer powers retain a symbolic common scale; positive norm is certified explicitly.'}


def replay_target(c,h,a,b,blocks,local_certificate,U=4,t=1,V=0,bits=160):
    """Pair a freshly compiled trial upper with a matching all-Fock window."""
    U,t,V=_exact(U),_exact(t),_exact(V)
    if type(local_certificate) is not dict:raise ValueError('Physical local lower certificate required')
    length=local_certificate.get('sites')
    if type(length) is not int or length not in (2,4,6):raise ValueError('Supported local window required')
    if (_exact(local_certificate.get('U'))*length!=(length-1)*U
        or _exact(local_certificate.get('t'))!=t or _exact(local_certificate.get('V',0))!=V):
        raise ValueError('Translated window does not reproduce the specified target')
    compiled=compile_block(c,h,U,t,V)
    upper=enclose(compiled,a,b,blocks,bits)
    if upper['sites']<length:raise ValueError('Chain too short for the periodic window cover')
    local=replay_local(local_certificate)
    boundary_norm=2*t+abs(V)
    lo=F(local['lower'])/(length-1)-boundary_norm/upper['sites']
    hi=F(upper['upper_per_site'])
    if lo>hi:raise ValueError('Ground-energy interval is inconsistent')
    return {'accepted':True,'target':compiled['target'],'sites':upper['sites'],
            'physical_source_replay':compiled['source_upper_replay'],'upper_replay':upper,'lower_replay':local,
            'lower_per_site':str(lo),'upper_per_site':str(hi),'width_per_site':str(hi-lo),
            'closing_bond_norm_bound':str(boundary_norm),
            'scope':'Half-filled open nearest-neighbor U,t,V ground-energy interval. Fresh physical trial-state upper and all-Fock local lower; translated profile sums reproduce the target. Closing-bond norm uses2t+abs(V). No ground-spin assumption, target gate optimum, or generic molecular transfer is claimed.'}
