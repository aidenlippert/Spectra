"""Select spin-change-one paired cubic directions, retaining the earlier block."""
import argparse
from itertools import combinations
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy import linalg,sparse
from research.acceptance_channels_20260915.dense_t2 import tensors,contracted_entry,normal_word


def general_tensor(words,coefficients,modes):
    # The contraction helper uses separate index labels for its creator and
    # annihilator arrays. Map them back to physical modes after contracting.
    labelled=[tuple((c,2*k+(1-c)) for c,k in w) for w in words]
    return tensors(labelled,coefficients,modes)


def general_entry(left,right,modes):
    for word,value in contracted_entry(left,right,modes):
        yield tuple((c,k//2) for c,k in word),value


def moment_arrays(moments,modes):
    one=np.zeros((modes,modes));four=np.zeros((modes,)*4)
    known_one=np.zeros((modes,modes),dtype=bool);known_four=np.zeros((modes,)*4,dtype=bool)
    for a in range(modes):
        known_four[a,a,:,:]=True;known_four[:,:,a,a]=True
    for word,value in moments.items():
        if len(word)==2:
            a,b=[k for _,k in word];one[a,b]=one[b,a]=value;known_one[a,b]=known_one[b,a]=True
        elif len(word)==4:
            a,b,c,d=[k for _,k in word]
            for p,q,r,s in ((a,b,c,d),(c,d,a,b)):
                for u,v,sgn1 in ((p,q,1),(q,p,-1)):
                    for w,z,sgn2 in ((r,s,1),(s,r,-1)):
                        four[u,v,w,z]=sgn1*sgn2*value;known_four[u,v,w,z]=True
    return one,four,known_one,known_four


def operator_moment_matrix(words,moments,modes):
    one,G,known_one,known=moment_arrays(moments,modes)
    triples=np.array([[k for _,k in w] for w in words]);p,q,r=triples.T
    M=np.zeros((len(words),len(words)))
    for i,(a,b,c) in enumerate(triples):
        def add(mask,indices,sign=1):
            ids=np.flatnonzero(mask)
            if not len(ids):return
            indices=tuple(v[ids] if isinstance(v,np.ndarray) else v for v in indices)
            if not known[indices].all():raise ValueError('Required quartic moment not supplied')
            M[i,ids]+=sign*G[indices]
        add(p==a,(r,q,b,c))
        linear=((r==c)&(q==b)).astype(int)-((r==b)&(q==c)).astype(int)
        ids=np.flatnonzero(linear)
        if not known_one[a,p[ids]].all():raise ValueError('Required quadratic moment not supplied')
        M[i,ids]+=linear[ids]*one[a,p[ids]]
        add(r==c,(a,q,b,p),-1);add(r==b,(a,q,c,p),1)
        add(q==c,(a,r,b,p),1);add(q==b,(a,r,c,p),-1)
    if np.max(abs(M-M.T))>1e-9:raise ValueError('Paired moment matrix is not Hermitian')
    return (M+M.T)/2


def select(case,checkpoint,output):
    started=time.monotonic();prepared=case/'prepared'
    frame=json.loads((prepared/'frame.json').read_text());m=frame['modes']
    T=sparse.load_npz(prepared/'twirl.npz');S=np.load(prepared/'selected.npy')
    y=np.load(checkpoint)['y'];values=-T[S].T@(np.load(prepared/'scale.npy')*y)/np.load(prepared/'weights.npy')
    moments={tuple(map(tuple,w)):float(v) for w,v in zip(frame['rows'],values) if len(w)<=4}
    norm=moments[()]
    if not np.isfinite(norm) or norm<=0:raise ValueError('Positive finite normalization required')
    moments={w:v/norm for w,v in moments.items()}
    spin=lambda k:1 if k%2==0 else -1
    words=[((1,p),(0,q),(0,r)) for p in range(m) for q,r in combinations(range(m),2)
        if spin(p)-spin(q)-spin(r)==1]
    if len(words)>3000:raise ValueError('Declared diagnostic matrix limit exceeded')
    M=operator_moment_matrix(words,moments,m)
    ev,V=linalg.eigh(M,subset_by_index=(0,3),check_finite=False)
    rank=int(np.count_nonzero(ev<-1e-7))
    if not rank:raise ValueError('No violated spin-change-one candidate')
    integers=np.rint(V[:,:rank]*10**7).astype(np.int64)
    result={'kind':'spin_change_one_T2_directions','frame_sha256':hashlib.sha256((prepared/'frame.json').read_bytes()).hexdigest(),
        'checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),'source':str(checkpoint),
        'seconds':time.monotonic()-started,'diagnostic_matrix_dimension':len(words),
        'new_energy_result':False,'exact_dual_obstruction':False,'unavailable_higher_moments_used':False,
        'channels':[{'words':words,'basis_integers':integers.tolist(),'denominator':10**7,
            'dimension':rank,'unverified_eigenvalues':ev[:rank].tolist()}]}
    with output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='channels'}|{'eigenvalues':ev.tolist()}),flush=True)


def solve(case,old_path,new_path,checkpoint,tag,seconds):
    import sys
    from research.acceptance_channels_20260915.campaign import ROOT
    from research.acceptance_channels_20260915.rank_update import RankUpdatedNormal,factorize
    from research.gpu_acceleration_20260915 import solve as base
    from research.gpu_acceleration_20260915.solve import Operator as BaseOperator
    from research.interacting_scaling_20260915.dictionary import representative
    sys.path.insert(0,str(ROOT/'.venv-interacting-libs'))
    from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient
    frame_hash=hashlib.sha256((case/'prepared/frame.json').read_bytes()).hexdigest()
    documents=[json.loads(p.read_text()) for p in (old_path,new_path)]
    if any(d['frame_sha256']!=frame_hash or len(d['channels'])!=1 for d in documents):
        raise ValueError('Collective source family changed')
    if documents[1]['checkpoint_sha256']!=hashlib.sha256(checkpoint.read_bytes()).hexdigest():
        raise ValueError('New directions were selected from a different checkpoint')
    previous=json.loads((checkpoint.parent/'operator_extension.json').read_text())
    if previous['channel_file_sha256']!=hashlib.sha256(old_path.read_bytes()).hexdigest():
        raise ValueError('Parent directions do not match the inherited Gram')
    channels=[d['channels'][0] for d in documents]
    ranks=[c['dimension'] for c in channels]
    if not 1<=ranks[0]<=8 or not 1<=ranks[1]<=4:raise ValueError('Declared twelve-direction envelope exceeded')
    meta=json.loads((case/'prepared/frame.json').read_text());count=len(meta['blocks'])
    raw=np.load(checkpoint);arrays={k:raw[k] for k in raw.files}
    if arrays[f'Q_{count}'].shape!=(ranks[0],ranks[0]) or f'Q_{count+1}' in arrays:
        raise ValueError('One parent collective block required')
    arrays[f'Q_{count+1}']=np.zeros((ranks[1],ranks[1]))
    restart=case/(tag+'_restart.npz')
    if restart.exists():raise FileExistsError(restart)
    np.savez_compressed(restart,**arrays)
    def construct(folder):
        started=time.monotonic();op=BaseOperator(folder);maps=[]
        lookup={tuple(map(tuple,w)):i for i,w in enumerate(op.meta['rows'])}
        T=sparse.load_npz(op.prepared/'twirl.npz');S=np.load(op.prepared/'selected.npy')
        transform=(sparse.diags(op.scale)@T[S]).tocsc();op.meta['pairs']=[]
        for j,channel in enumerate(channels):
            words=[tuple(map(tuple,w)) for w in channel['words']];den=channel['denominator']
            integers=np.asarray(channel['basis_integers'],dtype=np.int64);rank=ranks[j]
            if type(den) is not int or den<=0 or integers.shape!=(len(words),rank) or np.max(abs(integers))>den:
                raise ValueError('Normalized collective direction coefficients required')
            size=op.meta['modes'] if j else op.meta['modes']//2
            basis=[(general_tensor if j else tensors)(words,integers[:,i],size) for i in range(rank)]
            columns=[]
            for a in range(rank):
                for b in range(rank):
                    coefficients={}
                    for word,value in (general_entry if j else contracted_entry)(basis[a],basis[b],size):
                        w,sign=normal_word(word)
                        if w is None:continue
                        if representative(w) not in lookup:raise ValueError('Missing coefficient row')
                        if w in lookup:
                            row=lookup[w];coefficients[row]=coefficients.get(row,0)+sign*value
                    rows=list(coefficients);values=[coefficients[r]/(2*den*den) for r in rows]
                    vector=sparse.csc_matrix((values,(rows,np.zeros(len(rows),dtype=int))),shape=(len(lookup),1))
                    columns.append((transform@vector).tocsc())
            M=sparse.hstack(columns,format='csr');maps.append(M);g=len(op.meta['groups'])
            op.meta['groups'].extend([{'name':f'selected_T2_{j}-','words':words},
                {'name':f'selected_T2_{j}+','words':[tuple((1-c,k) for c,k in reversed(w)) for w in words]}])
            op.meta['pairs'].append({'adjoint_order':list(range(len(words)))})
            op.meta['blocks'].append({'physical_group':g,'dimension':rank,'kind':'selected_paired'})
            op.V.append(integers/den);op.members.append([g,g+1]);op.ids.append(j)
            op.Q.append(np.zeros((rank,rank)));op.M.append(M);op.MT.append(M.T.tocsr())
        op.G=RankUpdatedNormal(op.G,sparse.hstack(maps,format='csc'))
        for j,M in enumerate(maps):sparse.save_npz(folder/tag/f'added_map_{j}.npz',M)
        record={'seconds':time.monotonic()-started,'ranks':ranks,'added_Gram_entries':sum(r*r for r in ranks),
            'operator_terms':[len(c['words']) for c in channels],'coefficient_rows_added':0,
            'complete_global_cubic_map_constructed':False,'normal_outer_product_materialized':False,
            'parent_channel_sha256':hashlib.sha256(old_path.read_bytes()).hexdigest(),
            'new_channel_sha256':hashlib.sha256(new_path.read_bytes()).hexdigest(),
            'source_checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest()}
        (folder/tag/'operator_extension.json').write_text(json.dumps(record,indent=2)+'\n')
        return op
    old_lu=base.splu;base.splu=lambda A:factorize(A,old_lu) if isinstance(A,RankUpdatedNormal) else old_lu(A)
    base.Operator=construct;base.Quotient=SparseQuotient
    base.run(case,tag,seconds,.03,restart,'cpu_evd',None,100,None,True)


if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('select');a.add_argument('case',type=Path);a.add_argument('checkpoint',type=Path);a.add_argument('output',type=Path)
    a=sub.add_parser('solve');a.add_argument('case',type=Path);a.add_argument('parent',type=Path);a.add_argument('channels',type=Path)
    a.add_argument('checkpoint',type=Path);a.add_argument('tag');a.add_argument('--seconds',type=float,default=300)
    args=p.parse_args()
    if args.action=='select':select(args.case.resolve(),args.checkpoint.resolve(),args.output.resolve())
    else:solve(args.case.resolve(),args.parent.resolve(),args.channels.resolve(),args.checkpoint.resolve(),args.tag,args.seconds)
