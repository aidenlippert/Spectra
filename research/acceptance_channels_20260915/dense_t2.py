"""Paired cubic directions via quartic contractions and a low-rank normal update."""
from itertools import combinations


def contracted_terms(D, E, F, spatial):
    """Emit integer numerator terms, before normal ordering mixed spin labels."""
    pairs=list(combinations(range(spatial),2))
    for p in range(spatial):
        for s in range(spatial):
            if D[p,s]: yield ((1,2*p),(0,2*s)), int(D[p,s])
    for i,(q,r) in enumerate(pairs):
        for j,(t,u) in enumerate(pairs):
            if E[i,j]: yield ((1,2*t+1),(1,2*u+1),(0,2*q+1),(0,2*r+1)), -int(E[i,j])
    for p in range(spatial):
        for w in range(spatial):
            for s in range(spatial):
                for v in range(spatial):
                    coefficient=-int(F[p*spatial+w,s*spatial+v])
                    if coefficient:
                        yield ((1,2*p),(1,2*v+1),(0,2*w+1),(0,2*s)),coefficient


def tensors(words, coefficients, spatial):
    import numpy as np
    from operator import index
    if type(spatial) is not int or spatial < 2 or len(words) != len(coefficients):
        raise ValueError('Consistent word and coefficient dimensions required')
    pairs=list(combinations(range(spatial),2));lookup={v:i for i,v in enumerate(pairs)}
    # Accumulate with Python integers before allocating fixed-width arrays.
    # In particular, repeated words must not overflow before the dot-product check.
    accumulated={}
    for word,coefficient in zip(words,coefficients):
        if len(word)!=3 or [c for c,_ in word]!=[1,0,0]: raise ValueError('Mixed cubic words required')
        p,q,r=[i for _,i in word]
        if p%2 or q%2!=1 or r%2!=1 or q==r or min(p,q,r)<0 or max(p,q,r)>=2*spatial:
            raise ValueError('Require one alpha creator and two distinct beta annihilators')
        p,q,r=p//2,q//2,r//2
        try:
            if isinstance(coefficient,(bool,np.bool_)):raise TypeError
            value=index(coefficient)
        except TypeError as error:
            raise ValueError('Exact integer contraction coefficients required') from error
        if q>r:q,r,value=r,q,-value
        key=p,q,r
        accumulated[key]=accumulated.get(key,0)+value
    maximum=max(map(abs,accumulated.values()),default=0)
    if 2*int(maximum)**2*max(len(pairs),spatial)>=2**63:
        raise ValueError('Integer contraction would overflow; no unchecked arithmetic')
    C=np.zeros((spatial,len(pairs)),dtype=np.int64)
    A=np.zeros((spatial,spatial,spatial),dtype=np.int64)
    for (p,q,r),value in accumulated.items():
        C[p,lookup[(q,r)]]=value
        A[p,q,r]=value;A[p,r,q]=-value
    return C,A.reshape(spatial*spatial,spatial)


def contracted_entry(left,right,spatial):
    C,A=left;D,B=right
    return contracted_terms(C@D.T+D@C.T,C.T@D+D.T@C,A@B.T+B@A.T,spatial)


def normal_word(word):
    # All emitted terms are already creation-before-annihilation ordered.
    creators=[i for c,i in word if c];annihilators=[i for c,i in word if not c]
    if len(set(creators))<len(creators) or len(set(annihilators))<len(annihilators):return None,0
    inversions=sum(a>b for seq in (creators,annihilators) for k,a in enumerate(seq) for b in seq[k+1:])
    return tuple((1,i) for i in sorted(creators))+tuple((0,i) for i in sorted(annihilators)),(-1)**inversions


def embed_parent_gram(basis,parent,old,rank):
    import numpy as np
    inherited=np.asarray(parent['basis_integers'],dtype=np.int64)
    count=inherited.shape[1]
    if not np.array_equal(basis[:,:count],inherited) or old.shape!=(count,count) or count>=rank:
        raise ValueError('Nested collective basis or parent Gram changed')
    result=np.zeros((rank,rank));result[:count,:count]=old
    return result


def run(case,channel_path,source,tag,seconds,parent_path=None):
    import hashlib,json,sys,time
    import numpy as np
    from scipy import sparse
    from scipy.sparse.linalg import LinearOperator
    from research.gpu_acceleration_20260915 import solve as base
    from research.gpu_acceleration_20260915.solve import Operator as BaseOperator
    from research.acceptance_channels_20260915.campaign import ROOT
    from research.interacting_scaling_20260915.dictionary import representative
    from research.acceptance_channels_20260915.rank_update import RankUpdatedNormal, factorize
    sys.path.insert(0,str(ROOT/'.venv-interacting-libs'))
    from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient
    definition=json.loads(channel_path.read_text())
    if len(definition['channels'])!=1 or 'basis_integers' not in definition['channels'][0]:
        raise ValueError('One declared dense collective block required')
    channel=definition['channels'][0];words=[tuple(map(tuple,w)) for w in channel['words']]
    den=channel['denominator']
    if type(den) is not int or den<=0:raise ValueError('Positive integer denominator required')
    supplied=channel['basis_integers']
    if not supplied or not all(isinstance(row,list) and 1<=len(row)<=12 and
            all(type(v) is int and abs(v)<=den for v in row) for row in supplied):
        raise ValueError('Normalized integer direction columns required')
    V=np.asarray(supplied,dtype=np.int64);rank=V.shape[1]
    if not 1<=rank<=12 or len(V)!=len(words):raise ValueError('Declared rank envelope exceeded')
    meta=json.loads((case/'prepared/frame.json').read_text())
    if definition['frame_sha256']!=hashlib.sha256((case/'prepared/frame.json').read_bytes()).hexdigest():
        raise ValueError('Channel source family changed')
    raw=np.load(source); arrays={k:raw[k] for k in raw.files}
    block_key=f'Q_{len(meta["blocks"])}'
    if parent_path:
        if channel.get('parent_file_sha256')!=hashlib.sha256(parent_path.read_bytes()).hexdigest():
            raise ValueError('Nested collective parent hash failed')
        parent=json.loads(parent_path.read_text())['channels'][0]
        if parent['words']!=channel['words'] or parent['denominator']!=den or block_key not in arrays:
            raise ValueError('Parent dictionary or Gram is missing')
        arrays[block_key]=embed_parent_gram(V,parent,arrays[block_key],rank)
    else:
        if block_key in arrays or channel.get('parent_dimension',0):
            raise ValueError('Existing collective block requires its declared parent')
        arrays[block_key]=np.zeros((rank,rank))
    restart=case/(tag+'_restart.npz')
    if restart.exists():raise FileExistsError(restart)
    np.savez_compressed(restart,**arrays)
    def construct(folder):
        start=time.monotonic();op=BaseOperator(folder)
        s=op.meta['modes']//2
        lookup={tuple(map(tuple,w)):i for i,w in enumerate(op.meta['rows'])}
        T=sparse.load_npz(op.prepared/'twirl.npz');selected=np.load(op.prepared/'selected.npy')
        transform=(sparse.diags(op.scale)@T[selected]).tocsc()
        tensor_basis=[tensors(words,V[:,i],s) for i in range(rank)]
        columns=[]
        for i in range(rank):
            for j in range(rank):
                coefficients={}
                for word,value in contracted_entry(tensor_basis[i],tensor_basis[j],s):
                    w,sign=normal_word(word)
                    if w is None:continue
                    if representative(w) not in lookup:raise ValueError('Missing required quartic row')
                    if w in lookup:
                        row=lookup[w];coefficients[row]=coefficients.get(row,0)+sign*value
                rows=list(coefficients);values=[coefficients[r]/(2*den*den) for r in rows]
                vector=sparse.csc_matrix((values,(rows,np.zeros(len(rows),dtype=int))),shape=(len(lookup),1))
                columns.append((transform@vector).tocsc())
        M=sparse.hstack(columns,format='csr');g=len(op.meta['groups'])
        op.meta['groups'].extend([{'name':'dense_selected_t2-','words':words},
            {'name':'dense_selected_t2+','words':[tuple((1-c,k) for c,k in reversed(w)) for w in words]}])
        op.meta['pairs']=[{'adjoint_order':list(range(len(words)))}]
        op.meta['blocks'].append({'physical_group':g,'dimension':rank,'kind':'dense_selected_paired'})
        op.V.append(V/den);op.members.append([g,g+1]);op.ids.append(0)
        op.Q.append(np.zeros((rank,rank)));op.M.append(M);op.MT.append(M.T.tocsr())
        op.G=RankUpdatedNormal(op.G,M.tocsc())
        sparse.save_npz(folder/tag/'added_map.npz',M)
        record={'kind':'contracted_dense_paired_extension','seconds':time.monotonic()-start,
            'operator_terms':len(words),'collective_directions':rank,'added_Gram_entries':rank*rank,
            'added_map_nonzeros':M.nnz,'coefficient_rows_added':0,
            'normal_update_outer_product_materialized':False,
            'complete_global_cubic_map_constructed':False,
            'source_checkpoint_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'parent_channel_sha256':hashlib.sha256(parent_path.read_bytes()).hexdigest() if parent_path else None,
            'parent_dimension_preserved':channel.get('parent_dimension',0),
            'channel_file_sha256':hashlib.sha256(channel_path.read_bytes()).hexdigest(),
            'same_model_upper_and_ideals':True,'parent_discovery_and_preparation_additional':True}
        (folder/tag/'operator_extension.json').write_text(json.dumps(record,indent=2)+'\n')
        return op
    old_splu=base.splu
    base.splu=lambda matrix:factorize(matrix,old_splu) if isinstance(matrix,RankUpdatedNormal) else old_splu(matrix)
    base.Operator=construct;base.Quotient=SparseQuotient
    base.run(case,tag,seconds,.03,restart,'cpu_evd',None,100,None,True)


if __name__=='__main__':
    import argparse
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('channels',type=Path)
    p.add_argument('source',type=Path);p.add_argument('tag');p.add_argument('--seconds',type=float,default=300)
    p.add_argument('--parent',type=Path)
    a=p.parse_args();run(a.case.resolve(),a.channels.resolve(),a.source.resolve(),a.tag,a.seconds,
        a.parent.resolve() if a.parent else None)
