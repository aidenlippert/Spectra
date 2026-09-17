"""Direct four-word paired PSD blocks with joint reconstruction constraints."""
from fractions import Fraction
from experiments.marginal_symbolic import add, adj, canonical, mono, product, scale


def paired_entry(left, right):
    p, q = mono(left), mono(right)
    value = scale(add(product(adj(p), q), product(adj(q), p),
        product(p, adj(q)), product(q, adj(p))), Fraction(1, 2))
    if max(map(len, value), default=0) > 4 or canonical(adj(value)) != value:
        raise ValueError('Paired Gram entry failed exact quartic cancellation')
    return value


def run(case, channel_path, source, tag, seconds):
    import hashlib
    import json
    from pathlib import Path
    import sys
    import time
    import numpy as np
    from scipy import sparse
    from research.gpu_acceleration_20260915 import solve as base
    from research.gpu_acceleration_20260915.solve import Operator as BaseOperator
    from research.acceptance_channels_20260915.campaign import ROOT
    from research.interacting_scaling_20260915.dictionary import representative
    sys.path.insert(0, str(ROOT/'.venv-interacting-libs'))
    from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient
    channels = json.loads(channel_path.read_text())
    if not channels['channels'] or len(channels['channels']) > 8:
        raise ValueError('Require one to eight declared channels')
    if channels['frame_sha256'] != hashlib.sha256((case/'prepared/frame.json').read_bytes()).hexdigest():
        raise ValueError('Channels refer to a different source dictionary')
    source_raw = np.load(source)
    source_count = len(json.loads((case/'prepared/frame.json').read_text())['blocks'])
    restart = case/(tag+'_restart.npz')
    if restart.exists(): raise FileExistsError(restart)
    restart_arrays = {name: source_raw[name] for name in source_raw.files}
    for i, channel in enumerate(channels['channels']):
        n = len(channel['words'])
        if n not in (4,8,16,32,64): raise ValueError('Unsupported declared block size')
        restart_arrays[f'Q_{source_count+i}'] = np.zeros((n,n))
    np.savez_compressed(restart, **restart_arrays)
    def construct(folder):
        start = time.monotonic()
        op = BaseOperator(folder)
        rows = [tuple(map(tuple,w)) for w in op.meta['rows']]
        lookup = {w:i for i,w in enumerate(rows)}
        T = sparse.load_npz(op.prepared/'twirl.npz')
        selected = np.load(op.prepared/'selected.npy')
        transform = (sparse.diags(op.scale)@T[selected]).tocsc()
        op.meta['pairs'] = []
        additions = []
        for k, channel in enumerate(channels['channels']):
            words = [tuple(map(tuple,w)) for w in channel['words']]
            if len({tuple(sorted(canonical(mono(w)))) for w in words}) != len(words):
                raise ValueError('Repeated channel word')
            n = len(words)
            columns = []
            for i in range(n):
                for j in range(n):
                    polynomial = paired_entry(words[i], words[j])
                    if any(representative(w) not in lookup for w in polynomial):
                        raise ValueError('Channel needs coefficient rows absent from the supplied family')
                    kept = [(lookup[w],float(c)) for w,c in polynomial.items() if w in lookup]
                    vector = sparse.csc_matrix(([c for _,c in kept],([i for i,_ in kept],np.zeros(len(kept),dtype=int))),shape=(len(rows),1))
                    columns.append((transform@vector).tocsc())
            M = sparse.hstack(columns,format='csr')
            g = len(op.meta['groups'])
            op.meta['groups'].extend([
                {'name':f'selected_paired_{k}-','words':words},
                {'name':f'selected_paired_{k}+','words':[tuple((1-c,i) for c,i in reversed(w)) for w in words]}])
            op.meta['pairs'].append({'adjoint_order':list(range(n))})
            op.meta['blocks'].append({'physical_group':g,'dimension':n,'kind':'selected_paired_full_Gram'})
            op.V.append(np.eye(n));op.members.append([g,g+1]);op.ids.append(k)
            op.Q.append(np.zeros((n,n)));op.M.append(M);op.MT.append(M.T.tocsr())
            op.G = (op.G+M@M.T).tocsc()
            sparse.save_npz(folder/tag/f'added_map_{k}.npz',M)
            additions.append({'block':k,'Gram_entries':n*n,'map_nonzeros':M.nnz,
                'coefficient_rows_touched':int(np.count_nonzero(np.diff(M.indptr)))})
        record = {'kind':'direct_paired_channel_extension','seconds':time.monotonic()-start,
            'same_original_model_upper_and_ideals':True,
            'parent_preparation_reused_and_charged':True,
            'source_checkpoint':str(source),'source_checkpoint_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'channel_file_sha256':hashlib.sha256(channel_path.read_bytes()).hexdigest(),
            'new_global_cubic_map_constructed':False,'coefficient_rows_added':0,
            'blocks':additions,'added_Gram_entries':sum(a['Gram_entries'] for a in additions),
            'positive_pairing':'Two adjoint factors share a full PSD Gram; all four-word cross terms retained'}
        (folder/tag/'operator_extension.json').write_text(json.dumps(record,indent=2)+'\n')
        return op
    base.Operator = construct
    base.Quotient = SparseQuotient
    base.run(case,tag,seconds,.03,restart,'cpu_evd',None,100,None,True)


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('channels',type=Path)
    p.add_argument('source',type=Path);p.add_argument('tag');p.add_argument('--seconds',type=float,default=300)
    a=p.parse_args();run(a.case.resolve(),a.channels.resolve(),a.source.resolve(),a.tag,a.seconds)
