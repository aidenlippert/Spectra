"""Compare a constructed map with direct polynomial evaluation on the same dual."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from research.acceptance_channels_20260915.dense_t2 import tensors,contracted_entry,normal_word
from research.interacting_scaling_20260915.dictionary import representative


def run(case,checkpoint,channels,mapfile,output):
    frame=json.loads((case/'prepared/frame.json').read_text());T=sparse.load_npz(case/'prepared/twirl.npz')
    S=np.load(case/'prepared/selected.npy');y=np.load(checkpoint)['y']
    raw=-T[S].T@(np.load(case/'prepared/scale.npy')*y)/np.load(case/'prepared/weights.npy')
    moments={tuple(map(tuple,w)):float(v) for w,v in zip(frame['rows'],raw)}
    c=json.loads(channels.read_text())['channels'][0];V=np.array(c['basis_integers'],dtype=np.int64)
    words=[tuple(map(tuple,w)) for w in c['words']];den=c['denominator'];s=frame['modes']//2
    basis=[tensors(words,V[:,i],s) for i in range(V.shape[1])]
    expected=np.zeros((len(basis),len(basis)))
    for i,a in enumerate(basis):
        for j,b in enumerate(basis):
            for word,value in contracted_entry(a,b,s):
                w,sign=normal_word(word)
                if w is not None:expected[i,j]+=sign*value/(2*den*den)*moments[representative(w)]
    supplied=(-sparse.load_npz(mapfile).T@y).reshape(expected.shape)
    error=float(np.max(abs(expected-supplied)))
    record={'maximum_absolute_difference':error,'expected':expected.tolist(),'map':supplied.tolist(),
        'consistent_to_1e_9':error<1e-9,'source_checkpoint':str(checkpoint)}
    with output.open('x') as stream:json.dump(record,stream,indent=2);stream.write('\n')
    print(json.dumps(record),flush=True)
    if error>=1e-9:raise ValueError('Constructed map differs from direct polynomial evaluation')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('case','checkpoint','channels','mapfile','output'):p.add_argument(name,type=Path)
    a=p.parse_args();run(a.case,a.checkpoint,a.channels,a.mapfile,a.output)
