"""Exact spin-resolved diagonal overlap obstruction after charge-law closure."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,sys
from coherent_projector_overlap import ROOT,reconstruct,matrix_digest
from full_overlap_telescope import hole5,spin5
from experiments.marginal_local_hubbard_block import _reflection

def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    directory=parser.parse_args().directory.resolve();cp=directory/'range_two_family_limit_certificate.json'
    rp=directory/'overlap/full_overlap_replay.json';receipt=json.loads(rp.read_text());family=json.loads(cp.read_text())
    if not receipt.get('accepted') or not receipt.get('stationary_extension_refuted'):raise ValueError('Accepted full-overlap obstruction required')
    for name,h in receipt['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=h:raise ValueError('Stale full-overlap source')
    denominator,left,right,difference,_=reconstruct(family['mixture'])
    if matrix_digest(difference)!=receipt['difference_integer_matrix_sha256']:raise ValueError('Reconstructed density matrix differs')
    diagonal={s:difference.get((s,s),0) for s in range(1024)}
    if sum(diagonal.values()):raise ValueError('Diagonal overlap difference has nonzero trace')
    positive={s for s,a in diagonal.items() if a>0};negative={s for s,a in diagonal.items() if a<0}
    if not positive:raise ValueError('No diagonal obstruction')
    for action in (hole5,spin5):
        if {action(s)[0] for s in positive}!=positive:raise ValueError('Positive projector is not PH/spin invariant')
    if {_reflection(s,5)[0] for s in positive}!=negative:raise ValueError('Reflection does not exchange diagonal signs')
    probabilities=[F(sum(matrix.get((s,s),0) for s in positive),denominator) for matrix in (left,right)]
    variation=probabilities[0]-probabilities[1]
    if any(not 0<=p<=1 for p in probabilities) or variation!=F(sum(abs(a) for a in diagonal.values()),2*denominator):raise ValueError('Positive diagonal projector probabilities failed')
    y=lambda s:F(int(s in positive)-int(s in negative),2)
    direct=F(0)
    for item in family['mixture']:
        vector={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in vector.values())
        direct+=F(item['weight'])*sum((a*a*(y(s&1023)-y(s>>2)) for s,a in vector.items()),F(0))/norm
    if direct!=variation:raise ValueError('Direct spin-word telescope differs from diagonal overlap')
    charge=defaultdict(int)
    for s,a in diagonal.items():charge[tuple(((s>>(2*i))&3).bit_count()-1 for i in range(5))]+=a
    if any(charge.values()):raise ValueError('Coarse charge-law overlap did not close')
    files={Path(__file__).resolve(),cp,rp,Path(sys.modules['coherent_projector_overlap'].__file__).resolve(),Path(sys.modules['full_overlap_telescope'].__file__).resolve(),Path(sys.modules['full_overlap_density'].__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name and str(Path(name).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(name).resolve())
    result={'accepted':True,'coarse_charge_word_overlap_exactly_zero':True,'spin_resolved_word_overlap_refuted':True,
            'positive_projector_states':sorted(positive),'negative_difference_states':len(negative),
            'left_probability':str(probabilities[0]),'right_probability':str(probabilities[1]),
            'exact_diagonal_total_variation':str(variation),'diagonal_total_variation_float':float(variation),
            'exact_original_mixture_telescope_expectation':str(direct),'telescope_operator_norm_bound':1,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
            'scope':'Positive diagonal projector detects unequal overlapping spin-resolved word laws even though every coarse charge-word overlap is zero. Y=(Ppositive-Pnegative)/2 and T=Yleft-Yright obey ||T||<=1 and stationary translated sumzero. No energy improvement or general representability conclusion.'}
    (directory/'overlap/diagonal_overlap_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'accepted':True,'positive_states':len(positive),'total_variation':float(variation),'charge_overlap_zero':True}))

if __name__=='__main__':main()
