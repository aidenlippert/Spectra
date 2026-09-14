"""Exact replay of the bounded six-site nonuniform window proposal."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_local_energy_chain import replay_many
from experiments.marginal_window_family_bound import replay as replay_family
OUT=ROOT/'results/marginal_graded_hubbard8/weighted_window_family'


def main():
    OUT.mkdir(exist_ok=True)
    window={'kind':'local_hubbard_block_v1','sites':6,'U':'10/3','t':'1',
            'onsite_profile':['3/20','83/25','653/100','653/100','83/25','3/20'],
            'hopping_profile':['5/12','5/4','5/3','5/4','5/12'],'lower':'-3.05817'}
    (OUT/'six_site_certificate.json').write_text(json.dumps(window,indent=2)+'\n')
    old=json.loads((ROOT/'results/marginal_graded_hubbard8/local_energy_chain/certificate.json').read_text())
    chains=[]
    for N in (12,64,1000000):
        r=N%6
        blocks=[old['local']['physical6']]+([old['local'][f'physical{r}']] if r else [])
        chains.append({'kind':'local_energy_chain_v1','sites':N,'block_length':6,
                       'U':'4','t':'1','blocks':blocks,'overlap':window})
    start=time.monotonic();result=replay_many(chains);result['seconds']=time.monotonic()-start
    sources=['experiments/marginal_local_hubbard_block.py','experiments/marginal_local_energy_chain.py',
             'experiments/marginal_window_family_bound.py',
             'results/marginal_graded_hubbard8/discovery/six_weighted_window.py']
    result['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}
    result['scope_note']='Exact acceptance of this rational six-site profile; no optimality claim for the larger profile family.'
    (OUT/'six_site_independent_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    family_certificate=json.loads((OUT/'certificate.json').read_text())
    family=replay_family(family_certificate)
    if not family['overlap_consistency']['three_site_reductions_equal']:
        raise ValueError('Candidate pair data are not stationary across the block')
    candidate=F(family['arbitrary_three_site_boundary_correction_upper'])/3
    local_lower=F(result['chains'][0]['overlap_replay']['lower'])/5
    vector={int(s):a for s,a in family_certificate['upper_vector'].items()}
    norm=sum(a*a for a in vector.values())
    fillings=[F(sum(a*a*(((s>>(2*i))&3).bit_count()) for s,a in vector.items()),norm) for i in range(4)]
    if fillings!=[1]*4 or local_lower<=candidate:
        raise ValueError('Exact stationary pair-data separation failed')
    separation={'accepted':True,'candidate_bulk_energy_per_site':str(candidate),
        'candidate_site_fillings':list(map(str,fillings)),
        'certified_periodic_energy_lower_per_site':str(local_lower),
        'strict_separation_per_site':str(local_lower-candidate),'gap_float':float(local_lower-candidate),
        'source_sha256':result['source_sha256'],
        'scope':'Fresh exact family and six-site local replay. Repeating the limiting four-site witness adjacent-pair data on a periodic Hubbard chain withN>6 and mean filling1 violates the certified global positive-operator energy inequality. The data remain physical on the original four-site block.'}
    (OUT/'periodic_pair_marginal_separator.json').write_text(json.dumps(separation,indent=2)+'\n')
    for row in result['chains']:
        print(json.dumps({'sites':row['sites'],'lower_density':float(F(row['lower_per_site'])),
                          'upper_density':float(F(row['upper_per_site'])),'width_density':float(F(row['width_per_site']))}),flush=True)
    print(json.dumps({'accepted':True,'seconds':result['seconds']}),flush=True)

if __name__=='__main__':main()
