"""Standard-library replay of this bounded hidden-basis research campaign."""
from pathlib import Path
import hashlib,json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import decode
from research.certificate_scaling.hidden_density_basis import replay,replay_rank
from research.certificate_scaling.commutant_obstruction import obstruction

def main():
    start=time.monotonic();root=Path('results/certificate_scaling/hidden_density_basis')
    remote=Path('results/lambda_runs/hidden_density_basis');controls=[];ranks=[];bounds=[]
    paths=list(root.glob('M*.json'))+list(remote.glob('hidden_campaign_*/control_*/*[0-9].json'))
    for p in sorted(paths):
        cert=json.loads(p.read_text());rec=replay(cert)
        assert rec['width']=='0'
        controls.append({'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),**rec})
    for p in sorted(remote.glob('hidden_campaign_*/molecule_h*.json')):
        receipt=json.loads(p.read_text());fixture=Path(receipt['fixture'])
        assert hashlib.sha256(fixture.read_bytes()).hexdigest()==receipt['fixture_sha256']
        f=json.loads(fixture.read_text());h=decode(f['hamiltonian'],f['modes'],4)
        r=replay_rank(h,receipt);assert r['rules_out_real_density_basis']
        ranks.append({'fixture':str(fixture),**r})
    for p in sorted(root.glob('quotient_obstruction_h*.json')):
        saved=json.loads(p.read_text());atoms=saved['modes']//2
        fixture=Path(f'results/certificate_scaling/active_space_ladder/h{atoms}/fixture.json')
        assert hashlib.sha256(fixture.read_bytes()).hexdigest()==saved['fixture_sha256']
        f=json.loads(fixture.read_text());h=decode(f['hamiltonian'],f['modes'],4)
        checked=obstruction(h,f['modes'],complex_unitary=True,quotient_number_ideal=True)
        for k in ['inertia','imaginary_antisymmetric_inertia','total_negative_inertia','frobenius_distance_squared_lower','frobenius_distance_lower']:
            assert saved[k]==checked[k],(p,k)
        bounds.append({'fixture':str(fixture),**checked});print(json.dumps({'replayed':str(p),'bound':checked['frobenius_distance_lower'],'seconds':checked['total_seconds']}),flush=True)
    assert len(controls)==7 and len(ranks)==4 and len(bounds)==4
    sources=['hidden_density_basis.py','commutant_obstruction.py','hidden_campaign_replay.py','test_hidden_density_basis.py','test_commutant_obstruction.py','fermionic_ratio_chain.py']
    hashes={s:hashlib.sha256(Path('research/certificate_scaling',s).read_bytes()).hexdigest() for s in sources}
    result={'stdlib_only':True,'controls':controls,'rank_certificates':ranks,'quotient_obstructions':bounds,'source_hashes':hashes,'replay_seconds':time.monotonic()-start}
    (root/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'controls_replayed':len(controls),'rank_certificates_replayed':len(ranks),'quotient_obstructions_replayed':len(bounds),'seconds':result['replay_seconds']}),flush=True)

if __name__=='__main__':main()
