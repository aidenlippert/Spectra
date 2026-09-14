"""Fresh exact proof that spin certificates exceed the older hopping family."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,sys
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_range_two_family_limit import replay


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);parser.add_argument('old_directory',type=Path);args=parser.parse_args();folder=args.directory.resolve();old=args.old_directory.resolve()
    cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json';ep=folder/'profile_joint_r1_2_certificate.json';er=folder/'range_two_replay.json';op=old/'range_two_family_limit_certificate.json'
    current=json.loads(cp.read_text());receipt=json.loads(rp.read_text());energy=json.loads(ep.read_text());previous=json.loads(op.read_text())
    if not receipt['accepted'] or not receipt.get('spectator_hopping') or current['kind']!='joint_diagonal_range2_family_limit_v9' or previous['kind']!='joint_diagonal_range2_family_limit_v8':raise ValueError('Matched spectator family and older spin family required')
    for path in [cp,ep,er]:
        if receipt['source_sha256'][str(path.relative_to(ROOT))]!=hashlib.sha256(path.read_bytes()).hexdigest():raise ValueError('Current matching receipt is stale')
    for key in ['half_vector','charged_vector','ratio','theta_half','theta_joint','diagonal_shapes','W','range_two_density_profile']:
        if current[key]!=previous[key]:raise ValueError('Compared fixed families differ beyond added spectator moments')
    if not energy.get('spectator_hopping'):raise ValueError('Actual spectator energy correction required')
    # Replay the older physical mixture using current exact production code.
    # No reliance on an old source hash or numerical ceiling value.
    old_result=replay(previous);ceiling=F(old_result['periodic_family_upper']);lower=F(receipt['accepted_periodic_lower']);escape=lower-ceiling
    if escape<=0:raise ValueError('No strict escape from the older family ceiling')
    files={Path(__file__).resolve(),cp,rp,ep,er,op}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'new_periodic_lower':str(lower),'old_periodic_family_ceiling':str(ceiling),'strict_escape':str(escape),'strict_escape_float':float(escape),'old_family_fresh_replay':old_result,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Strict separation from the entire fixed spin+hopping+complete-charge correction family, including unrestricted coefficients and reflected mean-correct profiles, at the same projector sources/ratio/ceilings/sparse span. Fresh exact replay of the old dual plus the independently accepted new lower. Does not bound families with changed supports or sources and is not a general representability theorem.'}
    (folder/'family_escape.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'accepted':True,'strict_escape_float':float(escape)}))


if __name__=='__main__':main()
