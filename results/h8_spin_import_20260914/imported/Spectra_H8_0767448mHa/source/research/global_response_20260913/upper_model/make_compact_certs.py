import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from research.global_response_20260913.slater import make_certificate, check

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'results/global_response_20260913/upper_model'

def main():
    trials=json.loads((OUT/'compact_trials.json').read_text())
    specs=[
      ('h6','results/certificate_scaling/active_space_ladder/h6/fixture.json',3,3,None),
      ('h6_1p6','results/response_consistency_20260913/fresh_h6_1p6/fixture.json',3,3,None),
      ('h8','results/certificate_scaling/active_space_ladder/h8/fixture.json',4,4,None),
      ('h6_1p73','results/global_response_20260913/fresh_h6_1p73/fixture.json',3,3,None),
    ]
    out=[]
    for name,path,na,nb,spin in specs:
        d=json.loads((ROOT/path).read_text()); t=next((x for x in trials if x['case']==name),None)
        if t is None: aa=bb=['0']*(d['modes']//4)
        else:
            n=d['modes']//4; vals=t['tangents_1e8']; aa=[str(x)+'/100000000' for x in vals[:n]]; bb=[str(x)+'/100000000' for x in vals[2*n:3*n]]
        cert=make_certificate(d,aa,bb,na,nb,spin); check(d,cert)
        (OUT/(name+'_slater_certificate.json')).write_text(json.dumps(cert,indent=2)+'\n')
        out.append({'case':name,'certificate':cert,'fixture':path})
    # CH2 model: common zero rotations provide exact pure-spin controls.
    d=json.loads((ROOT/'results/ch2_validation_20260913/fixture.json').read_text()); n=d['modes']//4
    for label,na,nb,spin in (('ch2_s0',3,3,0),('ch2_s1',4,2,1)):
        z=['0']*n; cert=make_certificate(d,z,z,na,nb,spin); check(d,cert)
        (OUT/(label+'_slater_certificate.json')).write_text(json.dumps(cert,indent=2)+'\n'); out.append({'case':label,'certificate':cert,'fixture':'results/ch2_validation_20260913/fixture.json'})
    (OUT/'compact_certificates.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{'case':x['case'],'upper_Ha':x['certificate']['upper_Ha'],'S2':x['certificate']['target_spin']} for x in out],indent=2))
if __name__=='__main__': main()
