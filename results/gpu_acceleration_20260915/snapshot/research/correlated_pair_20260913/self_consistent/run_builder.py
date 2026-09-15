import argparse,json
from pathlib import Path
from .builder import build_from_path, paired_cubic_certificate

def main():
    p=argparse.ArgumentParser(); p.add_argument('--fixture',type=Path,required=True); p.add_argument('--out',type=Path,required=True); p.add_argument('--paired',action='store_true')
    a=p.parse_args()
    if a.paired:
        f=json.loads(a.fixture.read_text()); cert,rec=paired_cubic_certificate(f); a.out.mkdir(parents=True,exist_ok=True); (a.out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n'); (a.out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n'); print(json.dumps(rec,indent=2))
    else: print(json.dumps(build_from_path(a.fixture,a.out),indent=2))
if __name__=='__main__': main()
