"""Independent standard-library physical-mixture replay of the family ceiling."""
from pathlib import Path
import argparse,hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_joint_family_limit import replay


def main():
    parser=argparse.ArgumentParser();parser.add_argument('certificate',type=Path);a=parser.parse_args()
    start=time.monotonic();result=replay(json.loads(a.certificate.read_text()))
    files={Path(__file__).resolve(),a.certificate.resolve()}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result.update(seconds=time.monotonic()-start,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)})
    a.certificate.with_name('family_limit_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'}),flush=True)


if __name__=='__main__':main()
