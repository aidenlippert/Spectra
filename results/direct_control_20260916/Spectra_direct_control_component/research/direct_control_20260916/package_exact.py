"""Package the exact component with its actually imported local dependencies."""
import argparse,hashlib,json,shutil,sys,zipfile
from pathlib import Path
from . import collective_control,collective_obstruction,transfer_collective,tensor_operator
from . import test_tensor_operator,test_collective_control


def run(destination):
    repo=Path(__file__).resolve().parents[2];destination=destination.resolve()
    if destination.exists():raise FileExistsError(destination)
    destination.mkdir(parents=True)
    # Preserve original module bytes; no embedded replacement verifier or fake
    # import stubs are used. Namespace packages require no generated init files.
    sources=set()
    for module in tuple(sys.modules.values()):
        filename=getattr(module,'__file__',None)
        if filename:
            path=Path(filename).resolve()
            if repo in path.parents and path.suffix=='.py':sources.add(path)
    notes=repo/'research/direct_control_20260916'
    for path in sources:
        target=destination/path.relative_to(repo);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
    shutil.copyfile(notes/'bundle_replay.py',destination/'replay.py')
    for name in ('REPORT.md','DERIVATION.md','PROTOCOL.md','ARCHIVE_AUDIT.md'):
        if (notes/name).exists():shutil.copyfile(notes/name,destination/name)
    root=repo/'results/direct_control_20260916'
    cases=[('h8_original',root/'imported/Spectra_control_reduction/inputs/fixture.json',root/'imported/Spectra_control_reduction/inputs/state.json',root/'collective_control_bound64.json')]
    for name in ('h6_asymmetric','water_asymmetric'):
        old=repo/'results/transfer_solver_20260915/cases'/name
        cases.append((name,old/'fixture.json',old/'mps/state.json',root/'collective_transfer'/(name+'.json')))
    for name,fixture,state,artifact in cases:
        target=destination/'cases'/name;target.mkdir(parents=True)
        for src,out in ((fixture,'fixture.json'),(state,'state.json'),(artifact,'certificate.json')):shutil.copyfile(src,target/out)
    (destination/'README.md').write_text('# Exact strong-control component\n\nRun `python3 -B -S replay.py --out ../fresh-direct-control-replay` from this directory. The output directory must be new and outside the bundle. Python standard library only.\n\nThis certifies a separate strong-control regime on H8 and two transfer inputs. It does not solve the original amplitude-0.5 task or the general many-body problem. Read REPORT.md and DERIVATION.md before interpreting the numbers. Source Hamiltonian/MPS discovery is inherited. No original full-sector embeddings or teacher trajectories are included.\n')
    files={str(p.relative_to(destination)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(destination.rglob('*')) if p.is_file()}
    (destination/'MANIFEST.json').write_text(json.dumps({'kind':'sha256_integrity_manifest_v1','files_sha256':files},indent=2)+'\n')
    zip_path=destination.with_suffix('.zip')
    if zip_path.exists():raise FileExistsError(zip_path)
    with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(destination.rglob('*')):
            if p.is_file():z.write(p,arcname=str(Path(destination.name)/p.relative_to(destination)))
    print(json.dumps({'bundle':str(zip_path),'bytes':zip_path.stat().st_size,'sha256':hashlib.sha256(zip_path.read_bytes()).hexdigest(),'source_files':len(sources)}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--destination',type=Path,required=True);a=p.parse_args();run(a.destination)
