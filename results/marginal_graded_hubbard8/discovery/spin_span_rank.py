"""Exact independence of four spin corrections from the previous family span."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_spin_telescope import LABELS,actions
from experiments.marginal_local_hubbard_block import _sector
from experiments.marginal_polynomial_sos import integer_psd
from experiments.marginal_charged_projectors import charged_vectors


def main():
    cp=ROOT/'results/marginal_graded_hubbard8/spin_telescope/W_zero/polished/profile_joint_r1_2_certificate.json'
    c=json.loads(cp.read_text());charged,_=charged_vectors(c['joint']['vector'])
    sectors={_sector(int(s),6) for s,a in c['vector'].items() if a}|{_sector(s,6) for vector in charged for s,a in vector.items() if a}
    if sectors!={(3,3),(2,3),(3,2),(3,4),(4,3)}:raise ValueError('Expected fixed half/charged source sectors')
    data=[actions({key:1}) for key in LABELS]
    gram=[[sum(a*right[s].get(t,0) for s,row in enumerate(left) if _sector(s,6)==(1,1) for t,a in row.items() if (s^t).bit_count()==4) for right in data] for left in data]
    psd=integer_psd(gram)
    if psd['rank']!=4:raise ValueError('Four independent spin directions required')
    files={Path(__file__).resolve(),cp}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'sector':[1,1],'projector_sectors':[list(key) for key in sorted(sectors)],'changed_bits':4,'double_spin_flip_gram':gram,'exact_psd':psd,'labels':list(LABELS),'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Four independent spin telescope directions on double-spin-flip entries in the two-electron (1,1) sector. All previous diagonal corrections vanish off-diagonal; hopping changes only two bits. The fixed half/charged projectors occupy sectors of total particle number6 or5/7 and vanish here. Therefore these four directions are independent modulo the previous fixed-source certificate span, including penalty variations; no independence claim for arbitrary changed projector sources.'}
    (ROOT/'results/marginal_graded_hubbard8/spin_telescope/spin_rank.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'accepted':True,'gram':gram,'rank':4}))


if __name__=='__main__':main()
