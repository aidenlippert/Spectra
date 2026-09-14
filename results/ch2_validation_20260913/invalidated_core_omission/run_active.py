import json, time
from pathlib import Path
import numpy as np
from pyscf import gto, scf, fci, ao2mo

OUT=Path(__file__).parents[2]/'results/ch2_validation_20260913'
geom='C 0 0 0; H 0 0 1.117; H 0 1.047 -0.389'
mol=gto.M(atom=geom,basis='sto-3g',unit='Angstrom',charge=0,spin=0,verbose=0)
t=time.monotonic(); mf=scf.RHF(mol).run(conv_tol=1e-12,verbose=0)
mo=mf.mo_coeff; norb=6; core=mo[:,0]; active=mo[:,1:7]; h=active.T@mf.get_hcore()@active
eri=ao2mo.restore(1,ao2mo.kernel(mol,active),norb)
solver=fci.direct_spin1.FCI(mol); solver.conv_tol=1e-12
e,ci=solver.kernel(h,eri,norb,(3,3),nroots=20)
rows=[]
for i,x in enumerate(np.atleast_1d(e)):
    ss,mult=fci.spin_op.spin_square(np.asarray(ci[i]),norb,(3,3))
    rows.append({'root':i,'energy_elec':float(x),'S2':float(ss),'spin':round((np.sqrt(1+4*ss)-1)/2,8),'multiplicity':int(round(mult))})
ssolver=fci.direct_spin0.FCI(mol); ssolver.conv_tol=1e-12
es,_=ssolver.kernel(h,eri,norb,(3,3),nroots=1)
tesolver=fci.direct_spin1.FCI(mol); tesolver.conv_tol=1e-12
et,_=tesolver.kernel(h,eri,norb,(4,2),nroots=1)
sing={'root':'spin0_solver','energy_elec':float(np.atleast_1d(es)[0]),'S2':0.0,'spin':0.0,'multiplicity':1}
trip={'root':'Ms1_solver','energy_elec':float(np.atleast_1d(et)[0]),'S2':2.0,'spin':1.0,'multiplicity':3}
core_const=float(2*core@mf.get_hcore()@core + 2*mol.intor('int1e_ovlp').trace()*0)
manifest={'geometry':geom,'basis':'STO-3G','active_orbitals':'RHF canonical spatial orbitals 1:7; one spatial core orbital frozen','active_electrons':'6 (3 alpha, 3 beta)','observable':'Delta_ST=E_T-E_S','norb':norb,'active_dimension_Ms0':400,'nuclear_repulsion':float(mol.energy_nuc()),'rhf_total_energy':float(mf.e_tot),'core_orbital':'0','core_constant_onebody_only':core_const,'roots':rows,'singlet':sing,'triplet':trip,'gap_ET_minus_ES':trip['energy_elec']-sing['energy_elec'],'wall_seconds':time.monotonic()-t,'pyscf_version':__import__('pyscf').__version__,'status':'corrected_neutral_CH2_one_core_six_active_numerical_control'}
(OUT/'active_space_result.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
