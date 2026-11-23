"""All three ligand types must build from ONE session with no chemistry warnings."""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
from maxi_eval import ligands as L, structures as S, tmol_core as TC, protocol as P
T = L.build_templates(Path('ligands'))
rep = S.extract_replicates(Path('work/inputs/4BS0/tsa/4BS0_tsa_A1305.pdb'),'4BS0','tsa','crystal',T)[0]
prot = rep.protein[rep.protein.element!='H']
sess = TC.Session.create(Path('ligands'), prot, T)
for arm,resname in (('tsa','6NT'),('sub','H5J'),('prod','6VP')):
    cplx = prot + rep.ligand if arm=='tsa' else P.swap_complex(prot, rep.ligand, resname, T)[0]
    pose = TC.build_pose(sess, cplx)
    lb = pose.max_n_blocks-1
    print(arm, 'ligand block:', pose.block_type(0, lb).name, 'n_atoms', pose.block_type(0, lb).n_atoms)
