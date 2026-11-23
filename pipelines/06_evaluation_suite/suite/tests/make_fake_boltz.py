"""Fabricate Boltz-shaped outputs (4BS0 + each of 6NT/H5J/6VP) to test the three-ligand flow. Values are NOT predictions."""
import sys, warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np
import biotite.structure.io.pdbx as pdbx
from maxi_eval import ligands as L, structures as S, protocol as P
root = Path(sys.argv[1]); T = L.build_templates()
rep = S.extract_replicates(Path('work/inputs/4BS0/tsa/4BS0_tsa_A1305.pdb'),'4BS0','tsa','crystal',T)[0]
prot = rep.protein[rep.protein.element!='H']
for name, code in (('4bs0_6NT','6NT'),('4bs0_H5J','H5J'),('4bs0_6VP','6VP')):
    cplx = prot + rep.ligand if code=='6NT' else P.swap_complex(prot, rep.ligand, code, T)[0]
    d = root/'boltz_results_fake'/'predictions'/name; d.mkdir(parents=True, exist_ok=True)
    f = pdbx.CIFFile(); pdbx.set_structure(f, cplx, data_block=name); f.write(str(d/f'{name}_model_0.cif'))
    print('wrote', d/f'{name}_model_0.cif', code)
