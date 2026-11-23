import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np
from maxi_eval import ligands as L, structures as S, geometry as G
T = L.build_templates(Path('ligands'))
# Independent check of syn/anti: which side of the C-O axis (in the carboxylate plane) is H on?
#  syn  = same side as the other carboxylate O  ;  anti = opposite side
for f in ['4bs0','7k4q','8RD5']:
    r = S.extract_replicates(Path('/mnt/work/projects/Maxi')/f'{f}.cif', f.upper(),'tsa','crystal',T)[0]
    m, site = G.all_geometry(r)
    b = site['base']; P = G.Atoms.of(r.protein)
    C,O,Oo = [P.xyz(b['resid'],b[k]) for k in ('C','O','O_other')]
    ly = G.lig_xyz(r); H,_ = G.abstraction_h(r, site['roles'], ly)
    u = G._u(O-C); p = G._u((Oo-C)-((Oo-C)@u)*u)       # in-plane unit vec perpendicular to C-O, towards O'
    HO = H-O; side = HO@p; along = HO@u
    print(f, 'proximal', b['O'], '| H-O component toward O\' side: %+.2f A  along C->O: %+.2f A  => %s' % (side, along, 'SYN' if side>0 else 'ANTI'),
          '| tors_syn %.0f' % m['geo.tors_syn'])
    # also N3 (heavy atom) -- ligand-independent of ideal-H assumption
    X = ly[site['roles']['abstraction']]; XO = X-O
    print('     heavy-atom N3 side: %+.2f  (same convention)' % (XO@p), ' angle C-O...N3 %.0f' % G.angle(C,O,X), ' dihedral O\'-C-O...N3 %.0f' % G.dihedral(Oo,C,O,X))
