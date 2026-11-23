import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, torch
from maxi_eval import ligands as L, structures as S, geometry as G, tmol_core as TC, cst as C
T = L.build_templates(Path('ligands'))
rep = S.extract_replicates(Path('/mnt/work/projects/Maxi/7k4q.cif'),'7K4Q','tsa','crystal',T)[0]
site = G.detect_site(rep); aa = rep.protein + rep.ligand
sess = TC.Session.create(Path("ligands"), rep.protein, T)
pose = TC.build_pose(sess, aa)
bm = TC.block_map(pose)
cons = C.load_reference(Path('config/theozyme_reference.yaml'), 'hg3_lineage')
# make the constraints bite: target the Gln50-bearing group's tighter geometry
pose_c, applied, skipped = C.add_to_pose(pose, site, cons, bm)
print('applied', applied, 'skipped', skipped)
r = TC.energy_report(sess, pose_c, site['pocket_resids'])
print('E_constraint (tmol)  = %.4f' % r['E_constraint'])
# numpy expectation on the same coordinates (tmol H)
prot, lig = TC.pose_to_replicate_arrays(pose_c)
rep2 = S.Replicate('x','7K4Q','tsa','tmol',prot,lig)
m,_ = G.all_geometry(rep2, site=site)
ev = C.evaluate(m, [c for c in cons if c.key in applied])
print('numpy flat-penalty   = %.4f' % ev['cst.penalty'], '  (per-constraint z):', {k:round(v,2) for k,v in ev.items() if k.endswith('.z')})
# does the constraint set survive relax, and does relaxation reduce it?
p1 = TC.relax(sess, pose_c, pocket_radius=8.0, repeats=1)
print('constraint_set after relax:', p1.constraint_set is not None)
r1 = TC.energy_report(sess, p1, site['pocket_resids'])
print('E_constraint after constrained pocket relax = %.4f' % r1['E_constraint'])
# CA restraints
p2 = TC.add_ca_restraints(pose)
print('CA-restrained pose E_constraint:', TC.energy_report(sess, p2)['E_constraint'])
