import warnings, copy; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np
from maxi_eval import ligands as L, structures as S, geometry as G, tmol_core as TC, cst as C
T = L.build_templates(Path('ligands'))
rep = S.extract_replicates(Path('/mnt/work/projects/Maxi/7k4q.cif'),'7K4Q','tsa','crystal',T)[0]
site = G.detect_site(rep); aa = rep.protein + rep.ligand
sess = TC.Session.create(Path("ligands"), rep.protein, T); pose = TC.build_pose(sess, aa); bm = TC.block_map(pose)
cons = C.load_reference(Path('config/theozyme_reference.yaml'), 'hg3_lineage')
shift = {'geo.d_O_X':+0.6, 'geo.d_O_H':-0.5, 'geo.ang_O_H_X':-35, 'geo.ang_C_O_H':+30, 'geo.tors_syn':+90, 'geo.tors_C_O_H_X':-100, 'geo.tors_O_H_X_C3A':+150}
cons2 = [C.Constraint(c.key, c.target+shift.get(c.key,0), c.tol, c.kind, c.force_k, c.atoms) for c in cons if c.key in shift]
pose_c, applied, skipped = C.add_to_pose(pose, site, cons2, bm)
r = TC.energy_report(sess, pose_c, site['pocket_resids'])
prot, lig = TC.pose_to_replicate_arrays(pose_c)
m,_ = G.all_geometry(S.Replicate('x','7K4Q','tsa','tmol',prot,lig), site=site)
ev = C.evaluate(m, cons2)
print('tmol E_constraint = %.4f | numpy penalty = %.4f' % (r['E_constraint'], ev['cst.penalty']))
# relax with these violated constraints held: geometry should move toward shifted targets, E_cst drops
p1 = TC.relax(sess, pose_c, pocket_radius=8.0, repeats=1)
r1 = TC.energy_report(sess, p1, site['pocket_resids'])
prot, lig = TC.pose_to_replicate_arrays(p1)
m1,_ = G.all_geometry(S.Replicate('x','7K4Q','tsa','tmol',prot,lig), site=site)
ev1 = C.evaluate(m1, cons2)
print('after constrained relax: E_constraint %.4f  numpy penalty %.4f ; d_O_X %.2f -> %.2f (target %.2f)' % (r1['E_constraint'], ev1['cst.penalty'], m['geo.d_O_X'], m1['geo.d_O_X'], cons2[0].target))
from tmol.score import ScoreType
print('constraint weight after relax:', float(sess.sfxn.get_weight(ScoreType.constraint)), '(expected', sess.cst_weight, ')')
sess.sfxn.set_weight(ScoreType.constraint, sess.cst_weight)
r2 = TC.energy_report(sess, p1, site['pocket_resids'])
print('E_constraint after restoring weight: %.4f' % r2['E_constraint'])
