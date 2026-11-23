import sys, time, warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, torch
from maxi_eval import ligands as L, structures as S, geometry as G, tmol_core as TC
T = L.build_templates(Path('ligands'))
rep = S.extract_replicates(Path('/mnt/work/projects/Maxi/7k4q.cif'),'7K4Q','tsa','crystal',T)[0]
site = G.detect_site(rep)
aa = rep.protein + rep.ligand
t=time.time(); sess = TC.Session.create(Path("ligands"), rep.protein, T); print('session', round(time.time()-t,1),'s device', sess.device)
pose = TC.build_pose(sess, aa)
r0 = TC.energy_report(sess, pose, site['pocket_resids'])
print('RAW  E_total %.1f  per_res %.3f  int.total %.2f  int.hbond %.2f rep_intra %.1f rep_inter %.1f' % (r0['E_total'], r0['E_total_per_res'], r0['int.total'], r0.get('int.hbond',0), r0['rep_intra_total'], r0['rep_inter_total']))
mode = sys.argv[1] if len(sys.argv)>1 else 'pocket'
t=time.time()
if mode=='pocket':
    p1 = TC.relax(sess, pose, pocket_radius=8.0, move_ligand=True, repeats=1)
else:
    p1 = TC.relax(sess, TC.add_ca_restraints(pose), repeats=1)
print('relax', mode, round(time.time()-t,1),'s')
r1 = TC.energy_report(sess, p1, site['pocket_resids'])
print('RELAXED E_total %.1f (cst %.1f) per_res %.3f  int.total %.2f  int.hbond %.2f rep_intra %.1f rep_inter %.1f' % (r1['E_total'], r1['E_constraint'], r1['E_total_per_res'], r1['int.total'], r1.get('int.hbond',0), r1['rep_intra_total'], r1['rep_inter_total']))
prot, lig = TC.pose_to_replicate_arrays(p1)
rep2 = S.Replicate('x','7K4Q','tsa','tmol',prot,lig)
m, _ = G.all_geometry(rep2, site={**site})
print({k: round(v,2) for k,v in m.items() if k in ('geo.d_O_H','geo.d_O_X','geo.ang_O_H_X','geo.tors_syn','geo.h_source_model','geo.donor_best_d')})
