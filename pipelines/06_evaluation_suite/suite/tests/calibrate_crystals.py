"""Compute geometry for every Maxi crystal/NMR ligand copy; print table + dump JSON."""
import json, sys, warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
from maxi_eval import ligands as L, structures as S, geometry as G
T = L.build_templates(Path('ligands'))
rows = []
for f, src in [('29SB','nmr'),('4bs0','crystal'),('5RGA','crystal'),('7k4q','crystal'),('8RD5','crystal')]:
    reps = S.extract_replicates(Path('/mnt/work/projects/Maxi')/f'{f}.cif', f.upper(), 'tsa', src, T)
    for r in reps:
        m, site = G.all_geometry(r)
        m['rep']=r.rep_id; m['n_atoms']=len(r.protein)
        m['base']=f"{site['base']['resname']}{site['base']['resid']}" if site['base'] else None
        m['wall']=f"{site['wall']['resname']}{site['wall']['resid']}" if site['wall'] else None
        m['donors']=';'.join(f"{d['resname']}{d['resid']}:{d['atom']}@{d['d']:.2f}" for d in site['donors'])
        rows.append(m)
df = pd.DataFrame(rows).set_index('rep')
df.to_csv('work/crystal_geometry.csv')
pd.set_option('display.width',250); pd.set_option('display.max_columns',50); pd.set_option('display.max_rows',100)
cols = ['base','geo.d_O_H','geo.d_O_X','geo.ang_O_H_X','geo.ang_C_O_H','geo.tors_syn','geo.oop_H','geo.tors_C_O_H_X','geo.tors_O_H_X_C3A']
print(df[cols].round(2).to_string())
print(df[['wall','geo.wall_d_centroid','geo.wall_plane_angle','geo.lig_rel_sasa','geo.enclosure','geo.n_contacts','geo.contact_min','geo.base_sasa_free_frac']].round(2).to_string())
print(df[['donors','geo.n_donors_3p5','geo.donor_best_d','geo.donor_best_is_bb']].round(2).to_string())
print(df[[c for c in df.columns if 'wsite' in c or 'xtal_w' in c]].to_string())
