"""Run: $PY_TMOL tests/test_unit.py   (plain asserts; no pytest needed)"""
import math, warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np
from maxi_eval import geometry as G, ligands as L, structures as S, cst as C, field as F

# --- vector geometry
assert abs(G.angle([1,0,0],[0,0,0],[0,1,0]) - 90) < 1e-6
assert abs(abs(G.dihedral([1,0,0],[0,0,0],[0,0,1],[-1,0,1])) - 180) < 1e-6
assert abs(G.dihedral([1,0,0],[0,0,0],[0,0,1],[1,0,1])) < 1e-6
assert abs(G.wrap180(190) + 170) < 1e-9 and abs(G.wrap180(-190) - 170) < 1e-9
assert G.circ_std_deg([179, -179]) < 3            # wraps around +-180
assert G.circ_std_deg([0, 90, 180, 270]) > 60
# ideal H on a planar centre lies in the plane, 1 A away, opposite the neighbours
X = np.array([0.,0,0]); h = G.ideal_h(X, np.array([-1.,1,0]), np.array([-1.,-1,0]))
assert abs(np.linalg.norm(h)-1) < 1e-6 and h[0] > 0.9 and abs(h[2]) < 1e-9
# --- ligand renaming + superposition
T = L.build_templates(Path('ligands'))
assert set(T) == {'6NT','H5J','6VP'}
# every role atom named in ROLES / CORRESPONDENCE must exist in the CCD-derived template
for k, roles in L.ROLES.items():
    names = set(T[k].names)
    need = [roles[r] for r in ('abstraction','oxyanion','oxyanion_partner','ring_partner','ring_attached','nitro_N')] + roles['nitro_O'] + roles['benzo'] + ([roles['abstraction_H']] if roles['abstraction_H'] else [])
    assert set(need) <= names, (k, set(need)-names)
for k, corr in L.CORRESPONDENCE_TO_TSA.items():
    assert set(corr) <= set(T[k].names) and set(corr.values()) <= set(T['6NT'].names), k
# graph-matching rename survives shuffling for all three ligands, and H5J / 6VP (isomers, same C7N2O3) are NOT confused
for k, t in T.items():
    el = [a.GetSymbol() for a in t.mol.GetAtoms()]; perm = np.random.RandomState(3).permutation(len(el))
    nm = L.match_names(t, [el[i] for i in perm], t.coords()[perm])
    assert all(n is None if el[i]=='H' else n == t.names[i] for n, i in zip(nm, perm)), k
for a_, b_ in (('H5J','6VP'),('6VP','H5J'),('6NT','H5J')):
    el = [x.GetSymbol() for x in T[a_].mol.GetAtoms()]
    assert all(n is None for n in L.match_names(T[b_], el, T[a_].coords())), (a_, b_)
for k in ('H5J','6VP'):
    _, _, rmsd = L.place_onto_tsa(T[k], dict(zip(T['6NT'].names, T['6NT'].coords())))
    assert rmsd < 0.25, (k, rmsd)
# --- field: a +1 charge 2 A along +x gives E(0) pointing -x with magnitude 14.4/4 V/A
E = F.field_vector(np.zeros((1,3)), np.array([[2.,0,0]]), np.array([1.0]), 'vac')[0]
assert abs(E[0] + 14.3996/4) < 1e-3 and abs(E[1]) < 1e-9
# --- evaluate(): circular wrap, pass flags
cons = [C.Constraint('geo.tors_syn', 170.0, 30.0, 'circular'), C.Constraint('geo.d_O_X', 2.5, 0.3, 'linear')]
ev = C.evaluate({'geo.tors_syn': -175.0, 'geo.d_O_X': 3.2}, cons)
assert abs(ev['cst.tors_syn.dev'] - 15.0) < 1e-6 and ev['cst.tors_syn.pass'] == 1.0
assert ev['cst.d_O_X.pass'] == 0.0 and ev['cst.n_pass'] == 1.0 and ev['cst.penalty'] > 0
# --- enzdes .cst round trip on a real crystal copy (7K4Q): parser values vs geometry module
blocks = C.parse_enzdes_cst(Path('tests/example_theozyme.cst').read_text())
assert len(blocks) == 1 and blocks[0]['ligand']['atoms'] == ['N3','C3A','N2'] and blocks[0]['protein']['res3'] == ['ASP']
rep = S.extract_replicates(Path('/mnt/work/projects/Maxi/7k4q.cif'),'7K4Q','tsa','crystal',T)[0]
out = C.evaluate_enzdes_cst(rep, blocks)
site = G.detect_site(rep); m,_ = G.all_geometry(rep, site=site)
assert abs(out['cst1.distanceAB.value'] - m['geo.d_O_X']) < 1e-4          # same heavy-atom distance
assert out['cst1.resid'] == site['base']['resid']
print('all unit tests passed;  cst1 deviations on 7K4Q:', {k: round(v,1) for k,v in out.items() if k.endswith('.dev')})
