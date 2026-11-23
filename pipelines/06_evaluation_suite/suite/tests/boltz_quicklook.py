import json, warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import pandas as pd, numpy as np
from maxi_eval import ligands as L, structures as S, geometry as G, cst as C
pd.set_option('display.width',220); pd.set_option('display.max_columns',40)
man = json.load(open('work/inputs/manifest.json')); T = L.build_templates(Path('ligands'))
cons = C.load_reference(Path('config/theozyme_reference.yaml'),'hg3_lineage')
ip = pd.read_csv('work/ipsae/ipsae_models.csv')
print(ip[['rep_id','status','ipsae.p2l','ipsae.l2p','boltz.iptm','boltz.ligand_iptm','boltz.complex_iplddt','site.plddt_cat','site.plddt_lig','site.plddt_protein','site.pae_cat_to_lig','site.pae_lig_to_cat']].round(2).to_string(index=False))
rows=[]
for entry,e in man['entries'].items():
    if not entry.endswith('_boltz'): continue
    r = e['replicates'][0]
    rep = S.extract_replicates(Path(r['pdb']),entry,'tsa','boltz',T)[0]
    m, site = G.all_geometry(rep)
    ev = C.evaluate(m, cons)
    rows.append(dict(entry=entry, base=(f"{site['base']['resname']}{site['base']['resid']}" if site['base'] else None),
        d_base_X=round(m.get('geo.d_O_X',np.nan),2), ang_OHX=round(m.get('geo.ang_O_H_X',np.nan),0), tors_syn=round(m.get('geo.tors_syn',np.nan),0),
        donor=(f"{site['donors'][0]['resname']}{site['donors'][0]['resid']}@{site['donors'][0]['d']:.2f}" if site['donors'] else '-'),
        wall=(f"{site['wall']['resname']}{site['wall']['resid']}" if site['wall'] else '-'),
        rel_sasa=round(m['geo.lig_rel_sasa'],2), enclosure=round(m['geo.enclosure'],2), contact_min=round(m['geo.contact_min'],2),
        cst_pass=f"{int(ev['cst.n_pass'])}/{int(ev['cst.n_total'])}"))
print(pd.DataFrame(rows).to_string(index=False))
