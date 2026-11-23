# Quick check of RFD3 outputs: length, loop fraction (biotite P-SEA), ligand SASA fraction, residues within 4.5 A of ligand.
import sys,glob,gzip,numpy as np,biotite.structure as s,biotite.structure.io.pdbx as pdbx
from biotite.structure import sasa
for f in sorted(glob.glob(sys.argv[1]+"/*.cif.gz")):
    with gzip.open(f,"rt") as h: a=pdbx.get_structure(pdbx.CIFFile.read(h),model=1,use_author_fields=False)
    a=a[a.element!="H"]; lig=a.hetero; p=a[~lig]; l=a[lig]
    ca=p[p.atom_name=="CA"]; sse=s.annotate_sse(p) if len(ca)>5 else []
    loop=np.mean(np.array(sse)=="c") if len(sse) else np.nan
    sa=sasa(a,vdw_radii="ProtOr",point_number=200); fr=np.nanmean(sa[lig]/np.maximum(sasa(l,vdw_radii="ProtOr",point_number=200),1e-6)) if lig.any() else np.nan
    cell=s.CellList(p,5.0); near=set(p.res_id[np.unique(np.concatenate([cell.get_atoms(c,4.5) for c in l.coord]))]) if lig.any() else set()
    print(f.split("/")[-1][:40],"L=%d loop=%.2f lig_SASA_frac=%.2f n_res_4.5A=%d"%(len(ca),loop,fr,len(near)))
