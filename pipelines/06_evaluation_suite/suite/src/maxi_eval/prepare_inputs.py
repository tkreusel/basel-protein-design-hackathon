"""Stage 0: standardise every input into work/inputs/<entry>/<arm>/*.pdb + manifest.json.

Crystals / NMR: one standardized PDB per ligand copy / model (protein chain A, ligand chain B, canonical names, no H).
Boltz refolds : one per model_k; the PAE / confidence / pLDDT files are recorded for ipSAE2.
"""
from __future__ import annotations
import argparse, json, re, warnings
from pathlib import Path
import yaml
warnings.filterwarnings("ignore")
from . import ligands as L, structures as S, geometry as G


ARM_TOKENS = {"tsa": "tsa", "sub": "sub", "prod": "prod", "6NT": "tsa", "H5J": "sub", "6VP": "prod",
              "analogon": "tsa", "analog": "tsa", "analogue": "tsa", "substrate": "sub", "product": "prod", "apo": "apo"}


def infer_arm_from_cif(cif: Path) -> str | None:
    """Arm from the ligand residue name actually present in the file (6NT / H5J / 6VP)."""
    import biotite.structure.io.pdbx as pdbx
    st = pdbx.get_structure(pdbx.CIFFile.read(cif), model=1, include_bonds=False)
    for rn in dict.fromkeys(st.res_name):
        if str(rn) in L.RESNAME_TO_ARM:
            return L.RESNAME_TO_ARM[str(rn)]
    return None


def discover_boltz(root: Path, rx: re.Pattern, name_map: dict | None, known_entries: dict, suffix: str = "_boltz"):
    """yield (entry, arm, k, cif, pae, conf, plddt).

    Name -> (entry, arm):  name_map CSV > regex `<entry>(__|_)<arm|CCD>` > whole name as entry.  The arm is always cross-checked
    against the ligand found in the CIF (that is what was actually predicted). `entry` is matched case-insensitively to the
    configured entries (7k4q -> 7K4Q) and gets `suffix` so crystal and refold results stay separate.
    """
    low = {k.lower(): k for k in known_entries}
    for cif in sorted(root.rglob("*_model_*.cif")):
        m = re.match(r"^(?P<name>.+)_model_(?P<k>\d+)\.cif$", cif.name)
        if not m:
            continue
        name, k = m["name"], int(m["k"])
        arm = None
        if name_map and name in name_map:
            entry, arm = name_map[name]
        else:
            mm = rx.match(name)
            entry, tok = (mm["entry"], mm["arm"]) if mm else (name, None)
            arm = ARM_TOKENS.get(tok) if tok else None
            entry = low.get(entry.lower(), entry)
            if not entry.endswith(suffix):
                entry += suffix
        found = infer_arm_from_cif(cif)
        if found is None and arm == "apo":
            found = "apo"            # protein-only prediction: no ligand to infer from, the name token decides
        if found is None:
            print(f"  [skip] {cif.name}: no 6NT/H5J/6VP ligand in the file")
            continue
        if arm and arm != found:
            print(f"  [warn] {cif.name}: name says arm '{arm}' but the ligand is {L.ARM_TO_RESNAME.get(found)} ({found}); using '{found}'")
        d = cif.parent
        pae, conf, pl = d / f"pae_{name}_model_{k}.npz", d / f"confidence_{name}_model_{k}.json", d / f"plddt_{name}_model_{k}.npz"
        yield entry, found, k, cif, (pae if pae.exists() else None), (conf if conf.exists() else None), (pl if pl.exists() else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/entries.yaml")
    ap.add_argument("--out", default="work/inputs")
    ap.add_argument("--ligand-dir", default="ligands")
    ap.add_argument("--boltz-root", default=None, help="overrides config")
    ap.add_argument("--name-map", default=None)
    a = ap.parse_args()
    cfg = yaml.safe_load(Path(a.config).read_text())
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    T = L.build_templates(Path(a.ligand_dir))
    man = {"entries": {}}

    def add(entry, rep, pdbdir, extra=None):
        pdbdir.mkdir(parents=True, exist_ok=True)
        pdb = pdbdir / (rep.rep_id.replace(":", "_") + ".pdb")
        S.write_complex_pdb(rep, pdb)
        rec = dict(rep_id=rep.rep_id, arm=rep.arm, source=rep.source, pdb=str(pdb.resolve()), src=rep.meta["file"],
                   model=rep.meta["model"], n_res=int(len(set(rep.protein.res_id))))
        if rep.arm != "apo":   # catalytic site from geometry (see geometry.detect_site) for every ligand arm; designs can override in entries.yaml
            site = G.detect_site(rep, man["entries"][entry].get("catalytic_overrides"))
            rec["site"] = {"base": site["base"], "donors": site["donors"][:4], "wall": ({k: v for k, v in site["wall"].items() if k != "ring"} if site["wall"] else None),
                           "pocket_resids": site["pocket_resids"]}
        rec.update(extra or {})
        man["entries"][entry]["replicates"].append(rec)

    all_entries = dict(cfg["entries"])
    for k, e in (cfg.get("extra_entries") or {}).items():
        if e.get("enabled"):
            all_entries[k] = e
    for entry, e in all_entries.items():
        man["entries"][entry] = dict(kind=e["kind"], variant=e.get("variant"), kcat_KM=e.get("kcat_KM"), replicates=[], catalytic_overrides=e.get("catalytic_overrides"))
        reps = S.extract_replicates(Path(e.get("dir", cfg["reference_dir"])) / e["file"], entry, "tsa", e["kind"], T)
        if e.get("max_replicates"):
            reps = reps[: e["max_replicates"]]
        print(f"{entry}: {len(reps)} TSA replicate(s) from {e['file']}")
        for r in reps:
            add(entry, r, out / entry / "tsa")

    root = a.boltz_root or (cfg.get("boltz_refolds") or {}).get("root")
    if root:
        rx = re.compile(cfg["boltz_refolds"]["name_regex"])
        nm = None
        if a.name_map:
            import csv
            nm = {r["name"]: (r["entry"], r["arm"]) for r in csv.DictReader(open(a.name_map))}
        n = 0
        used_keys: set = set()
        for entry, arm, k, cif, pae, conf, pl in discover_boltz(Path(root), rx, nm, all_entries, cfg["boltz_refolds"].get("entry_suffix", "_boltz")):
            from .labels import load_labels, label_for
            man["entries"].setdefault(entry, dict(kind="boltz", variant=None, kcat_KM=(label_for(entry, load_labels(a.config)) or {}).get("kcat_KM"),
                                                  replicates=[], catalytic_overrides=None))
            reps = [S.extract_apo(cif, entry)] if arm == "apo" else S.extract_replicates(cif, entry, arm, "boltz", T)
            if not reps:
                print(f"  [warn] {cif.name}: no {arm} ligand parsed (wrong ligand? atoms?)"); continue
            k_eff = k
            while (entry, arm, k_eff) in used_keys:       # same protein+arm predicted again (e.g. a second batch): keep BOTH as replicates
                k_eff += 1
            if k_eff != k:
                print(f"  [note] {entry} {arm}: model {k} already registered -> registered as replicate m{k_eff} ({cif.parent.parent.parent.name})")
            used_keys.add((entry, arm, k_eff))
            r = reps[0]; r.rep_id = f"{entry}:{arm}:boltz_m{k_eff}"; r.meta["model"] = k_eff; k = k_eff
            add(entry, r, out / entry / arm, dict(boltz_cif=str(cif), boltz_pae=str(pae) if pae else None,
                                                  boltz_conf=str(conf) if conf else None, boltz_plddt=str(pl) if pl else None, boltz_model=k))
            n += 1
        print(f"Boltz refolds registered: {n}")
    (out / "manifest.json").write_text(json.dumps(man, indent=1))
    print("manifest:", out / "manifest.json")

if __name__ == "__main__":
    main()
