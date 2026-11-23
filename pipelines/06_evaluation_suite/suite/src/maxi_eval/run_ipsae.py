"""ipSAE2 (protein<->ligand) + site-focused Boltz confidence for every Boltz refold in the manifest.

For each Boltz model:
  ipsae.*     directional protein->ligand ipSAE from ipSAE2 (`ipsae` = protein-aligned, no PAE cut-off on ligand, d0 = 2 A)
  boltz.*     ipTM / ligand_iptm / complex_iplddt from Boltz's own confidence JSON
  site.*      mean PAE between the catalytic residues and the ligand (both directions) and pLDDT of the catalytic
              residues / ligand atoms -> 'does the model trust THIS active-site geometry', not the whole interface
Per (entry, arm): across k models -> mean and sd (seed/sample consistency), and ligand-pose RMSD spread.
Arms tsa / sub / prod are scored separately; the summary stage forms ipsae_tsa - ipsae_sub.

Runs in the RFD4/analysis venv (numpy, pandas); calls the `ipsae` executable from the ipsae2 env.
"""
from __future__ import annotations
import argparse, json, os, subprocess, tempfile
from pathlib import Path
import numpy as np, pandas as pd


def run_ipsae(cif, pae, conf, out_tsv, ipsae_bin):
    cmd = [ipsae_bin, str(cif), str(pae)] + ([str(conf)] if conf else []) + ["--output-file", str(out_tsv)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"ipsae failed: {r.stderr[-400:]}")
    return pd.read_csv(out_tsv, sep="\t")


def prot_lig_row(df: pd.DataFrame) -> pd.Series:
    sel = df[(df.type1 == "prot") & (df.type2 == "lig")]
    if sel.empty:
        sel = df[(df.type1 == "lig") & (df.type2 == "prot")]
        if sel.empty:
            raise RuntimeError("no protein-ligand row in ipSAE2 output (ligand beyond 5 A contact? use --all-pairs)")
    return sel.iloc[0]


def site_confidence(pae_npz, plddt_npz, n_res, n_lig, cat_resids: list[int], lig_names: list[str] | None, abstraction: str | None) -> dict:
    out = {}
    z = np.load(pae_npz)
    pae = z[z.files[0]]
    tok_prot = np.array([r - 1 for r in cat_resids if 0 < r <= n_res])
    lig_tok = np.arange(n_res, n_res + n_lig)
    if len(tok_prot) and pae.shape[0] >= n_res + n_lig:
        out["site.pae_cat_to_lig"] = float(pae[np.ix_(tok_prot, lig_tok)].mean())   # error of ligand when aligned on catalytic residues
        out["site.pae_lig_to_cat"] = float(pae[np.ix_(lig_tok, tok_prot)].mean())   # error of catalytic residues when aligned on ligand
        if lig_names and abstraction in lig_names:
            li = n_res + lig_names.index(abstraction)
            out["site.pae_cat_to_abstr"] = float(pae[tok_prot, li].mean())
    if plddt_npz and Path(plddt_npz).exists():
        p = np.load(plddt_npz); pl = p[p.files[0]]
        out["site.plddt_cat"] = float(pl[tok_prot].mean()) if len(tok_prot) else float("nan")
        out["site.plddt_lig"] = float(pl[lig_tok].mean()) if len(lig_tok) and pl.shape[0] >= n_res + n_lig else float("nan")
        out["site.plddt_protein"] = float(pl[:n_res].mean())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="work/inputs")
    ap.add_argument("--out", default="work/ipsae")
    ap.add_argument("--ipsae-bin", default=os.environ.get("IPSAE_BIN", "ipsae"))
    a = ap.parse_args()
    man = json.loads((Path(a.inputs) / "manifest.json").read_text())
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    # catalytic residues come from the entry's TSA replicate(s) (same sequence => same numbering)
    for entry, e in man["entries"].items():
        tsa = [r for r in e["replicates"] if r["source"] == "boltz" and r.get("site") and r["site"].get("base")]
        tsa.sort(key=lambda r: r["arm"] != "tsa")   # prefer the TSA arm, else any arm with a detected base
        cat = []
        if tsa and tsa[0].get("site"):
            s = tsa[0]["site"]
            cat = [x["resid"] for x in (s["base"], s["wall"]) if x] + [d["resid"] for d in s["donors"][:1] if not d["backbone"]]
        for r in e["replicates"]:
            if r["source"] != "boltz" or not r.get("boltz_pae") or r["arm"] == "apo":
                continue
            tag = r["rep_id"].replace(":", "_")
            try:
                df = run_ipsae(r["boltz_cif"], r["boltz_pae"], r.get("boltz_conf"), out / f"{tag}.tsv", a.ipsae_bin)
                row = prot_lig_row(df)
                rec = {"rep_id": r["rep_id"], "entry": entry, "arm": r["arm"], "model": r["boltz_model"], "status": "ok",
                       "ipsae.ipsae": float(row.ipsae), "ipsae.p2l": float(row.ipsae12), "ipsae.l2p": float(row.ipsae21),
                       "ipsae.iPAE12": float(row.iPAE12), "ipsae.iPAE21": float(row.iPAE21),
                       "ipsae.n_scored_contact12": float(row.n_scored_contact12), "ipsae.n_scored_contact21": float(row.n_scored_contact21)}
                if r.get("boltz_conf"):
                    c = json.load(open(r["boltz_conf"]))
                    rec.update({f"boltz.{k}": c[k] for k in ("confidence_score", "ptm", "iptm", "ligand_iptm", "complex_plddt", "complex_iplddt", "complex_ipde") if k in c})
                n_res = r["n_res"]
                import biotite.structure.io.pdb as pdb
                lig = pdb.PDBFile.read(r["pdb"]).get_structure(model=1)
                lig = lig[(lig.chain_id == "B") & (lig.element != "H")]
                rec.update(site_confidence(r["boltz_pae"], r.get("boltz_plddt"), n_res, len(lig), cat, list(lig.atom_name), None))
                rows.append(rec)
            except Exception as exc:  # fail loud, continue
                rows.append({"rep_id": r["rep_id"], "entry": entry, "arm": r["arm"], "status": "failed", "error": f"{type(exc).__name__}: {exc}"})
                print("FAILED", r["rep_id"], exc)
    df = pd.DataFrame(rows)
    df.to_csv(out / "ipsae_models.csv", index=False)
    if not df.empty and (df.status == "ok").any():
        ok = df[df.status == "ok"]
        num = ok.select_dtypes("number").columns.drop("model", errors="ignore")
        g = ok.groupby(["entry", "arm"])[list(num)].agg(["mean", "std", "count"])
        g.columns = [f"{a}|{b}" for a, b in g.columns]
        g.reset_index().to_csv(out / "ipsae_by_entry_arm.csv", index=False)
    print(f"{len(rows)} Boltz models processed -> {out}")

if __name__ == "__main__":
    main()
