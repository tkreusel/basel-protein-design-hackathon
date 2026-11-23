"""Writes rfd4/configs/inference_experiment/kemp3_*.yaml (steered arms, round 3). Run with any python3."""
import os
D = os.path.expanduser("~/HACKATHON_ETH/rfd4/configs/inference_experiment")
def terms(sc, cat, sasa, apol, pocket, nfar, nnear, lsasa, clash, tm):
    t = {"kemp_base.base_sc_polar_4p5": (-sc, 4.0, 2.0), "kemp_base.base_cation_6": (-cat, 2.0, 1.0), "kemp_base.base_sasa_rel": (-sasa, 0.5, 0.5),
         "kemp_base.base_apolar_frac": (apol, 1.0, 1.0), "kemp_pocket2.volume_ratio_band_error": (-pocket, 0.6, 0.3),
         "kemp_base.nitro_hb_far": (nfar, 1.0, 1.0), "kemp_base.nitro_hb_near_base": (-nnear, 2.0, 2.0),
         "kemp_pocket2.ligand_sasa_frac": (-lsasa, 0.2, 0.2), "kemp_base.lig_clash_atoms": (-clash, 3.0, 3.0), "kemp_base.pp_clash_frac": (-clash, 0.01, 0.01)}
    s = ", ".join(f"{k}: {{weight: {w}, clip: {{max: {c}}}, normalize_by: {n}}}" for k, (w, c, n) in t.items() if w != 0)
    if tm: s += f", tmol_lig_ddg.interface_ddg: {{weight: {-tm}, normalize_by: 20.0, clip: {{min: -80.0, max: 20.0}}}}"
    return s
W1 = dict(sc=1.0, cat=2.0, sasa=1.0, apol=0.5, pocket=1.0, nfar=1.0, nnear=1.0, lsasa=1.0, clash=1.0)          # balanced
W2 = dict(sc=1.5, cat=3.0, sasa=2.0, apol=1.0, pocket=0.5, nfar=0.5, nnear=1.0, lsasa=0.5, clash=1.0)          # base-heavy
arms = {"w1": (W1, 0, "kemp3_rfd4_buried.yaml"), "w2": (W2, 0, "kemp3_rfd4_buried.yaml"), "w1t": (W1, 1.0, "kemp3_rfd4_buried.yaml"), "w2t": (W2, 1.0, "kemp3_rfd4_buried.yaml"),
        "hg_w1": (W1, 0, "kemp3_rfd4_hgnotrp.yaml"), "hg_w1t": (W1, 1.0, "kemp3_rfd4_hgnotrp.yaml")}
for name, (w, tm, spec) in arms.items():
    r = terms(**w, tm=tm); m = "kemp_base, kemp_pocket2" + (", tmol_lig_ddg" if tm else "")
    cps = "\n".join(f'    "{k}": {{{"use_clean_prediction: true, " if k != 60 else ""}metrics: [{m}], rewards: {{{r}}}}}' for k in (20, 35, 50, 60))
    open(f"{D}/kemp3_{name}.yaml", "w").write(f"""# @package _global_
# Round 3 steered arm "{name}": base desolvation/no polar-or-cation near base + crystal-sized pocket + low ligand SASA + nitro H-bond (donor far from base) + no clashes{' + TMol dG' if tm else ''}.
defaults:
  - inference_on_user_inputs
  - _self_
inputs: ${{oc.env:PROJECT_ROOT}}/../rfd4_runs/inputs/{spec}
out_dir: ${{oc.env:PROJECT_ROOT}}/../rfd4_runs/out/20261121_k3_{name}
run_name_: k3_{name}
search:
  n_branch: 3
  beam_width: 2
  checkpoints:
{cps}
""")
print("written", list(arms))
