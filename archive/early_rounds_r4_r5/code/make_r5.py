"""Round 5 (channel arm; reuses the r4 base weights): writes rfd4/configs/inference_experiment/r4_<motif>_<arm>.yaml (steered arms) and runs/<campaign>/spec/*.
Weights: side-chain polar atoms near the base penalised 4x harder than r3; TMol weight halved; motif rigidity (base + oxyanion partner RMSD) rewarded;
pocket target 0.5 +- 0.1; nitro H-bond is CONDITIONED in the specs (C_HB acceptor), only 'donor near base' is penalised; search 4 branches x beam 3, checkpoints 30/42/52/60."""
import os, shutil
H = os.path.expanduser("~/HACKATHON_ETH"); D = H + "/rfd4/configs/inference_experiment"
M = {"hgwt": ("20261121_r4_hg317-tips_with-trp", "kemp4_hg317-tips_with-trp.yaml", "kemp_chem_hg"),
     "hgnt": ("20261121_r4_hg317-tips_no-trp", "kemp4_hg317-tips_no-trp.yaml", "kemp_chem_hg"),
     "gdm": ("20261121_r4_gdm1483-glu-ser", "kemp4_gdm1483-glu-ser.yaml", "kemp_chem")}
def rewards(chem, tm, channel=True):
    t = {"kemp_base.base_sc_polar_4p5": (-4.0, 3.0, 1.0), "kemp_base.base_cation_6": (-3.0, 2.0, 1.0), "kemp_base.base_sasa_rel": (-2.0, 0.5, 0.5),
         "kemp_base.base_apolar_frac": (0.5, 1.0, 1.0), "kemp_pocket3.volume_ratio_band_error": (-1.0, 0.6, 0.3), "kemp_pocket3.ligand_sasa_frac": (-0.5, 0.2, 0.2),
         "kemp_base.nitro_hb_near_base": (-1.0, 2.0, 2.0), "kemp_base.nitro_hb_far": (0.5, 1.0, 1.0),
         "kemp_base.lig_clash_atoms": (-1.0, 3.0, 3.0), "kemp_base.pp_clash_frac": (-1.0, 0.01, 0.01),
         f"{chem}.base_rmsd": (-2.0, 1.0, 0.25), f"{chem}.hb_rmsd": (-1.5, 1.0, 0.25)}
    t["kemp_pocket3.ligand_sasa_frac"] = (-0.25, 0.2, 0.2)       # relaxed: an open channel exposes some ligand surface
    if channel: t["kemp_channel.channel_band_error"] = (-2.0, 1.0, 0.3)   # bottleneck probe radius within 1.4-1.9 A (water passes, not wider)
    s = ", ".join(f"{k}: {{weight: {w}, clip: {{max: {c}}}, normalize_by: {n}}}" for k, (w, c, n) in t.items())
    if tm: s += f", tmol_lig_ddg.interface_ddg: {{weight: {-tm}, normalize_by: 20.0, clip: {{min: -80.0, max: 20.0}}}}"
    return s
for mk, (camp, spec, chem) in M.items():
    cd = f"{H}/runs/{camp}"; os.makedirs(f"{cd}/spec", exist_ok=True); os.makedirs(f"{cd}/logs", exist_ok=True); shutil.copy(f"{H}/inputs/{spec}", f"{cd}/spec/{spec}")
    for arm, tm in (("steer_base_channel", 0),):
        met = f"kemp_base, kemp_pocket3, kemp_channel, {chem}" + (", tmol_lig_ddg" if tm else ""); r = rewards(chem, tm)
        cps = "\n".join(f'    "{k}": {{{"use_clean_prediction: true, " if k != 60 else ""}metrics: [{met}], rewards: {{{r}}}}}' for k in (30, 42, 52, 60))
        txt = f"""# @package _global_
# {camp} / {arm}  (round 5, see runs/{camp}/README.md)
defaults:
  - inference_on_user_inputs
  - _self_
inputs: ${{oc.env:PROJECT_ROOT}}/../inputs/{spec}
out_dir: ${{oc.env:PROJECT_ROOT}}/../runs/{camp}/structures/{arm}
run_name_: {mk}_{arm}
search:
  n_branch: 4
  beam_width: 3
  checkpoints:
{cps}
"""
        open(f"{D}/r5_{mk}_{arm}.yaml", "w").write(txt); open(f"{cd}/spec/{arm}.experiment.yaml", "w").write(txt)
print("ok")
