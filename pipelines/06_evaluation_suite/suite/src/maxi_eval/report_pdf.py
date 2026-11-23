"""PDF report builder (reportlab) used by make_deliverables."""
from __future__ import annotations
import datetime
from pathlib import Path
from xml.sax.saxutils import escape
import numpy as np, pandas as pd


def build_pdf(path: Path, d: dict):
    import matplotlib
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.utils import ImageReader
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether, PageBreak, PageTemplate, Paragraph, Spacer, Table, TableStyle)

    fd = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
    pdfmetrics.registerFont(TTFont("DV", str(fd / "DejaVuSans.ttf"))); pdfmetrics.registerFont(TTFont("DVB", str(fd / "DejaVuSans-Bold.ttf")))
    pdfmetrics.registerFont(TTFont("DVI", str(fd / "DejaVuSans-Oblique.ttf")))
    pdfmetrics.registerFontFamily("DV", normal="DV", bold="DVB", italic="DVI", boldItalic="DVB")

    W, H = landscape(A4); MG = 28; FW = W - 2 * MG
    INK, MUTED, ACC = colors.HexColor("#1f2933"), colors.HexColor("#5f6b7a"), colors.HexColor("#0072B2")
    st = dict(
        title=ParagraphStyle("t", fontName="DVB", fontSize=20, leading=24, textColor=INK, spaceAfter=4),
        sub=ParagraphStyle("s", fontName="DV", fontSize=9, leading=12, textColor=MUTED, spaceAfter=8),
        h1=ParagraphStyle("h1", fontName="DVB", fontSize=13, leading=16, textColor=ACC, spaceBefore=8, spaceAfter=5),
        h2=ParagraphStyle("h2", fontName="DVB", fontSize=10, leading=13, textColor=INK, spaceBefore=6, spaceAfter=3),
        body=ParagraphStyle("b", fontName="DV", fontSize=8.6, leading=11.8, textColor=INK, spaceAfter=3),
        bul=ParagraphStyle("bl", fontName="DV", fontSize=8.6, leading=11.8, textColor=INK, leftIndent=11, bulletIndent=0, spaceAfter=2.5),
        small=ParagraphStyle("sm", fontName="DV", fontSize=7, leading=9, textColor=MUTED),
        cell=ParagraphStyle("c", fontName="DV", fontSize=6.6, leading=8),
        cellb=ParagraphStyle("cb", fontName="DVB", fontSize=6.6, leading=8, textColor=colors.white),
        tiny=ParagraphStyle("ti", fontName="DV", fontSize=5.6, leading=6.9, textColor=INK),
    )
    P = lambda t, s="body": Paragraph(t, st[s])
    E = lambda t: escape(str(t))
    story: list = []

    R, N, NH, key, M, lab, retest = d["R"], d["N"], d["NH"], d["key"], d["M"], d["lab"], d["retest"]
    TOK, TYPE, TOOL = d["glossary"]
    Ri = R.set_index("candidate")
    g = lambda c, col: (Ri.loc[c, col] if c in Ri.index and col in Ri.columns else np.nan)
    fam_n = lab.family.value_counts().to_dict()
    sel = lambda grp, df=N: float(df.set_index("group").loc[grp, "p_selection_aware"])
    f2 = lambda v: "n/a" if pd.isna(v) else f"{v:+.2f}"
    LB, FL, IT, ET = "tmol swap|int.lk_ball_iso|tsa-prod", "tmol swap|int.fa_lk|tsa-sub", "tmol swap|int.total|tsa-sub", "tmol swap|E_total|tsa-sub"
    PAE, IPTM, PR = "bz.sub-prod.cat_to_pocket_pae", "ipSAE2/Boltz|boltz.ligand_iptm|sub-prod", "pr.tsa_pk.int_sum"

    # ---------------------------------------------------------------- 1. title + summary
    story += [P("Kemp-eliminase in-silico evaluation: which metrics track activity?", "title"),
              P(f"Report generated {datetime.date.today().isoformat()} · {d['nlab']} proteins with literature kcat/KM ({fam_n.get('HG3', 0)} HG3 lineage, {fam_n.get('1A53', 0)} 1A53, {fam_n.get('KABLE', 0)} KABLE) · "
                f"Boltz-2 refolds with TSA (6NT), substrate (H5J), product (6VP) and apo · tools: ipSAE2, PLACER, tmol, PyRosetta · {len(R)} candidate metrics screened", "sub")]
    story += [P("1. Summary", "h1")]
    bullets = [
        f"<b>Strongest replicated signal: tmol ligand–pocket solvation/electrostatic energies, as differences between the TSA and the substrate/product complexes.</b> "
        f"Best candidate: lk_ball_iso, TSA − product, ρ = {f2(g(LB,'rho_all19_tmolmean'))} over all 19 proteins (within-family ρ {f2(g(LB,'rho_within_family'))}, leave-one-protein-out {f2(g(LB,'rho_loo_min'))} … {f2(g(LB,'rho_loo_max'))}); "
        f"fa_lk TSA − substrate ρ = {f2(g(FL,'rho_all19_tmolmean'))}; total ligand–pocket interaction TSA − substrate ρ = {f2(g(IT,'rho_all19_tmolmean'))}. "
        f"More favourable TSA-than-product/substrate interaction goes with higher kcat/KM.",
        f"<b>Within the HG3 lineage (n = 12) the same quantity gives ρ = {f2(g(LB,'rho_HG3_only'))}</b> (selection-aware p = {sel('ALL', NH):.3f}; tmol set p = {sel('set:ARM-tmol', NH):.3f}). "
        f"Across all 19 proteins the whole-pool selection-aware p is {sel('ALL'):.2f} (tmol arm set p = {sel('set:ARM-tmol'):.3f}), i.e. suggestive, not conclusive.",
        f"<b>It replicates:</b> two independent tmol relax runs gave ρ = {f2(g(LB,'rho_run1'))} and {f2(g(LB,'rho_run2'))} for this quantity (test–retest ICC 0.91). "
        f"The earlier headline (total energy TSA − substrate, ρ = −0.82 in run 1) did <b>not</b> reproduce: run 2 gave {f2(g(ET,'rho_run2'))}; it was relax noise plus selection. Reported values here are the mean of the two runs.",
        f"<b>Boltz confidence differences:</b> catalytic→pocket PAE, substrate − product ρ = {f2(g(PAE,'rho_all19_tmolmean'))} (Boltz-PAE set p = {sel('set:EXTRA-bz'):.3f}); ligand ipTM, substrate − product ρ = {f2(g(IPTM,'rho_all19_tmolmean'))} "
        f"(ipSAE2 set p = {sel('set:ARM-ipSAE2'):.3f}). These are deterministic given the structure but Boltz-seed variance is untested; ipTM drops to ρ ≈ 0.4 within HG3.",
        f"<b>No support for the classic preorganisation metrics.</b> PLACER apo 'pre-formed' fraction ρ ≈ +0.35, apo–holo pocket motion ≈ −0.1, tmol strain energy wrong-signed; PLACER (p = {sel('set:ARM-PLACER'):.2f}) and Rosetta (p = {sel('set:ARM-Rosetta'):.2f}) carry no significant signal. "
        f"Single-arm absolute values are weaker than arm differences throughout.",
        "<b>Reading guide.</b> ρ is Spearman rank correlation with log10 kcat/KM. Negative ρ for an energy difference means 'lower (more favourable) energy ↔ more active'. Every candidate comes from a screen of "
        f"{len(R)} metrics on only {d['nlab']} proteins, so p-values are corrected for that selection (permutation null of the best |ρ|); see section 3 and the caveats in section 8.",
    ]
    story += [Paragraph(b, st["bul"], bulletText="•") for b in bullets]
    story += [P("Files delivered alongside this PDF (in the same folder)", "h2"),
              Paragraph("<b>kemp_eval_dataset_key.csv</b>: one row per labelled protein; kcat/KM, family, mutation count, the ~125 most important metrics (top-40 correlates, preorganisation, all ipSAE2/Boltz values, run-1/run-2 values of the headline tmol terms).", st["bul"], bulletText="•"),
              Paragraph("<b>kemp_eval_dataset_full.csv</b>: all 19 proteins, every master column (~1,200), 545 extra features (prefix x:), all screened candidates (prefix cand:), iPAE ligand→protein.", st["bul"], bulletText="•"),
              Paragraph("<b>kemp_eval_correlations_all.csv</b>: every candidate with ρ, p, BH q, run-1/run-2 ρ, within-family ρ, leave-one-out range, AUROC, old-label ρ, HG3-only ρ and a plain-language description.", st["bul"], bulletText="•"),
              Paragraph("<b>kemp_eval_data_dictionary.csv</b>: the glossary of section 5.", st["bul"], bulletText="•"), PageBreak()]

    # ---------------------------------------------------------------- 2. what was measured
    story += [P("2. What was measured", "h1"),
              P("Every protein was refolded with Boltz-2 in complex with each of three ligands and without ligand: the <b>transition-state analogue 6NT (TSA)</b>, the <b>substrate 5-nitrobenzisoxazole H5J</b> and the <b>product 6VP</b> (neutral phenol; "
                "CCD codes, so tmol/Rosetta use the real chemistry). The catalytic site (base carboxylate, oxyanion donors, wall residues) is detected from the TSA geometry. Quantities are computed per arm and then as <b>differences between arms</b> (TSA − substrate, TSA − product, substrate − product)."),
              P("Two ways of comparing ligands are used. <b>swap</b> (same-pocket): all three ligands are placed on the same TSA-relaxed pocket by superposition, so differences reflect the ligand only. <b>native</b>: each ligand in its own independent Boltz refold, so differences "
                "also include pocket changes and prediction noise. The signal found is mostly in the swap comparison; where both exist the interaction-energy parts agree in sign."),
              P("Tools: <b>ipSAE2</b> (protein→ligand interface confidence from Boltz PAE; also ipTM/ligand ipTM and site PAE/pLDDT), <b>PLACER</b> (ensemble of 20 sampled all-atom models per structure; geometry/near-attack fractions, ligand-pose spread, pocket-crop only), "
                "<b>tmol</b> (Rosetta beta2016 energy in PyTorch: relax, per-term ligand–pocket energies, strain, burial, void, buried unsatisfied polars, electric-field projections), <b>PyRosetta</b> (dG_separated, shape complementarity, packstat, H-bonds, BUNS) and "
                "extra features mined from raw outputs (PLACER ensemble distributions, tmol per-residue energies, Boltz pocket PAE/pLDDT, pocket composition, mutation count vs HG3)."),
              P("Did we compute ipSAE / iPAE?", "h2"),
              P("<b>Yes.</b> From ipSAE2/Boltz, per arm (TSA, substrate, product) and as the three differences: ipSAE protein→ligand (p2l) and ligand→protein (l2p), interface PAE protein→ligand (iPAE12) and — added in this report — ligand→protein (iPAE21), "
                "Boltz ipTM, ligand ipTM, complex interface pLDDT, and site PAE (catalytic residues ↔ ligand) and pLDDT. <b>Not computed:</b> pDockQ/pDockQ2/LIS (ipSAE2 does not return them for ligand tokens). If 'ipde' meant a different quantity, tell me which."),
              ]
    ipc = R[R.candidate.str.startswith(("ipSAE2/Boltz|", "ip.", "ipAE21"))].copy()
    ipc = ipc.reindex(ipc.rho_all19_tmolmean.abs().sort_values(ascending=False).index).head(16)
    rows = [["candidate", "ρ", "p", "WF ρ", "meaning"]] + [[Paragraph(E(r.candidate), st["cell"]), f2(r.rho_all19_tmolmean), f"{r.p_all19:.3f}", f2(r.rho_within_family), Paragraph(E(r.description)[:230], st["cell"])] for r in ipc.itertuples()]
    story += [P("Best ipSAE2 / Boltz-confidence candidates (n = 19)", "h2")]
    story += [mk_table(rows, [210, 34, 34, 60, FW - 338], colors, Table, TableStyle, st), PageBreak()]

    # ---------------------------------------------------------------- 3. statistics
    story += [P("3. How the statistics work", "h1"),
              P("<b>ρ (Spearman)</b>: rank correlation between a metric and log10 kcat/KM across proteins (range −1…+1). <b>p</b>: ordinary two-sided p-value, which is <i>too optimistic</i> here because thousands of metrics were tried. "
                "<b>q (BH)</b>: Benjamini–Hochberg false-discovery rate over all screened candidates. <b>Selection-aware p</b> (the one to trust): the activity labels are permuted 5,000–10,000 times, the best |ρ| in the same candidate set is recorded each time, "
                "and p is the fraction of permutations whose best |ρ| is at least the observed best — it answers 'how often does pure noise produce a candidate this good when I pick the best of this many?'. It is computed for the whole pool and separately for each feature set (tool)."),
              P("<b>Within-family ρ</b>: ranks are formed inside each scaffold family (HG3, 1A53, KABLE) and pooled, which removes between-scaffold effects (the families sit at different baseline values). <b>LOPO</b>: range of ρ when each protein is dropped in turn (robustness to single points). "
                "<b>AUROC</b>: how well the metric separates proteins above/below a kcat/KM threshold (0.5 = chance; for metrics with negative ρ read 1 − AUROC). <b>ρ old labels</b>: same correlation using the earlier screening-table activities (10 proteins). "
                "<b>Partial given n_mut</b>: ρ after removing the mutation count versus HG3 (checks that a metric is not just tracking how evolved a variant is). <b>Run 1 / run 2</b>: ρ from each of the two independent tmol relax runs; the reported ρ uses their mean."),
              P("<b>Test–retest (tmol relax noise)</b>: tmol relaxation is stochastic. The whole tmol protocol was run twice independently on the 19 proteins; the table shows how reproducible each quantity is (retest ρ and ICC between runs) "
                "and the correlation with activity in each run and in their mean. Quantities with high retest and the same sign in both runs are credible; the total-energy difference TSA − substrate fell from −0.82 to −0.36."), ]
    if retest is not None and len(retest):
        rt = retest.copy().head(22)
        hdr = ["quantity", "n", "retest ρ", "ICC", "sd across proteins", "ρ vs kcat run 1", "ρ vs kcat run 2", "ρ vs kcat mean"]
        rows = [hdr] + [[f"{r.quantity}", int(r.n), f"{r.retest_rho:.2f}", f"{r.ICC:.2f}", f"{r.sd_across_proteins:.2f}", f2(r.rho_kcat_run1), f2(r.rho_kcat_run2), f2(r.rho_kcat_mean)] for r in rt.itertuples()]
        story += [P("tmol test–retest (two independent relax runs)", "h2"), mk_table(rows, [190, 24, 50, 40, 90, 85, 85, 85], colors, Table, TableStyle, st), PageBreak()]
    else:
        story += [PageBreak()]

    # ---------------------------------------------------------------- 4. results: per feature set + protein table + figures
    story += [P("4. Results", "h1"), P("Selection-aware test by feature set (all 19 proteins; tmol = mean of two runs)", "h2")]
    story += [sel_table(N, colors, Table, TableStyle, st, FW), P("Same, HG3 lineage only (n = 12)", "h2"), sel_table(NH, colors, Table, TableStyle, st, FW), PageBreak()]
    cols = [("tmol swap|int.lk_ball_iso|tsa-prod", "lk_ball_iso TSA−prod"), ("tmol swap|int.fa_lk|tsa-sub", "fa_lk TSA−sub"), ("tmol swap|int.total|tsa-sub", "int.total TSA−sub"),
            ("pr.tsa_pk.int_sum", "per-res int. sum (TSA)"), ("bz.sub-prod.cat_to_pocket_pae", "Boltz cat→pocket PAE sub−prod"), ("ipSAE2/Boltz|boltz.ligand_iptm|sub-prod", "ligand ipTM sub−prod"),
            ("placer.apo.frac_pass_core", "PLACER apo pre-formed"), ("pl2.sub_holo.nac30", "PLACER sub NAC30")]
    src = lambda c: (M[c] if c in M.columns else (key.set_index("protein")[c] if c in key.columns else pd.Series(dtype=float)))
    kk = key.set_index("protein")
    hdr = ["protein", "family", "kcat/KM", "n_mut"] + [t for _, t in cols]
    rows = [hdr]
    for p_ in kk.index:
        rows.append([p_, kk.loc[p_, "family"], f"{kk.loc[p_, 'kcat_KM']:.3g}", "" if pd.isna(kk.loc[p_, "n_mut_vs_HG3"]) else int(kk.loc[p_, "n_mut_vs_HG3"])] +
                    [("" if (p_ not in src(c).index or pd.isna(src(c)[p_])) else f"{src(c)[p_]:.3g}") for c, _ in cols])
    story += [P("Dataset overview: the main metrics per protein (sorted by kcat/KM; full values in the CSV files)", "h2"), mk_table(rows, [58, 36, 52, 32] + [(FW - 178) / len(cols)] * len(cols), colors, Table, TableStyle, st), PageBreak()]
    figs = [("fig1_best_correlations.png", "Best correlations with activity (all 19 proteins)"), ("fig4_HG3_only.png", "HG3 lineage only (n = 12)"),
            ("fig2_preorganisation.png", "Preorganisation-type metrics"), ("fig3_signal_by_feature_set.png", "Signal per feature set versus the noise level")]
    for fn, ttl in figs:
        fp = Path(d["figs"]) / fn
        if fp.exists():
            iw, ih = ImageReader(str(fp)).getSize(); sc = min(FW / iw, 455 / ih)
            story += [P(ttl, "h2"), Image(str(fp), width=iw * sc, height=ih * sc), PageBreak()]

    # ---------------------------------------------------------------- 5. glossary
    story += [P("5. Quantity glossary", "h1"),
              P("Candidate ids have the form <b>tool|metric|arm</b> (e.g. <i>tmol swap|int.lk_ball_iso|tsa-prod</i>) or a dotted column name (e.g. <i>pr.tsa_pk.int_sum</i>). "
                "<b>Arm</b> is tsa / sub / prod (the quantity in that complex) or A-B (difference, A minus B). tmol/Rosetta energies are in Rosetta energy units (REU); lower = more favourable. "
                "Field values are in MV/cm with the sign convention that positive favours the developing charge separation. Tokens below are the building blocks of every column name.")]
    story += [P("Arms and tools", "h2")]
    rows = [["token", "meaning"]] + [[k, v] for k, v in TYPE.items()] + [[k, v] for k, v in TOOL.items()]
    story += [mk_table(rows, [150, FW - 150], colors, Table, TableStyle, st)]
    items = list(TOK.items()); half = (len(items) + 1) // 2
    rows = [["token", "meaning", "token", "meaning"]]
    for i in range(half):
        a_ = items[i]; b_ = items[i + half] if i + half < len(items) else ("", "")
        rows.append([Paragraph(f"<b>{E(a_[0])}</b>", st["cell"]), Paragraph(E(a_[1]), st["cell"]), Paragraph(f"<b>{E(b_[0])}</b>", st["cell"]), Paragraph(E(b_[1]), st["cell"])])
    story += [P("Tokens", "h2"), mk_table(rows, [100, FW / 2 - 100, 100, FW / 2 - 100], colors, Table, TableStyle, st, plain=True), PageBreak()]

    # ---------------------------------------------------------------- 6. big correlation table
    top = int(d["top"])
    story += [P(f"6. Correlation table: top {top} of {len(R)} candidates by |ρ| (all 19 proteins; full table in kemp_eval_correlations_all.csv)", "h1"),
              P("ρ = Spearman vs log10 kcat/KM (tmol = mean of two runs); run 1/run 2 = same with each relax run; WF = within-family ρ; LOPO = range when any one protein is dropped; AUROC at kcat/KM ≥ 10⁵; HG3 = ρ within the HG3 lineage only (n = 12). "
                "Shaded ρ cells: blue = positive, orange = negative; darker = stronger. Check run 1/run 2 agreement for tmol rows before believing them.", "small")]
    hdr = ["#", "candidate", "meaning", "set", "n", "ρ", "q", "run1", "run2", "WF", "LOPO", "AUROC", "HG3"]
    wid = [17, 128, 262, 42, 14, 28, 24, 28, 28, 28, 46, 32, 28]
    rows = [hdr]; shade = []
    T = R.head(top)
    for i, r in enumerate(T.itertuples(), 1):
        rows.append([i, Paragraph(E(r.candidate), st["tiny"]), Paragraph(E(r.description)[:300], st["tiny"]), Paragraph(E(r.set), st["tiny"]), "" if pd.isna(r.n) else int(r.n), f2(r.rho_all19_tmolmean),
                     "" if pd.isna(r.q_BH_all19) else f"{r.q_BH_all19:.2f}", f2(r.rho_run1), f2(r.rho_run2), f2(r.rho_within_family),
                     "" if pd.isna(r.rho_loo_min) else f"{r.rho_loo_min:+.2f}…{r.rho_loo_max:+.2f}", "" if pd.isna(r.auroc_100000) else f"{r.auroc_100000:.2f}", f2(r.rho_HG3_only)])
        v = r.rho_all19_tmolmean
        if pd.notna(v):
            a_ = min(abs(v), 1) * 0.65
            shade.append(("BACKGROUND", (5, i), (5, i), colors.Color(1 - a_ * 0.9, 1 - a_ * 0.45, 1.0) if v > 0 else colors.Color(1.0, 1 - a_ * 0.55, 1 - a_ * 0.9)))
    t = Table(rows, colWidths=wid, repeatRows=1)
    t.setStyle(TableStyle([("FONT", (0, 0), (-1, -1), "DV", 5.6), ("FONT", (0, 0), (-1, 0), "DVB", 6), ("BACKGROUND", (0, 0), (-1, 0), ACC), ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                           ("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 1), (-1, -1), 0.25, colors.HexColor("#e1e5ea")), ("TOPPADDING", (0, 0), (-1, -1), 1.2), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.2),
                           ("LEFTPADDING", (0, 0), (-1, -1), 2), ("RIGHTPADDING", (0, 0), (-1, -1), 2)] + shade))
    story += [t, PageBreak()]

    # ---------------------------------------------------------------- 7. interpretation / caveats
    story += [P("7. What the leading quantities mean (cautiously)", "h1"),
              P("<b>lk_ball_iso / fa_lk differences.</b> These are the solvation terms of the ligand–pocket pair energy. A more favourable TSA-minus-product (or TSA-minus-substrate) value means that, in the same pocket, the pocket accommodates the charge-delocalised, "
                "partially desolvated transition-state-like ligand better than the neutral product/substrate. That is the expected signature of a pocket tuned to the transition state, and it tracks directed-evolution progress in the HG3 lineage. "
                "It is a force-field proxy, with small absolute magnitudes (sd ≈ 0.6 REU across proteins for lk_ball_iso), which is why averaging two relax runs matters."),
              P("<b>Boltz PAE / ipTM differences.</b> The model is more self-consistent about the catalytic geometry in the substrate complex than in the product complex for the more active enzymes. This is a confidence statement, not energetics; "
                "it may reflect how closely the sequence resembles evolved Kemp-eliminase data Boltz has seen, and a single seed was used."),
              P("<b>PLACER / Rosetta / classic preorganisation.</b> PLACER re-docking is a floor effect (cannot tell ligands apart; prmsd ≈ 3 Å for all arms), its apo 'pre-formed' fraction is a weak positive trend (ρ ≈ 0.35, n.s.), and its HG3-derived windows are invalid for 1A53/KABLE. "
                "Rosetta dG_separated does not distinguish arms. Pocket motion between apo and holo is near zero for every protein."),
              P("8. Caveats", "h1"),
              ]
    cav = ["n = 19 proteins (12 in one lineage, 4 + 3 in two others) with activities from mixed literature sources (some back-calculated from figures); the PDB column of the manual dataset lists 29SB under HG3.17, probably misplaced; KABLE1.4's sequence contains the letter 'O' (→ UNK residue).",
           "The candidate pool is large (2,174). Selection-aware p-values account for this, but the choice to average the two tmol runs was made after seeing both, and the HG3-only analysis is a subset; treat the HG3 p-value as a lead until tested on new proteins.",
           "Two high-leverage points (1A53-2, kcat/KM 4; KABLE2.5, 3×10⁷) and scaffold clustering inflate some correlations; within-family and LOPO columns are given for this reason. Dropping both leverage points keeps the top tmol rows near |ρ| ≈ 0.7.",
           "tmol relax noise (several REU) is comparable to some effect sizes; only quantities that agree between the two runs should be trusted. PLACER (20 samples) and Rosetta are not replicated; Boltz seed variance is untested.",
           "Correlation is not mechanism: a metric can track evolutionary distance from HG3 rather than catalysis (n_mut vs kcat/KM ρ ≈ 0.55 in HG3); the partial-ρ column controls for this.",
           "Next steps that would settle it: (1) run the remaining ~9 proteins of the 28-protein set and, ideally, held-out designs with measured activity; (2) 2–3 Boltz seeds per protein; (3) pre-register a small panel (lk_ball_iso TSA−product, fa_lk TSA−substrate, interaction total TSA−substrate, Boltz PAE sub−prod) and test it once."]
    story += [Paragraph(c, st["bul"], bulletText="•") for c in cav]

    # ---------------------------------------------------------------- build
    def deco(canvas, doc):
        canvas.saveState(); canvas.setFont("DV", 7); canvas.setFillColor(MUTED)
        canvas.drawString(MG, 14, "Kemp-eliminase in-silico evaluation: metric/activity correlations"); canvas.drawRightString(W - MG, 14, f"page {doc.page}"); canvas.restoreState()
    doc = BaseDocTemplate(str(path), pagesize=landscape(A4), leftMargin=MG, rightMargin=MG, topMargin=MG, bottomMargin=MG, title="Kemp-eliminase in-silico evaluation", author="")
    doc.addPageTemplates([PageTemplate(id="p", frames=[Frame(MG, MG, FW, H - 2 * MG, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)], onPage=deco)])
    doc.build(story)


def mk_table(rows, widths, colors, Table, TableStyle, st, plain=False):
    from reportlab.platypus import Paragraph
    ACC = colors.HexColor("#0072B2")
    rows = [[(Paragraph(str(c).replace("&", "&amp;").replace("<", "&lt;"), st["cellb"]) if isinstance(c, str) else c) for c in rows[0]]] + rows[1:]
    t = Table(rows, colWidths=widths, repeatRows=1)
    t.setStyle(TableStyle([("FONT", (0, 0), (-1, -1), "DV", 6.6), ("FONT", (0, 0), (-1, 0), "DVB", 6.8), ("BACKGROUND", (0, 0), (-1, 0), ACC), ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                           ("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 1), (-1, -1), 0.25, colors.HexColor("#e1e5ea")), ("TOPPADDING", (0, 0), (-1, -1), 1.6), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6),
                           ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f6f8fa")])]))
    return t


def sel_table(N, colors, Table, TableStyle, st, FW):
    from xml.sax.saxutils import escape
    from reportlab.platypus import Paragraph
    df = N.copy()
    rows = [["feature set", "candidates", "best |ρ|", "selection-aware p", "noise 95th pct", "best candidate"]]
    for r in df.itertuples():
        rows.append([r.group, int(r.n_candidates), f"{r.best_abs_rho:.2f}", f"{r.p_selection_aware:.3f}", f"{r.null95:.2f}", Paragraph(escape(str(r.best_candidate)), st["cell"])])
    t = mk_table(rows, [100, 60, 55, 90, 80, FW - 385], colors, Table, TableStyle, st)
    for i, r in enumerate(df.itertuples(), 1):
        if r.p_selection_aware < 0.05:
            t.setStyle(TableStyle([("BACKGROUND", (3, i), (3, i), colors.HexColor("#f8d7cf"))]))
    return t
