"""
Assemble the 7-panel Figure 5 composite

Extracted reproduction code (Claude Science lineage).
Source artifact version: 9b103f16-d49c-46e7-ad74-d9f286b67f5f
Conda environment: scrnaseq-perturbseq
"""

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
import os
import sys
import glob


META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


# Panel D: rebuild from lineage code (swap "module" -> "pathway")
import numpy as np
import pandas as pd
import anndata
from adjustText import adjust_text
from scipy.stats import mannwhitneyu
import scipy.cluster.hierarchy as sch
import os, tempfile
os.environ["NUMBA_CACHE_DIR"] = tempfile.mkdtemp()

np.random.seed(214)

atlas = anndata.read_h5ad("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/929e269b-6f1b-4fe9-9481-1d3dffb36182/vd08bb0dd_GWCD4i_zscore_subset.h5ad")
atlas.var_names = atlas.var["gene_name"].astype(str).values
sig = pd.read_csv("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/ad421b73-9ff3-4728-8bae-bc41272463c0/vf61e5db3_CD4_RvsNR_signature.csv").set_index("gene")

from pert2state_model.Perturb2StateModel import Perturb2StateModel


def mask_cis(Xp):
    Xm = Xp.copy()
    for g in [i for i in Xm.index if i in Xm.columns]:
        Xm.loc[g, g] = np.nan
    return Xm.fillna(Xm.mean())


def build_Xy(cond, sig_col):
    m = (atlas.obs["culture_condition"] == cond).values
    obs = atlas.obs[m]
    Xall = atlas.X[m].astype(np.float64)
    tgt = obs["target_contrast_gene_name"].astype(str).values
    keep = ~pd.Series(tgt).duplicated().values
    tgt = tgt[keep]
    Xall = Xall[keep]
    y = sig[sig_col].dropna()
    common = np.intersect1d(y.index.values, atlas.var_names.values)
    gi = pd.Index(atlas.var_names).get_indexer(common)
    X = pd.DataFrame(Xall[:, gi].T, index=common, columns=tgt).astype(np.float64)
    yv = y.loc[common].astype(np.float64)
    yv = (yv - yv.mean()) / yv.std()
    return mask_cis(X), yv


Xr, yr = build_Xy("Rest", "z_discovery_tumor")
mt = Perturb2StateModel(pca_transform=True, n_pcs=60, n_splits=5, n_repeats=1, positive=False)
mt.fit(Xr, yr, model_id="tumor_Rest")

coefs = mt.get_coefs().copy()

sig_gs = [g for g in sig.index[sig["padj_tumor"] < 0.01] if g in Xr.index]
df2cluster = Xr.T[sig_gs]
link = sch.linkage(df2cluster.values, method="ward")
ordered_regs = df2cluster.index[sch.leaves_list(link)].tolist()
cp = coefs.loc[ordered_regs].copy()
cp["rank"] = range(len(cp))
cp = cp.dropna(subset=["coef_mean"])

cp2 = cp.copy()
cp2.index = [i.split("_")[0] for i in cp2.index]

saga_module = ["USP22", "ATXN7L3", "SGF29", "ENY2", "TADA1", "TADA2B", "TADA3", "TAF5L", "TAF6L", "SUPT7L", "SUPT20H", "KAT2A", "KAT2B", "CCDC101"]
ifn_module = ["IRF1", "IRF9", "IFNAR1", "IFNAR2", "STAT1", "STAT2", "JAK1", "JAK2", "TYK2", "IRF2", "STAT3"]


def enr(module, name):
    present = [g for g in module if g in cp2.index]
    abs_mod = np.abs(cp2.loc[present, "coef_mean"])
    abs_all = np.abs(cp2["coef_mean"])
    if len(present) >= 3:
        _, p = mannwhitneyu(abs_mod, abs_all, alternative="greater")
    else:
        p = np.nan
    return present, p


saga_pl, saga_p = enr(saga_module, "SAGA/chromatin DUB pathway")
ifn_pl, ifn_p = enr(ifn_module, "Interferon pathway")

apply_figure_style()
top20 = pd.concat([cp2["coef_mean"].nsmallest(10), cp2["coef_mean"].nlargest(10)])
fig, ax = plt.subplots(figsize=(11, 5.2))
ax.errorbar(cp2["rank"], cp2["coef_mean"], yerr=cp2["coef_sem"], fmt=".", color="grey", markersize=3, alpha=0.4, elinewidth=1, capsize=0)
ax.errorbar(cp2.loc[top20.index, "rank"], cp2.loc[top20.index, "coef_mean"], yerr=cp2.loc[top20.index, "coef_sem"], fmt=".", color="#404040", markersize=7, elinewidth=2.2, capsize=0, alpha=0.95)
for mod, color, name, p in [(saga_pl, "#c51b8a", "SAGA/chromatin DUB pathway", saga_p), (ifn_pl, "#2c7fb8", "interferon pathway", ifn_p)]:
    ax.errorbar(cp2.loc[mod, "rank"], cp2.loc[mod, "coef_mean"], yerr=cp2.loc[mod, "coef_sem"], fmt=".", color=color, markersize=8, elinewidth=2.5, capsize=0, label=f"{name} (Mann-Whitney p = {p:.1g})", zorder=5)
ax.axhline(0, color="black", ls=":", alpha=1, lw=0.8)
ann = list(dict.fromkeys(list(top20.index) + saga_pl + ifn_pl))
texts = []
for g in ann:
    col = "#c51b8a" if g in saga_pl else ("#2c7fb8" if g in ifn_pl else "#404040")
    texts.append(ax.text(cp2.loc[g, "rank"], cp2.loc[g, "coef_mean"], g, fontsize=6.5, fontweight="bold", color=col, fontstyle="italic"))
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="->", color="black", alpha=0.5, lw=0.5))
ax.set_xlabel("Regulators (clustered by effect on R-vs-NR signature genes)")
ax.set_ylabel("Predicted regulator\neffect  ($w_r$)")
ax.set_title("Predicted regulator effects on the overall CD4 R-vs-NR signature (tumor model, Rest)", fontsize=9, loc="left")
ax.legend(frameon=False, fontsize=7, loc="upper left")
fig.savefig("panelD_regulator_effects.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# Panel G: rebuild from lineage code (swap "module" -> "pathway")
import json
from matplotlib.lines import Line2D

tg = json.load(open("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/16716cfb-d54e-4644-9a01-cba9ca73ca22/v7449d16d_chembl_targets.json"))
cands = {c["gene"]: c for c in json.load(open("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/b3eaf633-3a7d-42a9-9298-010a7b672a36/veeb6ed5c_candidates.json"))}
dd = json.load(open("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/c709b381-6a17-4e4a-a541-98d2a6ed8140/v048919e9_drug_details.json"))
ifn = json.load(open("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/86849dd1-e599-42b0-aacc-a569710831dd/v50898db3_ifnar_drugs.json"))


def clean_name(n):
    if not n:
        return None
    n = n.lower()
    for salt in [" phosphate", " hydrochloride", " citrate", " maleate", " hemihydrate", " malate", " sulfate"]:
        n = n.replace(salt, "")
    return n.strip()


def approved_leads(gene, n_show=3):
    ds = dd.get(gene, {}).get("drugs", [])
    names = []
    for d in ds:
        if d.get("max_phase") == 4 and d.get("name"):
            cn = clean_name(d["name"])
            if cn and cn not in names:
                names.append(cn)
    return names[:n_show], len([d for d in ds if d.get("max_phase") == 4])


jak2_leads, _ = approved_leads("JAK2", 3)
tyk2_leads, _ = approved_leads("TYK2", 3)
lead_map = {"JAK2": ", ".join(jak2_leads), "TYK2": ", ".join(tyk2_leads), "IFNAR1": "anifrolumab"}

rows = []
for g, c in cands.items():
    n_appr = len([d for d in dd.get(g, {}).get("drugs", []) if d.get("max_phase") == 4])
    tid = tg.get(g, {}).get("target_chembl_id")
    ttype = tg.get(g, {}).get("target_type")
    if g == "IFNAR1":
        n_appr = len([d for d in ifn.get("IFNAR1", []) if d.get("max_phase") == 4 and d.get("action") == "ANTAGONIST"])
    tier = ("approved drug" if n_appr > 0 else ("target, no directional drug" if tid else "no ChEMBL target"))
    rows.append({"gene": g, "wr": c["wr"], "rank": c["rank"], "ifn_module": c["ifn_module"],
                 "target_chembl_id": tid, "target_type": ttype, "n_approved": n_appr,
                 "lead_drugs": lead_map.get(g, ""), "tier": tier})
dt = pd.DataFrame(rows).sort_values("wr")

apply_figure_style()
IFN_COLOR = "#2c7fb8"
jak2_show = jak2_leads[:4]
d = dt[dt["wr"] < 0].sort_values("wr").head(20).iloc[::-1].reset_index(drop=True)
tier_color = {"approved drug": "#238b45", "target, no directional drug": "#feb24c", "no ChEMBL target": "#d9d9d9"}
fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 6.4), gridspec_kw={"width_ratios": [1.5, 1.0]})
y = np.arange(len(d))
for i, r in d.iterrows():
    axL.plot([0, abs(r["wr"])], [i, i], color="#bdbdbd", lw=1, zorder=1)
    axL.scatter(abs(r["wr"]), i, s=95, color=tier_color[r["tier"]], edgecolor=(IFN_COLOR if r["ifn_module"] else "white"),
                linewidth=(1.8 if r["ifn_module"] else 0.5), zorder=3)
axL.set_yticks(y)
axL.set_yticklabels([f'{g}' + (" ✦" if m else "") for g, m in zip(d["gene"], d["ifn_module"])], fontsize=7.5)
for tl, m in zip(axL.get_yticklabels(), d["ifn_module"]):
    if m:
        tl.set_color(IFN_COLOR)
        tl.set_fontweight("bold")
axL.set_xlabel("|predicted regulator effect $w_r$|  (inhibit → responder-like)")
axL.set_title("Model-nominated targets: inhibit to sensitize non-responders", fontsize=8.6, loc="left")
axL.legend(handles=[Line2D([], [], marker='o', ls='', mfc=tier_color["approved drug"], mec='none', ms=8, label="approved inhibitor exists"),
                    Line2D([], [], marker='o', ls='', mfc=tier_color["target, no directional drug"], mec='none', ms=8, label="ChEMBL target, no directional drug"),
                    Line2D([], [], marker='o', ls='', mfc=tier_color["no ChEMBL target"], mec='none', ms=8, label="no ChEMBL target"),
                    Line2D([], [], marker='o', ls='', mfc='white', mec=IFN_COLOR, mew=1.8, ms=8, label="interferon pathway (✦)")],
           frameon=False, fontsize=6.3, loc="lower right")
axL.margins(y=0.02)
axR.axis("off")
axR.set_title("Actionable interferon-axis drugs", fontsize=8.6, loc="left")
cards = [
    ("IFNAR1", "anti-IFNAR1 receptor antibody", "anifrolumab (approved, SLE)", "ANTAGONIST — blocks type-I IFN receptor", 1),
    ("TYK2", "JAK-family kinase (type-I IFN)", ", ".join(tyk2_leads) + "\n(+ brepocitinib Ph3)", "JAK-family INHIBITOR — 4 approved", 4),
    ("JAK2", "JAK-family kinase (IFN-γ/ISG)", ", ".join(jak2_show[:2]) + ",\n" + ", ".join(jak2_show[2:4]), "INHIBITOR — 12 approved", 12),
]
yc = 0.92
for gene, tclass, drugs, moa, n in cards:
    axR.add_patch(plt.Rectangle((0.02, yc - 0.24), 0.96, 0.24, transform=axR.transAxes, fc="#eef6fb", ec=IFN_COLOR, lw=1.2, zorder=1))
    axR.text(0.06, yc - 0.03, gene, transform=axR.transAxes, fontsize=11, fontweight="bold", color=IFN_COLOR, va="top")
    axR.text(0.06, yc - 0.085, tclass, transform=axR.transAxes, fontsize=6.6, fontstyle="italic", color="#333", va="top")
    axR.text(0.06, yc - 0.125, moa, transform=axR.transAxes, fontsize=6.4, color="#555", va="top")
    axR.text(0.06, yc - 0.165, drugs, transform=axR.transAxes, fontsize=6.8, color="#111", va="top")
    axR.text(0.94, yc - 0.03, f"{n}✓", transform=axR.transAxes, fontsize=9, fontweight="bold", color="#238b45", va="top", ha="right")
    yc -= 0.30
axR.text(0.02, 0.02, "✓ = approved inhibitors on this ChEMBL target (mechanism-of-action; drug names as retrieved)",
         transform=axR.transAxes, fontsize=5.6, color="#777", style="italic")
fig.savefig("panelG_drug_targets.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# Compose the final figure
apply_figure_style()
fig = plt.figure(figsize=(16, 26))
gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.0, 1.0, 0.92], hspace=0.11, wspace=0.05)


def add(pos, img, lab):
    ax = fig.add_subplot(pos)
    ax.imshow(mpimg.imread(img))
    ax.axis("off")
    ax.text(-0.02, 1.02, lab, transform=ax.transAxes, fontsize=17, fontweight="bold", va="bottom", ha="left")


add(gs[0, 0], "/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/ed5d7c05-1bad-452c-ad45-026ce6c5ddd0/v9ac8acef_panelA_schematic.png", "A")
add(gs[0, 1], "/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/2deb0527-5e89-4862-b186-fd7f0991ebdd/vb52c5cef_panelB_volcano.png", "B")
add(gs[1, 0], "/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/5998d4ab-5c9e-4c47-aab6-858cacf750e3/v0466839b_panelC_prediction.png", "C")
add(gs[1, 1], "panelD_regulator_effects.png", "D")
add(gs[2, 0], "/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/fffafb2f-b375-45ed-988f-3392daa1b78d/v94565d35_panelE_regulators.png", "E")
add(gs[2, 1], "/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/4bd403f6-1805-42df-9c98-38dad1956d9e/v465a64f7_panelF_ifn_volcano.png", "F")
add(gs[3, :], "panelG_drug_targets.png", "G")
fig.suptitle("Figure 5.  Perturb-seq regulators that reproduce the CD4⁺ T-cell ICB-resistance signature", fontsize=15, y=0.902, x=0.5)
fig.savefig("figure5_composite.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("recomposed")