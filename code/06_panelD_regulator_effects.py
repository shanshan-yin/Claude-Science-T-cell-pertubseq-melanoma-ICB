"""
Panel D: predicted regulator effects

Extracted reproduction code (Claude Science lineage).
Source artifact version: 776eddcb-2073-407e-bd11-280549cf9341
Conda environment: scrnaseq-perturbseq
"""

import os
import tempfile
os.environ["NUMBA_CACHE_DIR"] = tempfile.mkdtemp()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import anndata
from adjustText import adjust_text
from scipy.stats import mannwhitneyu
import scipy.cluster.hierarchy as sch

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


saga_pl, saga_p = enr(saga_module, "SAGA/chromatin DUB module")
ifn_pl, ifn_p = enr(ifn_module, "Interferon module")

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
print("saved Panel D with corrected title (overall R-vs-NR signature)")