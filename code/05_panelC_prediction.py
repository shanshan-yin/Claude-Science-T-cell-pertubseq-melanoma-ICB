"""
Panel C: prediction accuracy bars

Extracted reproduction code (Claude Science lineage).
Source artifact version: 0466839b-cd8d-4f73-beb1-4fbfecb85581
Conda environment: scrnaseq-perturbseq
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl

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


c = pd.read_csv("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/e5acbf1a-cde2-4af4-906b-b7d74905f034/v3e3b4a50_panelC_accuracy_tumor.csv")

apply_figure_style()
COL = "#c51b8a"
UB = 0.285
conds = ["Rest", "Stim 8 hr", "Stim 48 hr"]
r_over = c["r_overall"].values
sem = c["sem_overall"].values
fig, ax = plt.subplots(figsize=(6.4, 6.0))
x = np.arange(len(conds))
ax.bar(x, r_over, yerr=sem, width=0.6, color=COL, edgecolor="none", error_kw=dict(lw=1.1, capsize=3))
ax.axhline(UB, ls="--", lw=0.9, color="grey")
ax.text(len(conds)-0.5, UB+0.012, "within-tumor replication ceiling (r=0.29)", fontsize=6.4, color="grey", ha="right", va="bottom")
for xi, (v, e) in enumerate(zip(r_over, sem)):
    ax.text(xi, v+e+0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=7.5, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(conds)
ax.set_ylim(0, 0.7)
ax.set_ylabel("held-out prediction accuracy\n(Pearson r)")
ax.set_title("Perturb-seq predicts the melanoma-tumor CD4 R-vs-NR signature\n(overall signature, all 9,447 genes)", fontsize=8.4, loc="left")
fig.savefig("panelC_prediction.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("Panel C single-chart:", [f"{v:.3f}" for v in r_over])