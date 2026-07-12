"""
Panel B: R-vs-NR volcano

Extracted reproduction code (Claude Science lineage).
Source artifact version: b52c5cef-4149-4e97-b2f8-e087beda095c
Conda environment: scrnaseq-perturbseq
"""

import os, tempfile
os.environ["NUMBA_CACHE_DIR"] = tempfile.mkdtemp()

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

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


sig = pd.read_csv("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/ad421b73-9ff3-4728-8bae-bc41272463c0/vf61e5db3_CD4_RvsNR_signature.csv")

apply_figure_style()
s = sig.copy()
s["y"] = -np.log10(s["padj_tumor"].clip(lower=1e-300))
lfc = s["log2FC_tumor"]
sig_up_nr = (s["padj_tumor"] < 0.05) & (lfc > 0.5)
sig_up_r = (s["padj_tumor"] < 0.05) & (lfc < -0.5)
fig, ax = plt.subplots(figsize=(7.6, 6.0))
cut = -np.log10(0.05)
ax.scatter(lfc[~(sig_up_nr | sig_up_r)], s["y"][~(sig_up_nr | sig_up_r)], s=5, c="#dcdcdc", alpha=0.5, rasterized=True, linewidths=0)
ax.scatter(lfc[sig_up_nr], s["y"][sig_up_nr], s=10, c="#d94801", alpha=0.75, linewidths=0, label=f"up in non-responder (n={int(sig_up_nr.sum())})")
ax.scatter(lfc[sig_up_r], s["y"][sig_up_r], s=10, c="#08519c", alpha=0.75, linewidths=0, label=f"up in responder (n={int(sig_up_r.sum())})")
ax.axhline(cut, ls="--", lw=0.7, color="grey")
ax.axvline(0, ls="-", lw=0.5, color="grey")
ax.set_ylim(top=s["y"].max() * 1.05)
ax.set_xlabel("log₂FC  (Non-responder / Responder)")
ax.set_ylabel("−log₁₀ adjusted p")
ax.set_title("CD4⁺ T-cell R-vs-NR signature (melanoma tumor)", fontsize=8.6, loc="left")
ax.legend(frameon=False, fontsize=6.8, loc="upper left")
fig.savefig("panelB_volcano.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"clean Panel B: up-NR={int(sig_up_nr.sum())}, up-R={int(sig_up_r.sum())}, no gene/pathway labels")