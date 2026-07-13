"""
Assemble the 8-panel Figure 5 composite (A-H).

Tiles the eight committed panel PNGs in results/figures/ into the final
figure5_composite.png. The individual panels are produced by scripts 04-07
(and the schematic / ranked-regulator panels supplied alongside them).

Conda environment: scrnaseq-perturbseq
"""

import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import sys, glob, matplotlib.font_manager as fm
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
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.size": base,
        "axes.labelsize": base, "axes.titlesize": base,
        "legend.fontsize": secondary, "xtick.labelsize": tick, "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid), "legend.frameon": False,
        "figure.dpi": 200, "savefig.dpi": 300, "savefig.bbox": "tight",
        "axes.titleweight": "normal", "axes.titlelocation": "left",
        "axes.labelweight": "normal", "lines.linewidth": 1.2, "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    if font:
        mpl.rcParams["font.sans-serif"] = [font, "DejaVu Sans"]


FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "figures")

# Panel letter -> committed PNG filename
panels = {
    "A": "panelA_schematic.png",         # analysis schematic
    "B": "panelB_volcano.png",           # R-vs-NR signature volcano
    "C": "panelC_prediction.png",        # held-out prediction accuracy
    "D": "panelD_ranked_regulators.png", # ranked regulators by pathway
    "E": "panelD_regulator_effects.png", # predicted regulator effects
    "F": "panelE_saga_volcano.png",      # SAGA-pathway targets
    "G": "panelF_ifn_volcano.png",       # interferon-pathway targets
    "H": "panelG_drug_targets.png",      # nominated drug targets
}

apply_figure_style()
fig = plt.figure(figsize=(16, 30))
gs = fig.add_gridspec(5, 2, height_ratios=[1.0, 1.0, 0.78, 1.0, 0.85], hspace=0.12, wspace=0.05)
layout = [("A", (0, 0)), ("B", (0, 1)), ("C", (1, 0)), ("D", (1, 1)),
          ("E", (2, slice(None))), ("F", (3, 0)), ("G", (3, 1)), ("H", (4, slice(None)))]
for lab, pos in layout:
    ax = fig.add_subplot(gs[pos])
    ax.imshow(mpimg.imread(os.path.join(FIG_DIR, panels[lab])))
    ax.axis("off")
    ax.text(-0.02, 1.02, lab, transform=ax.transAxes, fontsize=20, fontweight="bold", va="bottom", ha="left")
fig.suptitle("Figure 5.  Perturb-seq regulators that reproduce the CD4\u207a T-cell ICB-resistance signature",
             fontsize=22, fontweight="bold", y=0.905, x=0.5)
fig.savefig(os.path.join(FIG_DIR, "figure5_composite.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("recomposed 8-panel figure")
