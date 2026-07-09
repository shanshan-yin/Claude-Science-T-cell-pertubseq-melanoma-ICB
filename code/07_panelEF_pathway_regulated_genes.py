"""
Panels E/F: pathway-regulated signature genes

Extracted reproduction code (Claude Science lineage).
Source artifact version: 7b6241b9-c074-41c3-a92d-36f0e5482ae4
Conda environment: scrnaseq-perturbseq
"""

import os, tempfile
os.environ["NUMBA_CACHE_DIR"] = tempfile.mkdtemp()

import numpy as np
import pandas as pd
import anndata as ad
import io
import urllib.request
import h5py
import matplotlib as mpl
from matplotlib import rcParams

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


# Load atlas
atlas = ad.read_h5ad("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/929e269b-6f1b-4fe9-9481-1d3dffb36182/vb2a694f5_GWCD4i_zscore_subset.h5ad")
if "gene_name" in atlas.var.columns:
    atlas.var_names = atlas.var["gene_name"].astype(str).values

# Load signature
sig = pd.read_csv("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/ad421b73-9ff3-4728-8bae-bc41272463c0/vf61e5db3_CD4_RvsNR_signature.csv")
if "gene" in sig.columns:
    sig = sig.set_index("gene")

np.random.seed(214)

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
    assert (X.index == yv.index).all()
    return mask_cis(X), yv

from pert2state_model.Perturb2StateModel import Perturb2StateModel

conditions = ["Rest", "Stim8hr", "Stim48hr"]
models_only = {}
for cond in conditions:
    X, y = build_Xy(cond, "z_discovery_tumor")
    m_fit = Perturb2StateModel(pca_transform=True, n_pcs=60, n_splits=5, n_repeats=1, positive=False)
    m_fit.fit(X, y, model_id=f"tumor_{cond}")
    models_only[("tumor", cond)] = m_fit

mt = models_only[("tumor", "Rest")]
coefs = mt.get_coefs().copy()

import scipy.cluster.hierarchy as sch

Xr, yr = build_Xy("Rest", "z_discovery_tumor")

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

from scipy.stats import mannwhitneyu

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

# Panel D top regulators
top_pos = cp2["coef_mean"].nlargest(15).index.tolist()
top_neg = cp2["coef_mean"].nsmallest(15).index.tolist()
panelD_regs = list(dict.fromkeys(top_pos + top_neg + saga_pl + ifn_pl))

# Remote atlas streaming to fetch adj_p
url = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad"

class HTTPRangeFile(io.RawIOBase):
    def __init__(self, url):
        self.url = url
        self.pos = 0
        self.size = int(urllib.request.urlopen(urllib.request.Request(url, method="HEAD")).headers["Content-Length"])
    def readable(self): return True
    def seekable(self): return True
    def seek(self, o, w=0):
        self.pos = {0: o, 1: self.pos + o, 2: self.size + o}[w]
        return self.pos
    def tell(self): return self.pos
    def readinto(self, b):
        n = len(b)
        if self.pos >= self.size: return 0
        end = min(self.pos + n - 1, self.size - 1)
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        for a in range(4):
            try:
                d = urllib.request.urlopen(req).read()
                break
            except Exception:
                if a == 3: raise
        b[:len(d)] = d
        self.pos += len(d)
        return len(d)

def decode_cat(grp):
    if isinstance(grp, h5py.Group) and "categories" in grp:
        return grp["categories"][:].astype(str)[grp["codes"][:]]
    return grp[:].astype(str)

f = h5py.File(HTTPRangeFile(url), "r")
tc = decode_cat(f["obs"]["target_contrast_gene_name"])
cc = decode_cat(f["obs"]["culture_condition"])
var_names = decode_cat(f["var"]["gene_name"]) if "gene_name" in f["var"] else decode_cat(f["var"]["_index"])
adjp = f["layers"]["adj_p_value"]
rows = {}
for reg in panelD_regs:
    idx = np.where((tc == reg) & (cc == "Rest"))[0]
    if len(idx):
        rows[reg] = adjp[int(idx[0]), :]
f.close()
adjp_all = pd.DataFrame(rows, index=var_names).T  # regs x all genes

n_reg_by = (adjp_all < 0.1).sum(axis=0)
res_regs = [r for r in top_pos + saga_pl if r in adjp_all.index]
resp_regs = [r for r in top_neg + ifn_pl if r in adjp_all.index]
n_by_res = (adjp_all.loc[res_regs] < 0.1).sum(axis=0)
n_by_resp = (adjp_all.loc[resp_regs] < 0.1).sum(axis=0)

vv = pd.DataFrame({
    "log2FoldChange": sig["log2FC_tumor"],
    "padj": sig["padj_tumor"],
}).dropna()
vv["neglog10padj"] = -np.log10(vv["padj"].clip(lower=1e-300))
vv["n_reg"] = vv.index.map(n_reg_by.to_dict()).fillna(0)
vv["n_res_reg"] = vv.index.map(n_by_res.to_dict()).fillna(0)
vv["n_resp_reg"] = vv.index.map(n_by_resp.to_dict()).fillna(0)

reg_sig = vv[(vv["n_reg"] >= 1) & (vv["padj"] < 0.05) & (vv["log2FoldChange"].abs() > 0.5)]

# SAGA-specific panel
saga_regs_present = [r for r in saga_pl if r in adjp_all.index]
n_by_saga = (adjp_all.loc[saga_regs_present] < 0.1).sum(axis=0)
vv["n_saga"] = vv.index.map(n_by_saga.to_dict()).fillna(0)
saga_reg_sig = vv[(vv["n_saga"] >= 1) & (vv["padj"] < 0.05) & (vv["log2FoldChange"].abs() > 0.5)]

# Fetch KD z-scores for SAGA-regulated genes
rest = atlas[atlas.obs["culture_condition"] == "Rest"].copy()
reg_arr = rest.obs["target_contrast_gene_name"].astype(str).values
X_rest = rest.X
X_rest = np.asarray(X_rest.todense()) if hasattr(X_rest, "todense") else np.asarray(X_rest)
gidx = {g: i for i, g in enumerate(rest.var_names)}

saga_regs = ['USP22', 'ATXN7L3', 'SGF29', 'TADA1', 'TADA2B', 'TADA3', 'TAF5L', 'TAF6L', 'SUPT7L', 'SUPT20H']

def per_gene_kd(regs, genes):
    rows_mask = np.isin(reg_arr, regs)
    return {g: float(np.nanmean(X_rest[rows_mask, gidx[g]])) for g in genes if g in gidx}

saga_df = saga_reg_sig[["log2FoldChange", "padj", "n_saga"]].copy()
saga_df.index.name = "gene"
saga_df = saga_df.reset_index()
saga_df["SAGA_KDz"] = saga_df["gene"].map(per_gene_kd(saga_regs, saga_df["gene"].tolist()))

saga_df.to_csv("panelE_regulated_genes.csv", index=False)
print("saved panelE_regulated_genes.csv:", saga_df.shape)