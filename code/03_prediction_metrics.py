"""
Held-out prediction accuracy per condition

Extracted reproduction code (Claude Science lineage).
Source artifact version: 129199bf-e50d-44e3-b63b-2f7eab431380
Conda environment: scrnaseq-perturbseq
"""

import os, tempfile
os.environ["NUMBA_CACHE_DIR"] = tempfile.mkdtemp()

import anndata, numpy as np, pandas as pd
from pert2state_model.Perturb2StateModel import Perturb2StateModel

np.random.seed(214)

atlas = anndata.read_h5ad("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/929e269b-6f1b-4fe9-9481-1d3dffb36182/va6baf95d_GWCD4i_zscore_subset.h5ad")
atlas.var_names = atlas.var["gene_name"].astype(str).values
sig = pd.read_csv("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/ad421b73-9ff3-4728-8bae-bc41272463c0/vf61e5db3_CD4_RvsNR_signature.csv").set_index("gene")

def mask_cis(Xp):
    Xm = Xp.copy()
    for g in [i for i in Xm.index if i in Xm.columns]:
        Xm.loc[g, g] = np.nan
    return Xm.fillna(Xm.mean())

def build_Xy(cond, sig_col):
    m = (atlas.obs["culture_condition"] == cond).values
    obs = atlas.obs[m]; Xall = atlas.X[m].astype(np.float64)
    tgt = obs["target_contrast_gene_name"].astype(str).values
    keep = ~pd.Series(tgt).duplicated().values
    tgt = tgt[keep]; Xall = Xall[keep]
    y = sig[sig_col].dropna()
    common = np.intersect1d(y.index.values, atlas.var_names.values)
    gi = pd.Index(atlas.var_names).get_indexer(common)
    X = pd.DataFrame(Xall[:, gi].T, index=common, columns=tgt).astype(np.float64)
    yv = y.loc[common].astype(np.float64); yv = (yv - yv.mean()) / yv.std()
    return mask_cis(X), yv

conditions = ["Rest", "Stim8hr", "Stim48hr"]
sig_cols = {"tumor": "z_discovery_tumor", "blood": "z_blood_bladder"}
evals = []

for sname, scol in sig_cols.items():
    for cond in conditions:
        X, y = build_Xy(cond, scol)
        m = Perturb2StateModel(pca_transform=True, n_pcs=60, n_splits=5, n_repeats=1, positive=False)
        m.fit(X, y, model_id=f"{sname}_{cond}")
        se = m.summarize_eval().iloc[0]
        evals.append({
            "signature": sname, "condition": cond,
            "n_genes": X.shape[0], "n_perturb": X.shape[1],
            "test_pearson": se["test_pearson"], "test_pearson_se": se["test_pearson_se"],
            "test_spearman": se["test_spearman"], "test_r2": se["test_r2"],
            "train_pearson": se["train_pearson"]
        })

ev_df = pd.DataFrame(evals)
os.makedirs("data/models", exist_ok=True)
ev_df.to_csv("prediction_metrics.csv", index=False)
print(ev_df.to_string(index=False))