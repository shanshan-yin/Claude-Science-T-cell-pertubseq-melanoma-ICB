"""
Build CD4 R-vs-NR signature (per-cohort pseudobulk DESeq2)

Extracted reproduction code (Claude Science lineage).
Source artifact version: f61e5db3-a9cf-450b-be08-f83138ef1198
Conda environment: scrnaseq-perturbseq
"""

import anndata
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
import os

os.makedirs("data/signature", exist_ok=True)

def pseudobulk(ad, min_cells=10):
    samples = ad.obs["sample_id"].value_counts()
    keep = samples[samples >= min_cells].index.tolist()
    a = ad[ad.obs["sample_id"].isin(keep)]
    rows = {}; meta = []
    for sid in keep:
        m = a.obs["sample_id"].values == sid
        rows[sid] = np.asarray(a.X[m].sum(0)).ravel()
        o = a.obs[m].iloc[0]
        meta.append({"sample_id": sid, "response": str(o["response"]), "patient": str(o["patient"]),
                     "therapy": str(o.get("therapy", "NA")), "n_cells": int(m.sum())})
    pb = pd.DataFrame(rows, index=a.var_names).T
    cd = pd.DataFrame(meta).set_index("sample_id").loc[pb.index]
    return pb, cd

def run_deseq(pb, cd, cohort):
    counts = pb.round().astype(int)
    keep = (counts > 0).mean(0) >= 0.25
    counts = counts.loc[:, keep]
    cd = cd.copy()
    cd["condition"] = cd["response"].map({"Responder": "R", "Non-responder": "NR"})
    inf = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(counts=counts, metadata=cd, design="~condition", inference=inf, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "NR", "R"], inference=inf, quiet=True)
    st.summary()
    res = st.results_df.copy()
    res["zscore"] = res["log2FoldChange"] / res["lfcSE"]
    res["cohort"] = cohort
    return res

sf = anndata.read_h5ad("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/35bb6cd1-0570-4042-815c-7d77b3c443ba/vdde0852e_SadeFeldman_CD4.h5ad")
bl = anndata.read_h5ad("/Users/syin/.claude-science/orgs/c8e50097-6b7c-40ca-b8eb-60f27d93ee37/artifacts/proj_50aa59908860/72e5e4af-1525-4869-8a9a-af0a0a241a85/va64923e7_Bladder_CD4.h5ad")

pb_sf, cd_sf = pseudobulk(sf)
pb_bl, cd_bl = pseudobulk(bl)

res_sf = run_deseq(pb_sf, cd_sf, "SadeFeldman_melanoma")
res_bl = run_deseq(pb_bl, cd_bl, "Bladder_atezolizumab")

# Within-tumor split for replication upper bound
rng = np.random.RandomState(214)
cd = cd_sf.copy()
pats = cd.groupby("response")["patient"].unique()
half1 = set(); half2 = set()
for resp, pl in pats.items():
    pl = list(pl); rng.shuffle(pl)
    half1 |= set(pl[::2]); half2 |= set(pl[1::2])
m1 = cd["patient"].isin(half1); m2 = cd["patient"].isin(half2)

def deseq_z(counts, cd):
    counts = counts.round().astype(int)
    keep = (counts > 0).mean(0) >= 0.25; counts = counts.loc[:, keep]
    cd = cd.copy(); cd["condition"] = cd["response"].map({"Responder": "R", "Non-responder": "NR"})
    inf = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(counts=counts, metadata=cd, design="~condition", inference=inf, quiet=True); dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "NR", "R"], inference=inf, quiet=True); st.summary()
    return (st.results_df["log2FoldChange"] / st.results_df["lfcSE"]).rename("z")

z1 = deseq_z(pb_sf.loc[m1.values], cd_sf.loc[m1.values])
z2 = deseq_z(pb_sf.loc[m2.values], cd_sf.loc[m2.values])

sig = pd.DataFrame({
    "z_discovery_tumor": res_sf["zscore"],
    "log2FC_tumor": res_sf["log2FoldChange"],
    "lfcSE_tumor": res_sf["lfcSE"],
    "padj_tumor": res_sf["padj"],
})
sig["z_tumor_splitA"] = z1
sig["z_tumor_splitB"] = z2
sig["z_blood_bladder"] = res_bl["zscore"]
sig.index.name = "gene"
sig = sig.sort_values("z_discovery_tumor", ascending=False)
sig.to_csv("CD4_RvsNR_signature.csv")
print("signature genes:", len(sig))