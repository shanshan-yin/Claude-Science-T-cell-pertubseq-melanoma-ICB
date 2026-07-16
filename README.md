# Regulatory pathways of the CD4⁺ T-cell melanoma ICB non-responding state

Using single-cell RNA-seq and genome-scale Perturb-seq to find the regulatory
pathways in CD4⁺ T cells that drive the melanoma immune-checkpoint-blockade (ICB)
**non-responding** state — and to ask whether those pathways are druggable.

Immune checkpoint blockade reactivates T cells to kill tumors and has transformed
cancer treatment, but many patients don't respond. Most work on ICB response has
focused on CD8⁺ killer T cells; emerging evidence shows the CD4⁺ helper T cell is
also important. Focusing on CD4⁺ helper cells, we compare responders and
non-responders in scRNA-seq of 32 melanoma patients to define a non-responder
signature, then integrate a genome-scale Perturb-seq atlas of ~2,600 gene knockdowns
in healthy-donor CD4⁺ T cells. Because each knockdown is a **direct perturbation**,
the model learns which regulators *causally* push cells toward responding or
non-responding — not merely what correlates.

## Explore the results online

Interactive, self-contained pages (no install required):

- **[Landing page](https://raw.githack.com/shanshan-yin/Claude-Science-T-cell-pertubseq-melanoma-ICB/main/websites/index.html)** — overview and links
- **[Full report](https://raw.githack.com/shanshan-yin/Claude-Science-T-cell-pertubseq-melanoma-ICB/main/websites/report.html)** — the 8-panel composite figure with narrative:
  signature, prediction accuracy, opposing regulatory pathways, and druggable targets
- **[Gene explorer](https://raw.githack.com/shanshan-yin/Claude-Science-T-cell-pertubseq-melanoma-ICB/main/websites/gene_explorer.html)** — look up any gene: its direction in
  the resistance signature, and whether it acts as a regulator or downstream target of
  a pathway


## Summary of findings

1. **By GO/pathway enrichment, the resistance signature is dominated by cell-cycle and interferon programs, and the response signature by a Wnt/memory program.**
   Genes up in non-responders are most strongly enriched for cell-cycle and
   proliferation (E2F targets, G2-M checkpoint, mitotic cell cycle; adjusted
   p ≈ 7×10⁻³⁹), with a secondary interferon-γ signaling signal (cellular
   response to IFN-γ; adjusted p ≈ 2×10⁻⁵). Genes up in responders are enriched
   for Wnt/β-catenin–TCF/LEF signaling — the stem/memory program. In the
   volcano, 218 genes are up in non-responders and 91 up in responders.
2. **The healthy-donor Perturb-seq atlas predicts the tumor signature strongly.**
   Held-out cross-validation Pearson r = **0.57** (Rest), 0.51 (Stim 8 hr),
   0.51 (Stim 48 hr) — all above the within-tumor replication ceiling (r = 0.29).
   A model trained on the resting-condition data reproduces the patient signature
   well, and the top knockdowns cluster into coherent pathways.
3. **Two opposing pathways drive the axis.** The **interferon pathway** (IRF1/IRF9,
   IFNAR1) reproduces the *responder* state on knockdown, and its downstream genes are
   up in non-responders — so inhibiting interferon should push cells toward responding.
   The **SAGA/chromatin-DUB pathway** (USP22/ATXN7L3/SGF29) reproduces *resistance*.
   The two converge on the same shared cell-cycle/proliferation genes with opposite
   regulatory sign.
4. **The responder-reproducing axis is druggable.** The interferon/JAK-STAT targets
   map to approved compounds (anifrolumab for IFNAR1; JAK/TYK2 inhibitors such as
   ruxolitinib, baricitinib, upadacitinib) — turning the analysis into a concrete,
   testable hypothesis to sensitize non-responders.

See [`results/reports/summary.md`](results/reports/summary.md) for the full write-up
and [`results/reports/CD4_Perturbseq_Melanoma_ICB_Summary.docx`](results/reports/CD4_Perturbseq_Melanoma_ICB_Summary.docx)
for the formatted report.

## Data

Two datasets underlie the analysis — a **clinical signature** cohort (scRNA-seq) and a
**perturbation atlas** (genome-scale Perturb-seq):

| Dataset | Accession | Compartment | Platform | Size | Role |
|---|---|---|---|---|---|
| Sade-Feldman melanoma | [GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575) | tumor CD4⁺ | Smart-seq2 TPM | 5,391 cells · 43 samples (16 R / 27 NR) | R-vs-NR clinical signature |
| GWCD4i Perturb-seq atlas | CZI (`GWCD4i.DE_stats`) | healthy-donor CD4⁺ | CRISPRi Perturb-seq | ~2,600 knockdowns × 3 conditions · 7,874 profiles | causal perturbation matrix |

The Sade-Feldman melanoma ICB cohort is public via GEO. The genome-scale CD4⁺
Perturb-seq atlas is a large public resource and is **not** redistributed here; the
scripts stream the required z-score layer directly. Raw inputs are excluded via
`.gitignore`.

## Repository layout

```
.
├── code/                               reproduction scripts (numbered by pipeline stage)
│   ├── 01_build_CD4_RvsNR_signature.py      per-cohort pseudobulk DESeq2 → per-gene z
│   ├── 02_fit_perturb2state_rank_regulators.py  fit Perturb2StateModel, rank regulators
│   ├── 03_prediction_metrics.py             held-out accuracy per condition
│   ├── 04_panelB_volcano.py                 Panel B — R-vs-NR volcano
│   ├── 05_panelC_prediction.py              Panel C — prediction accuracy bars
│   ├── 06_panelD_regulator_effects.py       Panel E — predicted regulator effects
│   ├── 07_panelEF_pathway_regulated_genes.py  Panels F/G — pathway-regulated genes
│   └── 10_figure5_composite.py              assemble the 8-panel composite (A–H)
├── results/
│   ├── figures/    figure5_composite.png + individual panels (A, B, C, D, E, F, G, H)
│   ├── tables/     signature, ranked regulators, prediction metrics, gene sets, metadata
│   └── reports/    summary.md + Word report
├── websites/       self-contained interactive pages (landing, report, gene explorer)
├── requirements.txt
└── README.md
```

## Reproducing

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install https://github.com/emdann/pert2state_model/archive/refs/heads/main.tar.gz
# scripts in code/ are numbered in pipeline order; run 01 → 10
```

The scripts were produced in a Python 3.11 conda environment; exact package versions
are pinned in `requirements.txt`.

## Method (brief)

- **Signature.** CD4⁺ cells isolated from the cohort, pseudobulked per patient, DESeq2
  non-responder-vs-responder; per-gene z-score = log₂FC / lfcSE.
- **Model.** `Perturb2StateModel` (ElasticNet + TruncatedSVD PCA, 5-fold CV) regressing
  the signature onto the perturbation matrix, per culture condition (Rest / Stim 8 hr /
  Stim 48 hr).

---

*Analysis and figures produced with Claude Science. Drug-target annotations were
retrieved from ChEMBL; gene-set enrichment used Enrichr.*
