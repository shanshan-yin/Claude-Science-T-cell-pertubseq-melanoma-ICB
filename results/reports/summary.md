# Mapping the CD4⁺ T-cell ICB-resistance state to perturb-seq regulators

*Using scRNAseq and Perturb-seq to investigate the novel regulatory pathways in CD4+ T cells to change melanoma ICB non-responding state.*

## Question

**Which CD4⁺ knockdowns reproduce the transcriptional state that distinguishes ICB
non-responders from responders — and are any of them druggable?**

## Method (mirrors the source notebooks)

1. **Signature (Panel B).** CD4⁺ T cells isolated from an ICB scRNA-seq cohort, pseudobulked
   per patient, DESeq2 non-responder-vs-responder → per-gene z-score = log₂FC / lfcSE
   (identical to the paper's definition).
2. **Perturbation matrix.** The genome-scale CD4⁺ perturb-seq atlas (`GWCD4i.DE_stats`),
   subset to the 7,874 knockdown×condition profiles passing the ≥10-DE-gene filter,
   zscore layer (log_fc/lfcSE), cis/on-target diagonal effects masked and mean-imputed.
3. **Model.** `Perturb2StateModel` (ElasticNet + TruncatedSVD PCA, 5-fold CV, n_pcs=60,
   scale_y, seed=214) regressing the signature (per-gene z, target) onto the perturbation
   matrix (features), per culture condition (Rest / Stim 8 hr / Stim 48 hr).

## Data

| Cohort | Accession | Compartment | Platform | CD4 cells | Samples (R/NR) |
|---|---|---|---|---|---|
| Sade-Feldman melanoma | GSE120575 | tumor | Smart-seq2 TPM | 5,391 | 43 (16 R / 27 NR) |

Yost BCC/SCC (GSE123813) was evaluated but deferred: it has clean CD4 subset clusters but
its R/NR labels live only in a supplementary table behind access that was not reachable this
session.

## Key findings

**1. The melanoma-tumor CD4 resistance signature recapitulates published ICB-resistance biology (Panel B).**
The CD38⁺CD39⁺ T-eee exhaustion program (mean z = +1.7; ENTPD1/CD39, CD38, PD-1, LAG3,
TIM-3) and a cytotoxic-CD4 program (mean z = +2.0; GZMB, PRF1, CX3CR1, IFNG) are all up in
non-responders, while responders show the LEF1/TCF7/FOXP1 stem-memory program. In the
volcano, 218 genes are up in non-responders and 91 up in responders.

**2. The perturb-seq atlas predicts the tumor signature strongly (Panels C, D).**
Held-out cross-validation Pearson r (observed vs predicted):

| | Rest | Stim 8 hr | Stim 48 hr |
|---|---|---|---|
| **Tumor CD4** | **0.574** | 0.512 | 0.506 |

All models exceed the within-tumor discovery-vs-replication upper bound (patient split,
r = 0.285) — the healthy-donor CD4 perturbation atlas contains regulators whose knockdown
reproduces the tumor resistance state better than one melanoma half-cohort predicts the other.

**3. Two opposing pathways drive the axis (Panels D–F).**
Two pathways dominate the model in opposite directions.

- **Reproduce the resistance state** (knockdown → non-responder-like): ELOF1, and the SAGA
  deubiquitinase pathway **ATXN7L3 / USP22 / SGF29**, plus LIN37 (DREAM complex), FUBP1,
  YTHDF2. Signed knockdown analysis shows the SAGA pathway natively restrains the
  proliferation program — losing it raises the resistance genes (Panel E).
- **Reproduce the response state** (knockdown → responder-like), top-10 by most-negative
  coefficient: WAC, **IRF9**, WDR82, UBR5, **IRF1**, MED24, EPC1, KMT5B, DOP1A, and CARMIL2
  (a known regulator of CD4 costimulation) — the interferon regulators IRF9/IRF1 rank
  high; IFNAR1 (rank 12) and BRD8 (rank 13) fall just outside the top 10. Native interferon
  signaling activates the resistance/proliferation program, so inhibiting this axis pushes
  cells toward the responder state (Panel F).

**4. The responder-reproducing axis is druggable with approved compounds (Panel G).**
Because a small-molecule inhibitor mimics a knockdown, the response-reproducing regulators
are the ones worth inhibiting to sensitize non-responders. Screening these against ChEMBL,
the druggable candidates all fall on the interferon / JAK-STAT axis: IFNAR1 (anifrolumab, an
approved antagonist antibody), TYK2 (approved JAK-family inhibitors) and JAK2 (12 approved
inhibitors including ruxolitinib, momelotinib, baricitinib). The SAGA pathway, by contrast,
is both the wrong therapeutic direction and essentially undruggable.

## Limitations

- **Single tumor cohort.** The tumor signature rests on one melanoma dataset (43 samples);
  replication is a within-cohort patient split, not a second independent tumor cohort. A
  second tumor cohort would strengthen the Panel C ceiling.
- **Healthy-donor atlas vs TME.** The perturb-seq atlas is from in-vitro-stimulated
  healthy-donor CD4 cells; it cannot contain microenvironmental signals (hypoxia, adenosine,
  chronic antigen) that shape the real resistance state. Nominated regulators are hypotheses
  about cell-intrinsic transcriptional control, not validated TME drivers.
- **Drug nomination is directional, not preclinical.** The ChEMBL matches identify approved
  compounds against the right targets and direction; they are a starting point for
  combination hypotheses, not evidence of efficacy.
- **Published-program provenance.** The CD38/CD39 T-eee gene set was verified in the source
  paper's full text (Mitra & Thompson 2023, CCR-23-0653). The cytotoxic-CD4 gene set is a
  domain-curated canonical CTL marker list (concept per Bae 2025); its genes were **not**
  extracted from that paper's tables.

## Files

- `figure5_composite.png` — 8-panel composite (A schematic, B volcano, C accuracy,
  D ranked regulators, E regulator effects, F SAGA volcano, G interferon volcano, H drug targets)
- `CD4_RvsNR_signature.csv` — tumor signature + discovery + split-replication z-scores
- `pert2state_models.pkl` — fitted models; `prediction_metrics.csv` — accuracy table
- `ranked_regulators.csv` — ElasticNet coefficients (tumor), all regulators
- `published_CD4_programs.csv`, `cd4_cells_metadata.csv`, `acquisition_manifest.csv`
