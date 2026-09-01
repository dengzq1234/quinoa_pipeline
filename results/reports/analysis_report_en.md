# Quinoa Salt Tolerance — Analysis Report

**Date:** 2026-09-01  
**Pipeline:** quinoa_pipeline (Steps 1–5, all R)  
**Core question:** Which genes respond differently to NaCl in salt-tolerant (ST) vs
non-tolerant (NT) quinoa accessions, and what are their biological functions?

---

## 0. Summary table

| # | Question | Analysis used | One-line answer | Key output |
|---|---|---|---|---|
| **M** | Is rentrez + blastp the right approach for functional annotation? | Replaced with biomaRt (Ensembl Plants) curated orthology | Better approach: biomaRt maps quinoa → Arabidopsis orthologs directly; no FASTA or local BLAST needed | `pipeline_update_summary.txt` |
| **Core** | Which genes differ between ST and NT under NaCl? What do they do? | DESeq2 interaction LRT (Step 1) + GO enrichment (Steps 3 / 5) | 2,551 DEGs; functions centred on **photosynthesis / chloroplast thylakoid / pigment (chlorophyll) metabolism / cell wall**; ST maintains and upregulates the photosynthetic apparatus; NT collapses | `degs_salt_tolerance/`, `05_cluster_go/` |
| **Q1** | Do three qPCR validation genes (PPR40 / AS-1 / trmH) respond to NaCl? | Simplified model `~ Tol + Treatment`, NaCl main effect (T30) | All three are **non-significant** (padj 0.73 / 0.91 / 0.90); \|log2FC\| < 0.3 | `04_cluster/Q1_qPCR_NaCl_main_effect.tsv` |
| **Q2** | Can differentially responding genes be grouped by expression profile? | `degPatterns` (DEGreport), top 1,000 DEGs | **4 clusters**, 947 genes; two main axes = ST-specific induction (C2, n=495) and NT collapse (C3, n=178) | `04_cluster/Q2_*` |
| **Q3** | Are there outlier samples (some replicates cluster poorly)? | Blind VST PCA across all samples; flag > 2 SD from group centroid (Tol × Treatment × Time) | **No outliers detected** (0 samples flagged) | `04_cluster/Q3_outlier_samples.tsv` |

---

## 1. Data and methods

- **Input:** `gene_count.csv` (58,884 genes × 73 samples, raw read counts) +
  `samples.csv` (`Treatment` = CONT/NaCl; `Tol` = ST/NT; `Time` = T0/T30).
- **Low-expression filter:** `rowSums(counts) ≥ 10` → ~44,500 expressed genes.
- **Grouping variable:** `Tol` (ST/NT biological class) — **not** `Tolerance`
  (accession IDs 14/227/287/460/577/670).
- **Differential expression model (Step 1, T30 samples only):**
  - Full model: `~ Tol + Treatment + Tol:Treatment`
  - LRT against reduced model: `~ Tol + Treatment`
  - Contrast extracted: `TolST.TreatmentNaCl`
    (positive = more strongly NaCl-induced in ST than NT)
  - Significance: `padj < 0.05` and `|log2FC| ≥ 1`
- **Functional annotation (Step 3):** quinoa NCBI gene ID → quinoa Ensembl ID →
  Arabidopsis ortholog via Ensembl Plants biomaRt (two-query approach to avoid
  the attribute-page error).
- **GO enrichment:** `clusterProfiler` + `org.At.tair.db`, all three ontologies
  (BP / MF / CC), BH correction.
- **Clustering (Step 4):** `degPatterns` (DEGreport) on top 1,000 DEGs by padj;
  `time = Treatment`, `col = Tol`.
- **Cluster GO enrichment (Step 5):** `enrichGO` on Arabidopsis orthologs per
  cluster; **universe = all expressed-gene orthologs (~8,218 Arabidopsis genes)**,
  not just DEGs — this is the statistically correct background.

---

## 2. Results

### 2.1 Annotation approach (M)

The original Python/bash pipeline used `rentrez` to fetch sequences and local
`blastp`/DIAMOND against Arabidopsis. **Replaced with biomaRt (Ensembl Plants)
curated ortholog mapping:**

- Direct quinoa → Arabidopsis gene pairs; no FASTA download, no local BLAST.
- `clusterProfiler` runs BP/MF/CC enrichment in one step.
- FASTA retrieval and DIAMOND BLAST retained as optional scripts
  (`scripts/optional/`); not part of the main pipeline.
- Original Python scripts archived under `scripts/legacy/`.

> A test run on 1,411 uncharacterised quinoa genes validated end-to-end
> feasibility (`gene_list_1411/`), recovering biologically plausible terms
> (photosynthesis, chloroplast thylakoid, glycoside hydrolase activity).
> **Note: this test used a separate gene list, not the DEG list**, and serves
> only as a methodological proof-of-concept. The current pipeline uses the
> full 2,551-gene DEG list.

### 2.2 Differential expression (Step 1)

| Metric | Value |
|---|---|
| Significant DEGs (padj < 0.05, \|log2FC\| ≥ 1) | **2,551** |
| More strongly NaCl-induced in ST (up in ST) | 2,337 |
| More strongly NaCl-suppressed in ST (down in ST) | 214 |

The strong directional asymmetry (2,337 vs 214) indicates that **salt tolerance
in ST accessions is driven primarily by gene activation, not repression.**

Outputs (in `rnaseq_analysis/results/`): full results table, PCA, volcano plot,
top-500 heatmap; DEG ID list `03_DEG_list_for_annotation.txt` used as input for
Steps 3 and 4.

### 2.3 DEG-level GO enrichment (Step 3)

637 of 2,551 DEGs (~25%) mapped to an Arabidopsis ortholog. GO enrichment
(119 significant terms; `degs_salt_tolerance/04_GO_enrichment_results.tsv`)
shows four main themes:

- **Cell wall / polysaccharide / xylan:** `polysaccharide metabolic process`,
  `plant-type secondary cell wall biogenesis`, `xylan metabolic process`,
  `apoplast` (CC, p.adj ≈ 5×10⁻⁹) — CESA4/7, IRX9/15, XTH, GUX2, FLA11.
- **Pigment / chlorophyll metabolism:** `pigment metabolic process`
  (p.adj ≈ 2×10⁻⁵), `porphyrin` — CHLM, FLU, ACSF, LIL3, BPG1.
- **Fatty acid / wax / phenylpropanoid:** `fatty acid biosynthetic process`
  (p.adj ≈ 2×10⁻⁶), `very-long-chain fatty acid`, cuticle/epicuticular wax
  — CER4/6/26, KCS, LACS2.
- **Chloroplast organisation / protein targeting:** `plastid organization`
  (p.adj ≈ 9×10⁻⁸), `protein localization to chloroplast`.
- **Stress response:** `response to water deprivation`, `response to cold`,
  `antioxidant activity`, `detoxification` — GST, peroxidase, SOD/FSD3.

### 2.4 Q1 — NaCl main effect for qPCR validation genes

Model `~ Tol + Treatment` (T30), contrast `Treatment_NaCl_vs_CONT`:

| gene_id | Gene | baseMean | log2FC | pvalue | padj | Significant |
|---|---|---|---|---|---|---|
| 110697655 | PPR40 | 120.4 | −0.275 | 0.589 | 0.732 | No |
| 110710750 | AS-1 | 241.2 | +0.110 | 0.833 | 0.905 | No |
| 110703716 | trmH | 63.4 | +0.175 | 0.818 | 0.895 | No |

**Conclusion:** At T30, pooling ST and NT, none of the three genes show a
significant transcriptional NaCl response; effect sizes are negligible
(\|log2FC\| < 0.3). Their stable expression across conditions makes them
candidates for RT-qPCR reference/housekeeping genes.

If qPCR data show changes, possible explanations are: (a) accession-specific
effects (interaction term, not main effect); (b) time-point differences
(T0 vs T30 or earlier); (c) post-transcriptional regulation.

### 2.5 Q2 — Expression pattern clustering

`degPatterns` on top 1,000 DEGs; 947 genes assigned to **4 clusters**:

| Cluster | n genes | Expression pattern |
|---|---|---|
| **1** | 257 | Crossover: ST↑ / NT↓ under NaCl |
| **2** | 495 | ST-specific strong induction (NT largely unchanged) |
| **3** | 178 | NT collapses under NaCl; ST remains stable |
| **4** | 17 | ST-specific downregulation (NT stable) |

C2 + C3 together account for 71% of clustered genes and represent opposite
poles of the same biological axis: **ST active induction / NT passive loss.**

Outputs: `04_cluster/Q2_degPatterns_clusters.png`,
`04_cluster/Q2_cluster_assignments.tsv`.

### 2.6 Q3 — Outlier sample detection

Blind VST PCA across all 73 samples; samples flagged if > 2 SD from their
group centroid (groups defined by Tol × Treatment × Time).

**Result: 0 outlier samples** (`Q3_outlier_samples.tsv` is empty). Data quality
is good; no sample removal is needed.

The visual impression of some replicates clustering loosely is expected —
these are field-collected accessions, not isogenic lines, so within-group
biological variance is real. If a reviewer raises this concern, a sensitivity
analysis (re-run DESeq2 excluding the visually atypical replicates and show
results are stable) is the strongest response.

Output: `04_cluster/Q3_PCA_outliers.png`.

### 2.7 Per-cluster GO enrichment (Step 5, genome-wide universe)

Universe = ~8,218 Arabidopsis orthologs of all expressed quinoa genes
(statistically correct; using only DEG orthologs as background would
artificially inflate enrichment).

| Cluster | Sig. terms | Top finding | Biological interpretation |
|---|---|---|---|
| **1** (crossover ST↑/NT↓) | 2 | `phototropism` (p.adj 0.019; RPT2/JK218/JK224); `response to lipid` (p.adj 0.043; LOX2/LTP3/annexin/GASA) | Divergence in light signalling and membrane/lipid signalling between ST and NT |
| **2** (ST-specific induction) | 40 | **`thylakoid` (CC, p.adj 4×10⁻³³)**; `photosynthesis` (BP, p.adj 2×10⁻¹⁷); `light harvesting`; `photosynthetic membrane`; `plastid stroma/envelope` | **ST actively upregulates chloroplast light-harvesting complexes and photosynthetic electron transport under NaCl** (LHCA1/2/3, LHCB2/3, CAB4, CP24, PSB\*, PSA\*, PETC, RBCS1A, SBPASE, PRK) |
| **3** (NT collapse) | 23 | `pigment metabolic process` (p.adj 3×10⁻⁵); `tetrapyrrole / porphyrin metabolic process`; `chlorophyll metabolic process`; `lipid / wax biosynthetic process`; `cell wall` (CC) | **NT loses chlorophyll/pigment biosynthesis capacity and cell wall integrity under NaCl** (AtCLA1, ACLA-1, GBP, PSB29, CER6, KCS, XTH9, GASA1, BGAL3) |
| **4** (ST downregulation) | — | < 5 orthologs; enrichment skipped | Insufficient gene count for interpretation |

---

## 3. Biological conclusions

1. **The core of salt tolerance is maintaining photosynthesis.** All major signals
   converge on the chloroplast: ST upregulates thylakoid light-harvesting complexes
   and the photosynthetic electron transport chain (Cluster 2); NT loses
   chlorophyll/tetrapyrrole biosynthesis capacity (Cluster 3). These are the two
   poles of a single physiological axis.

2. **ST tolerance is additive, not subtractive.** 2,337 genes are more strongly
   induced in ST vs only 214 more suppressed; ST-specific induction (C2, 495 genes)
   vastly outnumbers ST-specific repression (C4, 17 genes).

3. **A secondary structural axis: cell wall and cuticle.** Both DEG-level and
   Cluster 3 enrichments recover secondary cell wall biosynthesis (CESA/IRX/XTH)
   and very-long-chain fatty acid/wax genes (CER/KCS) — NT may simultaneously
   lose cell wall and epidermal barrier integrity under salt stress.

4. **The three qPCR genes (PPR40 / AS-1 / trmH) are not part of this main axis**
   (see Q1). For validation experiments, LHC genes from Cluster 2 and
   chlorophyll-synthesis genes from Cluster 3 are more compelling candidates.

---

## 4. Limitations and next steps

- **Annotation coverage:** only ~19% of expressed quinoa genes and ~25% of DEGs
  map to an Arabidopsis ortholog. GO results are biased toward conserved,
  well-characterised pathways; quinoa-specific or uncharacterised genes are
  systematically missed.
- **Q1 reports NaCl main effect only:** if the interaction term
  (`TolST.TreatmentNaCl`) for the three qPCR genes is needed, it can be
  extracted separately from the full DESeq2 results table.
- **Pipeline integration:** Steps 4, 4a, 5a, and 5 are currently run independently
  of `run_pipeline.R`.

### Suggested next steps

1. Sensitivity analysis for Q3: re-run DESeq2 excluding visually atypical
   replicates; confirm DEG overlap and GO stability.
2. Targeted visualisation of Cluster 2 / 3 genes (LHC and chlorophyll-synthesis
   genes): normalised count boxplots (ST/NT × CONT/NaCl) as main-figure candidates.
3. Consider eggNOG-mapper or the quinoa genome annotation (if available) to
   supplement ortholog coverage and recover missing genes.
4. Incorporate Steps 4/5/5a into `run_pipeline.R` and convert all paths to
   relative paths for portability.

---

## 5. File index

```
results/
├── degs_salt_tolerance/
│   ├── 04_quinoa_arabidopsis_orthologs_2551.tsv   # 637 quinoa→Arabidopsis pairs
│   ├── 04_GO_enrichment_results.tsv               # DEG-level GO (119 terms)
│   └── 04_GO_dotplot.png
├── 04_cluster/
│   ├── Q1_qPCR_NaCl_main_effect.tsv               # PPR40 / AS-1 / trmH
│   ├── Q2_degPatterns_clusters.png                # 4-cluster expression pattern plot
│   ├── Q2_cluster_assignments.tsv                 # 947 genes × cluster
│   ├── Q3_PCA_outliers.png
│   └── Q3_outlier_samples.tsv                     # empty = no outliers
└── 05_cluster_go/
    ├── background_orthologs_all.tsv               # GO universe (~8,218 Arabidopsis genes)
    ├── cluster1_GO_corrected.tsv / _dotplot.png   # 2 terms
    ├── cluster2_GO_corrected.tsv / _dotplot.png   # 40 terms (photosynthesis / thylakoid)
    ├── cluster3_GO_corrected.tsv / _dotplot.png   # 23 terms (pigment / cell wall)
    └── all_clusters_GO_corrected_combined.tsv     # 65 terms combined
```
