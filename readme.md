# Quinoa Salt Tolerance Pipeline

End-to-end RNA-seq analysis pipeline for identifying biological functions that
distinguish salt-tolerant (ST) from non-tolerant (NT) quinoa accessions under
NaCl stress.

**Biological question:** Which genes respond differently to NaCl in ST vs NT
accessions, and what do those genes do?

---

## Requirements

Run once to install all required R packages:

```r
source("scripts/00_install_packages.R")
```

---

## Project Structure

```
├── data/
│   ├── input/              # standalone gene ID lists for annotation-only runs
│   └── reference/          # Arabidopsis DIAMOND database (optional BLAST step)
├── results/
│   ├── 01_deseq2/          # DESeq2 output (PCA, volcano, heatmap, DEG list)
│   ├── 02_filter/          # filtered uncharacterized genes (optional)
│   ├── 03_annotation/      # orthologs + GO enrichment (DEG-level)
│   ├── 04_cluster/         # expression clustering (degPatterns) + outlier QC
│   └── 05_cluster_go/      # GO enrichment per expression cluster
├── scripts/
│   ├── 00_install_packages.R       # one-time setup
│   ├── 01_deseq2_analysis.R        # Step 1: differential expression ST vs NT
│   ├── 02_filter_gene_list.R       # Step 2 (optional): filter uncharacterized genes
│   ├── 03_annotate_and_go.R        # Step 3: biomaRt orthology + GO enrichment
│   ├── 04_cluster_analysis.R       # Step 4: expression clustering + outlier QC
│   ├── 05a_build_background_orthologs.R  # Step 5a: full background for GO
│   ├── 05_cluster_go_enrichment.R  # Step 5: GO enrichment per cluster
│   ├── optional/
│   │   ├── get_fasta.R             # fetch protein FASTAs from NCBI
│   │   └── blast_against_arabidopsis.sh  # DIAMOND BLAST vs Arabidopsis
│   └── legacy/                     # original Python pipeline (archived)
└── run_pipeline.R                  # entry point for Steps 1–3
```

---

## Running the full pipeline

Edit the two paths at the top of `run_pipeline.R`:

```r
rnaseq_data_dir <- "/path/to/rnaseq_analysis/data"   # gene_count.csv + samples.csv
results_base    <- "/path/to/quinoa_pipeline/results"
```

Then run from the `quinoa_pipeline/` directory:

```r
source("run_pipeline.R")          # Steps 1–3
source("scripts/04_cluster_analysis.R")        # Step 4 (standalone)
source("scripts/05a_build_background_orthologs.R")  # Step 5a (run once)
source("scripts/05_cluster_go_enrichment.R")   # Step 5 (standalone)
```

Steps 4, 5a, and 5 can also be run standalone — each has `if (!exists(...))` guards
for their input paths so they work independently or sourced from a master script.

---

## Pipeline steps

### Step 1 — DESeq2 differential expression (`01_deseq2_analysis.R`)

**Input:** `gene_count.csv` (58,884 genes × 73 samples), `samples.csv`

**What it does:**
- Filters to T30 samples (post-stress phase)
- Model: `~ Tol + Treatment + Tol:Treatment` with LRT
- Tests which genes respond differently to NaCl in ST vs NT accessions
- Uses `Tol` (ST/NT biological class), NOT `Tolerance` (accession ID)

**Outputs in `results/01_deseq2/`:**
- `01_PCA_all_samples.png` — QC PCA of all 73 samples
- `02_DESeq2_T30_TolxTreatment_full.tsv` — full results table with gene annotations
- `03_DEG_list_for_annotation.txt` — DEG gene IDs (input for Steps 3 and 4)
- `03a_DEGs_up_in_ST.tsv` — 2,337 genes more induced by NaCl in ST than NT
- `03b_DEGs_down_in_ST.tsv` — 214 genes more suppressed by NaCl in ST than NT
- `04_PCA_T30.png` — PCA of T30 samples only
- `05_volcano_T30_ST_vs_NT.png` — volcano plot
- `06_heatmap_top500.png` — heatmap of top 500 DEGs

---

### Step 2 — Filter uncharacterized genes (`02_filter_gene_list.R`) — optional

**Input:** gene ID list (one NCBI gene ID per line)

**What it does:** Queries NCBI via `rentrez` to fetch gene metadata and keeps
only uncharacterized protein-coding genes. Skipped by default in `run_pipeline.R`
because GO enrichment works directly on the full DEG list.

> Note: requires ~17 min for 2,500 genes due to NCBI API rate limits (0.4 s/gene).

**Output:** TSV with columns `gene_id, gene_name, gene_desc, gene_type,
refseq_genomic, refseq_mrna, refseq_peptide, related_acc`

---

### Step 3 — Arabidopsis ortholog mapping + GO enrichment (`03_annotate_and_go.R`)

**Input:** DEG gene ID list from Step 1 (or Step 2 if filtered)

**What it does:**
- Maps quinoa NCBI gene IDs → quinoa Ensembl IDs → Arabidopsis orthologs
  via Ensembl Plants (biomaRt, two-query approach to avoid attribute-page error)
- Runs GO enrichment across BP, MF, CC ontologies using `clusterProfiler`
  and `org.At.tair.db`

**Outputs in `results/03_annotation/`:**
- `04_quinoa_arabidopsis_orthologs.tsv` — quinoa → Arabidopsis gene mapping (637 pairs)
- `04_GO_enrichment_results.tsv` — full GO enrichment table
- `04_GO_dotplot.png` — dotplot of top 15 GO terms per ontology

---

### Step 4 — Expression clustering + outlier QC (`04_cluster_analysis.R`)

**Input:** `gene_count.csv`, `samples.csv`, DEG list from Step 1

**What it does (answers Raquel's three questions):**
- **Q1** — NaCl main effect (model: `~ Tol + Treatment`) for three qPCR validation
  genes: PPR40 (110697655), AS-1 (110710750), trmH (110703716)
- **Q2** — Expression clustering via `degPatterns` (DEGreport) on the top 1,000
  DEGs by padj. Identifies groups of genes with shared ST/NT × NaCl profiles.
- **Q3** — Outlier sample detection: PCA of all samples, flag any sample
  > 2 SD from its group centroid (Tol × Treatment × Time)

**Outputs in `results/04_cluster/`:**
- `Q1_qPCR_NaCl_main_effect.tsv` — log2FC, p-value, padj for the three qPCR genes
- `Q2_degPatterns_clusters.png` — expression profile plot (4 clusters found)
- `Q2_cluster_assignments.tsv` — per-gene cluster assignment (947 genes × 4 clusters)
- `Q3_PCA_outliers.png` — PCA with outlier samples circled
- `Q3_outlier_samples.tsv` — list of flagged samples

**Cluster summary (947 genes total):**

| Cluster | n genes | Expression pattern |
|---------|---------|-------------------|
| 1 | 257 | Crossover: ST↑ / NT↓ under NaCl |
| 2 | 495 | ST-specific strong induction under NaCl |
| 3 | 178 | NT collapses under NaCl; ST stable |
| 4 | 17  | ST-specific downregulation under NaCl |

---

### Step 5a — Build background ortholog table (`05a_build_background_orthologs.R`)

**Run once before Step 5.** Output is cached; re-running is a no-op if the file exists.

**What it does:**
- Downloads the full quinoa gene list from Ensembl Plants (biomaRt, no filter)
- Filters to all expressed genes (rowSums ≥ 10, n = 44,524)
- Downloads the complete quinoa → Arabidopsis ortholog table (no filter)
- Joins locally to produce a correct GO enrichment universe

> **Technical note:** biomaRt 2.66 uses an httr2 POST backend incompatible with
> Ensembl Plants (HTTP 405). This script patches `.submitQueryXML` to use GET
> and downloads full tables without row-level filters, then filters locally.

**Output:** `results/05_cluster_go/background_orthologs_all.tsv`  
8,401 expressed quinoa genes with Arabidopsis ortholog → 8,218 unique Arabidopsis
genes used as GO enrichment universe.

---

### Step 5 — GO enrichment per cluster (`05_cluster_go_enrichment.R`)

**Input:**
- Cluster assignments from Step 4 (`Q2_cluster_assignments.tsv`)
- DEG ortholog table from Step 3 (`04_quinoa_arabidopsis_orthologs_2551.tsv`)
- Background universe from Step 5a (`background_orthologs_all.tsv`)

**What it does:**
- For each cluster, retrieves the Arabidopsis ortholog genes
- Runs `enrichGO` (clusterProfiler, org.At.tair.db, keyType="TAIR", ont="ALL")
  using the full expressed-gene ortholog set as universe (not just DEGs)
- Saves per-cluster TSV and dotplot

**Outputs in `results/05_cluster_go/`:**
- `cluster1_GO_corrected.tsv` / `_dotplot.png` — 2 significant terms
- `cluster2_GO_corrected.tsv` / `_dotplot.png` — 40 significant terms
- `cluster3_GO_corrected.tsv` / `_dotplot.png` — 23 significant terms
- `all_clusters_GO_corrected_combined.tsv` — all 65 terms combined

**Key biological findings:**

| Cluster | Top GO finding | Interpretation |
|---------|---------------|----------------|
| 2 | Thylakoid (CC, p.adj = 4×10⁻³³); Photosynthesis (BP, p.adj = 2×10⁻¹⁷) | ST actively upregulates chloroplast light-harvesting machinery under NaCl |
| 3 | Pigment metabolic process; Porphyrin/tetrapyrrole metabolism; Cell wall | NT loses chlorophyll biosynthesis capacity and cell wall integrity under NaCl |
| 1 | Phototropism; Response to lipid | Divergent light and membrane signaling between ST and NT |

Clusters 2 and 3 together describe opposite poles of the same biological axis:
ST protects/induces photosynthetic capacity while NT simultaneously fails to
maintain pigment synthesis and structural integrity.

---

## Optional steps

Only needed if protein FASTA sequences are required (e.g. phylogenetics,
structural prediction):

```r
source("scripts/optional/get_fasta.R")
# bash scripts/optional/blast_against_arabidopsis.sh
```

---

## Notes

- `Tol` (ST/NT) is the correct biological grouping variable — not `Tolerance`
  (accession numbers 14, 227, 287, 460, 577, 670)
- GO enrichment uses Arabidopsis annotations since quinoa lacks a comprehensive
  GO database; ~19% of expressed quinoa genes have a mapped Arabidopsis ortholog
- Steps 4 and 5 are run independently of `run_pipeline.R` (no integration yet)
- Orthology uses curated Ensembl Plants mappings, not BLAST hits
- The biomaRt GET patch in Step 5a is specific to biomaRt ≥ 2.66 on
  Ensembl Plants; earlier versions or other marts do not need it
