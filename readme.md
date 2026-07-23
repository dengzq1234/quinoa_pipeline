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
│   ├── 02_filter/          # filtered uncharacterized genes (optional step)
│   └── 03_annotation/      # orthologs + GO enrichment
├── scripts/
│   ├── 00_install_packages.R       # one-time setup (all three steps)
│   ├── 01_deseq2_analysis.R        # Step 1: differential expression ST vs NT
│   ├── 02_filter_gene_list.R       # Step 2 (optional): filter uncharacterized genes
│   ├── 03_annotate_and_go.R        # Step 3: biomaRt orthology + GO enrichment
│   ├── optional/
│   │   ├── get_fasta.R             # fetch protein FASTAs from NCBI
│   │   └── blast_against_arabidopsis.sh  # DIAMOND BLAST vs Arabidopsis
│   └── legacy/                     # original Python pipeline (archived)
└── run_pipeline.R                  # single entry point for the full pipeline
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
source("run_pipeline.R")
```

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
- `03_DEG_list_for_annotation.txt` — DEG gene IDs (input for Step 3)
- `03a_DEGs_up_in_ST.tsv` — genes more induced by NaCl in ST than NT
- `03b_DEGs_down_in_ST.tsv` — genes more suppressed by NaCl in ST than NT
- `04_PCA_T30.png` — PCA of T30 samples only
- `05_volcano_T30_ST_vs_NT.png` — volcano plot
- `06_heatmap_top500.png` — heatmap of top 500 DEGs

---

### Step 2 — Filter uncharacterized genes (`02_filter_gene_list.R`) — optional

**Input:** gene ID list (one NCBI gene ID per line)

**What it does:** Queries NCBI via `rentrez` to fetch gene metadata and keeps
only uncharacterized protein-coding genes. Skipped by default in `run_pipeline.R`
because GO enrichment (Step 3) works directly on the full DEG list.

> Note: requires ~17 min for 2,500 genes due to NCBI API rate limits (0.4s/gene).

**Output:** TSV with columns `gene_id, gene_name, gene_desc, gene_type,
refseq_genomic, refseq_mrna, refseq_peptide, related_acc`

---

### Step 3 — Arabidopsis ortholog mapping + GO enrichment (`03_annotate_and_go.R`)

**Input:** gene ID list from Step 1 (or Step 2 if filtered)

**What it does:**
- Maps quinoa NCBI gene IDs → quinoa Ensembl IDs → Arabidopsis orthologs
  via Ensembl Plants (biomaRt) — curated orthology, no BLAST required
- Runs GO enrichment across BP, MF, CC ontologies using `clusterProfiler`
  and `org.At.tair.db` (Arabidopsis GO annotations)

**Outputs in `results/03_annotation/`:**
- `04_quinoa_arabidopsis_orthologs.tsv` — quinoa → Arabidopsis gene mapping
- `04_GO_enrichment_results.tsv` — full GO enrichment table (padj, gene ratio, etc.)
- `04_GO_dotplot.png` — dotplot of top 15 GO terms per ontology

---

## Optional steps

Only needed if protein FASTA sequences are required (e.g. phylogenetics,
structural prediction):

```r
source("scripts/optional/get_fasta.R")
bash scripts/optional/blast_against_arabidopsis.sh
```

---

## Notes

- `Tol` (ST/NT) is the correct biological grouping variable — not `Tolerance`
  (accession numbers 14, 227, 287, 460, 577, 670)
- Step 3 can be run standalone on any gene ID list — Steps 1 and 2 are not
  required if you already have a gene list
- Orthology uses curated Ensembl Plants mappings, not BLAST hits
- GO enrichment uses Arabidopsis annotations since quinoa lacks a comprehensive
  GO database
