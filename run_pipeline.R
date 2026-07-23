# ── Quinoa Salt Tolerance Pipeline ─────────────────────────────────────────────
# Full pipeline: DESeq2 differential expression → Arabidopsis annotation & GO enrichment
#
# Edit the two paths below, then source this file.

# ── Configuration ──────────────────────────────────────────────────────────────

# Directory containing gene_count.csv and samples.csv
rnaseq_data_dir <- "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/data"

# Root directory for all results
results_base <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results"

# ── Step 1: DESeq2 — ST vs NT differential expression (T30) ───────────────────
# Input : gene_count.csv + samples.csv
# Output: DEG list, PCA, volcano plot, heatmap
data_dir   <- rnaseq_data_dir
output_dir <- file.path(results_base, "01_deseq2")
source("scripts/01_deseq2_analysis.R")

# ── Step 2: Filter to uncharacterized genes (optional) ─────────────────────────
# Skips by default — enrichGO works on the full DEG list.
# Uncomment if you need to isolate uncharacterized protein-coding genes only.
# NOTE: requires ~17 min (NCBI API call per gene).
#
# input_file  <- file.path(results_base, "01_deseq2", "03_DEG_list_for_annotation.txt")
# output_file <- file.path(results_base, "02_filter", "filtered_genes.tsv")
# source("scripts/02_filter_gene_list.R")

# ── Step 3: Arabidopsis ortholog mapping + GO enrichment ───────────────────────
# Input : DEG gene ID list from Step 1
# Output: ortholog table, GO enrichment TSV, dotplot
input_file <- file.path(results_base, "01_deseq2", "03_DEG_list_for_annotation.txt")
output_dir <- file.path(results_base, "03_annotation")
source("scripts/03_annotate_and_go.R")

message("Pipeline complete. Results in: ", results_base)
