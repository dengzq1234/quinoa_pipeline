# ── Quinoa Annotation Pipeline ─────────────────────────────────────────────────
# Edit the paths below, then source this file to run the full pipeline.

# Input: one NCBI gene ID per line
input_file <- "data/input/gene_list_1411.txt"

# Output directory prefix
output_dir <- "results/gene_list_1411/"

# ── Step 1: Filter uncharacterized protein-coding genes ────────────────────────
source("scripts/01_filter_gene_list.R")

# ── Step 2: Arabidopsis ortholog mapping + GO enrichment ───────────────────────
source("scripts/02_annotate_and_go.R")

message("Pipeline complete. Results in: ", output_dir)
