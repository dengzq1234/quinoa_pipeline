###############################################################
### 04b_cluster_timecourse.R
### Expression clustering across ALL conditions and time points
###
### Extends 04_cluster_analysis.R (T30-only) by including T0,
### so each cluster's mean expression is shown across four points:
###   CONT-T0 → CONT-T30
###   NaCl-T0 → NaCl-T30
### split by ST / NT — this is the full trajectory Raquel asked for.
###
### DEGs: same top-1000 interaction DEGs from the T30 model
### VST:  computed on ALL samples (T0 + T30)
###
### Output: results/04_cluster/Q2_timecourse_clusters.png
###         results/04_cluster/Q2_timecourse_assignments.tsv
###############################################################

library(DESeq2)
library(DEGreport)
library(dplyr)
library(tibble)
library(ggplot2)

# ── Configuration ──────────────────────────────────────────────────────────────
if (!exists("data_dir"))   data_dir   <- "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/data"
if (!exists("deg_file"))   deg_file   <- "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/results/03_DEG_list_for_annotation.txt"
if (!exists("output_dir")) output_dir <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/04_cluster"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ──────────────────────────────────────────────────────────────────
message("Loading data...")
counts_raw <- read.csv(file.path(data_dir, "gene_count.csv"),
                       sep = ";", header = TRUE, stringsAsFactors = FALSE)
samples    <- read.csv(file.path(data_dir, "samples.csv"),
                       sep = ";", header = TRUE, stringsAsFactors = FALSE)

sample_cols  <- samples$sample_name
counts       <- as.matrix(counts_raw[, sample_cols])
mode(counts) <- "integer"
rownames(counts) <- counts_raw$gene_id

rownames(samples) <- samples$sample_name
samples$Treatment <- factor(samples$Treatment, levels = c("CONT", "NaCl"))
samples$Tol       <- factor(samples$Tol,       levels = c("NT", "ST"))
samples$Time      <- factor(samples$Time,       levels = c("T0", "T30"))

# Low-expression filter
keep     <- rowSums(counts) >= 10
counts_f <- counts[keep, ]
message("Expressed genes: ", nrow(counts_f))

# ── Step 1: Identify top 1,000 DEGs from the T30 interaction model ─────────────
# (same genes used in 04_cluster_analysis.R — keeps results comparable)
message("\n[Step 1] Re-running T30 interaction model to rank DEGs...")

s_T30      <- samples[samples$Time == "T30", ]
counts_T30 <- counts_f[, rownames(s_T30)]

dds_T30 <- DESeqDataSetFromMatrix(counts_T30, colData = s_T30,
                                   design = ~ Tol + Treatment + Tol:Treatment)
dds_T30 <- DESeq(dds_T30, test = "LRT", reduced = ~ Tol + Treatment)

deg_ids <- trimws(readLines(deg_file, warn = FALSE))
deg_ids <- deg_ids[nzchar(deg_ids)]

res_T30 <- results(dds_T30, name = "TolST.TreatmentNaCl")
top1000 <- as.data.frame(res_T30) %>%
  rownames_to_column("gene_id") %>%
  filter(gene_id %in% deg_ids, !is.na(padj)) %>%
  arrange(padj) %>%
  head(1000) %>%
  pull(gene_id)

message("Top 1,000 DEGs selected.")

# ── Step 2: VST on ALL samples (T0 + T30) ─────────────────────────────────────
# Use a model that accounts for Time, Treatment, and Tol.
# blind = FALSE uses the model to stabilise variance more accurately.
message("\n[Step 2] Computing VST on all samples (T0 + T30)...")

dds_all <- DESeqDataSetFromMatrix(counts_f, colData = samples,
                                   design = ~ Tol + Time + Treatment)
dds_all <- estimateSizeFactors(dds_all)
dds_all <- estimateDispersions(dds_all)
vst_all <- vst(dds_all, blind = FALSE)

# Extract expression matrix for top 1,000 DEGs across all samples
expr_mat <- assay(vst_all)[top1000, ]
message("Expression matrix: ", nrow(expr_mat), " genes × ", ncol(expr_mat), " samples")

# ── Step 3: Build combined condition factor ────────────────────────────────────
# degPatterns uses one column as the x-axis ("time").
# We create a four-level factor so the x-axis reads:
#   CONT-T0 → CONT-T30   and   NaCl-T0 → NaCl-T30
samples$Condition <- factor(
  paste(samples$Treatment, samples$Time, sep = "-"),
  levels = c("CONT-T0", "CONT-T30", "NaCl-T0", "NaCl-T30")
)

message("\nSample counts per condition × Tol:")
print(table(samples$Condition, samples$Tol))

# ── Step 4: degPatterns ────────────────────────────────────────────────────────
message("\n[Step 3] Running degPatterns (all time points)...")
message("  time  = Condition  (CONT-T0 / CONT-T30 / NaCl-T0 / NaCl-T30)")
message("  col   = Tol        (ST / NT)")

png_out <- file.path(output_dir, "Q2_timecourse_clusters.png")
png(png_out, width = 3600, height = 3000, res = 300)

clusters_tc <- degPatterns(
  expr_mat,
  metadata = samples,
  time     = "Condition",
  col      = "Tol",
  plot     = TRUE,
  reduce   = TRUE,
  minc     = 15
)

dev.off()
message("Saved: ", basename(png_out))

# ── Step 5: Save cluster assignments ──────────────────────────────────────────
tsv_out <- file.path(output_dir, "Q2_timecourse_assignments.tsv")
write.table(clusters_tc$df, tsv_out,
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved: ", basename(tsv_out))

n_cl <- length(unique(clusters_tc$df$cluster))
message("\nClusters found: ", n_cl)
message("Genes per cluster:")
print(table(clusters_tc$df$cluster))

message("\n=== Done ===")
message("Timecourse cluster plot: ", png_out)
message("Cluster assignments:     ", tsv_out)
