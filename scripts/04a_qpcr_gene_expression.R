###############################################################
### 04a_qpcr_gene_expression.R
### Test whether Raquel's three qPCR validation genes show
### significant NaCl response at T30.
###
### Genes:
###   PPR40  Phytozome: AUR62042953  NCBI: 110697655
###   AS-1   Phytozome: AUR62030471  NCBI: 110710750
###   trmH   Phytozome: AUR62044515  NCBI: 110703716
###
### Model: ~ Tol + Treatment  (NaCl main effect, T30 only)
### Output: results/04_cluster/Q1_qPCR_NaCl_main_effect.tsv
###############################################################

library(DESeq2)
library(dplyr)
library(tibble)

# ── Configuration ──────────────────────────────────────────────────────────────
if (!exists("data_dir"))   data_dir   <- "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/data"
if (!exists("output_dir")) output_dir <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/04_cluster"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ──────────────────────────────────────────────────────────────────
message("Loading count matrix and sample metadata...")
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

# Basic expression filter
keep     <- rowSums(counts) >= 10
counts_f <- counts[keep, ]

# Subset to T30
s_T30      <- samples[samples$Time == "T30", ]
counts_T30 <- counts_f[, rownames(s_T30)]

# ── DESeq2: NaCl main effect (no interaction term) ────────────────────────────
# This isolates the overall NaCl effect averaged across ST and NT.
# The interaction model (used in Step 1) tests ST-vs-NT differences;
# here we ask simply: does NaCl change these genes at all?
message("Running DESeq2 (~ Tol + Treatment) on T30 samples...")
dds_main <- DESeqDataSetFromMatrix(counts_T30, colData = s_T30,
                                    design = ~ Tol + Treatment)
dds_main <- DESeq(dds_main)

res_main <- results(dds_main, name = "Treatment_NaCl_vs_CONT", alpha = 0.05)

# ── Extract qPCR gene results ──────────────────────────────────────────────────
qpcr_ids   <- c("110697655", "110710750", "110703716")
qpcr_names <- c("PPR40",     "AS-1",      "trmH")

qpcr_res <- as.data.frame(res_main[qpcr_ids, ]) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    gene_name   = qpcr_names,
    phytozome   = c("AUR62042953", "AUR62030471", "AUR62044515"),
    significant = ifelse(!is.na(padj) & padj < 0.05, "YES", "NO")
  ) %>%
  select(gene_id, gene_name, phytozome,
         baseMean, log2FoldChange, lfcSE, pvalue, padj, significant)

message("\nNaCl main effect (T30, model: ~ Tol + Treatment):")
print(qpcr_res)

out_file <- file.path(output_dir, "Q1_qPCR_NaCl_main_effect.tsv")
write.table(qpcr_res, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
message("\nSaved: ", out_file)
