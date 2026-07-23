###############################################################
### DESeq2 Analysis — Quinoa Salt Tolerance (T30)
### Goal: identify genes whose NaCl response differs between
###       salt-tolerant (ST) and non-tolerant (NT) accessions
###############################################################

# ── Configuration ─────────────────────────────────────────────────────────────
if (!exists("data_dir"))   data_dir   <- "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/data"
if (!exists("output_dir")) output_dir <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/01_deseq2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Packages ───────────────────────────────────────────────────────────────────
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(ggrepel)

###############################################################
### PASO 1 — Lectura y preparación de datos ----
###############################################################

counts_raw <- read.csv(file.path(data_dir, "gene_count.csv"),
                       sep = ";", header = TRUE, stringsAsFactors = FALSE)
samples    <- read.csv(file.path(data_dir, "samples.csv"),
                       sep = ";", header = TRUE, stringsAsFactors = FALSE)

# Annotation columns embedded in gene_count.csv
ann_col_names <- c("gene_name", "gene_chr", "gene_start", "gene_end",
                   "gene_strand", "gene_length", "gene_biotype",
                   "gene_description", "Family")

# Gene annotation table
gene_ann <- counts_raw[, c("gene_id", intersect(ann_col_names, colnames(counts_raw)))]
rownames(gene_ann) <- counts_raw$gene_id

# Count matrix — only sample columns, no annotation
sample_cols <- samples$sample_name
missing <- setdiff(sample_cols, colnames(counts_raw))
if (length(missing) > 0) stop("Samples missing from count matrix: ", paste(missing, collapse = ", "))

counts <- as.matrix(counts_raw[, sample_cols])
mode(counts) <- "integer"
rownames(counts) <- counts_raw$gene_id

# Low-count filter: keep genes with >= 10 reads across all samples
keep    <- rowSums(counts) >= 10
counts_f <- counts[keep, ]
cat("Genes after filtering:", sum(keep), "of", nrow(counts), "\n")

# Factor metadata — Tol is the biological ST/NT variable (NOT Tolerance which is accession ID)
rownames(samples) <- samples$sample_name
samples$Treatment <- factor(samples$Treatment, levels = c("CONT", "NaCl"))
samples$Tol       <- factor(samples$Tol,       levels = c("NT", "ST"))   # NT is reference
samples$Time      <- factor(samples$Time,       levels = c("T0", "T30"))
samples$Tolerance <- factor(samples$Tolerance)

###############################################################
### PASO 2 — PCA on full dataset (QC) ----
###############################################################

dds_qc <- DESeqDataSetFromMatrix(
  countData = counts_f,
  colData   = samples,
  design    = ~ Tol + Treatment    # simple design for blind QC
)
vst_qc <- vst(dds_qc, blind = TRUE)

pca_data   <- plotPCA(vst_qc, intgroup = c("Tol", "Treatment", "Time"),
                      ntop = 1000, returnData = TRUE)
pct_var    <- round(100 * attr(pca_data, "percentVar"))

p_pca_all <- ggplot(pca_data, aes(PC1, PC2, color = Tol, shape = Treatment)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = name), size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  facet_wrap(~Time, labeller = label_both) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  scale_color_manual(values = c("NT" = "#0072B2", "ST" = "#D32F2F")) +
  theme_bw(base_size = 13) +
  ggtitle("PCA — All samples (blind VST)")

ggsave(file.path(output_dir, "01_PCA_all_samples.png"),
       p_pca_all, width = 12, height = 5, dpi = 300)
cat("PCA (all samples) saved.\n")

###############################################################
### PASO 3 — DESeq2 on T30 samples: Tol × Treatment interaction ----
###############################################################

# Filter metadata and counts to T30 only (post-stress phase)
s_T30      <- samples[samples$Time == "T30", ]
counts_T30 <- counts_f[, rownames(s_T30)]
cat("T30 samples:", nrow(s_T30), "\n")

dds_T30 <- DESeqDataSetFromMatrix(
  countData = counts_T30,
  colData   = s_T30,
  design    = ~ Tol + Treatment + Tol:Treatment
)

# LRT: test whether the ST/NT × NaCl interaction term is significant
# Genes with low padj responded DIFFERENTLY to NaCl in ST vs NT accessions
dds_T30 <- DESeq(dds_T30, test = "LRT", reduced = ~ Tol + Treatment)

cat("Coefficients in model:\n")
print(resultsNames(dds_T30))

###############################################################
### PASO 4 — Extract results for TolST:TreatmentNaCl ----
###############################################################

# Positive log2FC  → gene more induced (or less repressed) by NaCl in ST than NT
# Negative log2FC → gene more repressed (or less induced) by NaCl in ST than NT
res_T30 <- results(dds_T30, name = "TolST.TreatmentNaCl", alpha = 0.05)

res_df <- as.data.frame(res_T30) %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj)) %>%
  left_join(gene_ann, by = "gene_id") %>%
  arrange(padj)

n_sig <- sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, na.rm = TRUE)
cat("Significant DEGs (padj < 0.05, |log2FC| >= 1):", n_sig, "\n")

# Full results table
write.table(res_df,
            file.path(output_dir, "02_DESeq2_T30_TolxTreatment_full.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Full results saved.\n")

###############################################################
### PASO 5 — Export DEG list for annotation pipeline ----
###############################################################

degs <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1) %>%
  pull(gene_id)

writeLines(as.character(degs),
           file.path(output_dir, "03_DEG_list_for_annotation.txt"))
cat("DEG gene IDs saved:", length(degs), "genes →",
    file.path(output_dir, "03_DEG_list_for_annotation.txt"), "\n")

# Also save separate up/down lists
degs_up <- res_df %>%
  filter(padj < 0.05, log2FoldChange >= 1) %>%
  arrange(padj)

degs_down <- res_df %>%
  filter(padj < 0.05, log2FoldChange <= -1) %>%
  arrange(padj)

write.table(degs_up,
            file.path(output_dir, "03a_DEGs_up_in_ST.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(degs_down,
            file.path(output_dir, "03b_DEGs_down_in_ST.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Up in ST:", nrow(degs_up), "| Down in ST:", nrow(degs_down), "\n")

###############################################################
### PASO 6 — PCA on T30 only ----
###############################################################

vst_T30 <- vst(dds_T30, blind = FALSE)

pca_T30   <- plotPCA(vst_T30, intgroup = c("Tol", "Treatment"), ntop = 1000, returnData = TRUE)
pct_T30   <- round(100 * attr(pca_T30, "percentVar"))

p_pca_T30 <- ggplot(pca_T30, aes(PC1, PC2, color = Tol, shape = Treatment)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = name), size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  xlab(paste0("PC1: ", pct_T30[1], "% variance")) +
  ylab(paste0("PC2: ", pct_T30[2], "% variance")) +
  scale_color_manual(values = c("NT" = "#0072B2", "ST" = "#D32F2F")) +
  theme_bw(base_size = 13) +
  ggtitle("PCA — T30 samples")

ggsave(file.path(output_dir, "04_PCA_T30.png"),
       p_pca_T30, width = 8, height = 6, dpi = 300)
cat("PCA (T30) saved.\n")

###############################################################
### PASO 7 — Volcano plot ----
###############################################################

res_plot <- res_df %>%
  mutate(status = case_when(
    padj < 0.05 & log2FoldChange >= 1  ~ "Up in ST under NaCl",
    padj < 0.05 & log2FoldChange <= -1 ~ "Down in ST under NaCl",
    TRUE ~ "Not significant"
  ))

top_labels <- res_plot %>%
  filter(status != "Not significant") %>%
  arrange(padj) %>%
  head(20)

p_volcano <- ggplot(res_plot, aes(log2FoldChange, -log10(padj), color = status)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Up in ST under NaCl"   = "#D32F2F",
      "Down in ST under NaCl" = "#0072B2",
      "Not significant"       = "gray75"
    ),
    name = NULL
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data      = top_labels,
    aes(label = gene_name),
    size      = 3,
    color     = "black",
    box.padding   = 0.4,
    point.padding = 0.3,
    max.overlaps  = Inf
  ) +
  theme_bw(base_size = 13) +
  labs(
    title    = "DESeq2: ST vs NT × NaCl interaction (T30)",
    subtitle = "Positive FC = higher NaCl induction in ST than NT",
    x        = expression(Log[2]~Fold~Change~(TolST:TreatmentNaCl)),
    y        = expression(-Log[10]~padj)
  )

ggsave(file.path(output_dir, "05_volcano_T30_ST_vs_NT.png"),
       p_volcano, width = 10, height = 7, dpi = 300)
cat("Volcano plot saved.\n")

###############################################################
### PASO 8 — Heatmap of top 500 significant genes ----
###############################################################

top500 <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(500) %>%
  pull(gene_id)

expr_500    <- assay(vst_T30)[top500, ]
expr_scaled <- t(scale(t(expr_500)))

ann_col_hm <- s_T30[, c("Tol", "Treatment"), drop = FALSE]

hm_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

png(file.path(output_dir, "06_heatmap_top500.png"),
    width = 3800, height = 2400, res = 300)
pheatmap(
  expr_scaled,
  annotation_col  = ann_col_hm,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = FALSE,
  show_colnames   = TRUE,
  fontsize        = 10,
  fontsize_col    = 9,
  angle_col       = 90,
  main            = "Top 500 DEGs by LRT padj (T30, Tol × Treatment interaction)",
  color           = hm_colors,
  border_color    = NA,
  treeheight_row  = 60,
  treeheight_col  = 40
)
dev.off()
cat("Heatmap saved.\n")

###############################################################
### RESUMEN FINAL ----
###############################################################

cat("\n=== Summary ===\n")
cat("Total genes tested:       ", nrow(res_df), "\n")
cat("DEGs (padj<0.05, |FC|≥1):", n_sig, "\n")
cat("  Up in ST under NaCl:   ", nrow(degs_up), "\n")
cat("  Down in ST under NaCl: ", nrow(degs_down), "\n")
cat("DEG list for annotation:  ", file.path(output_dir, "03_DEG_list_for_annotation.txt"), "\n")
cat("Results directory:        ", output_dir, "\n")
cat("\nNext step: feed 03_DEG_list_for_annotation.txt into\n")
cat("  quinoa_pipeline/scripts/02_annotate_and_go.R\n")

