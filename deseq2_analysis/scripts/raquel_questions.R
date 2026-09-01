###############################################################
### Analysis to answer Raquel's three questions
### Q1: Main effect of NaCl for PPR40, AS-1, trmH
### Q2: Clustering analysis (degPatterns)
### Q3: Identify outlier samples from PCA
###############################################################

library(DESeq2)
library(DEGreport)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)

data_dir   <- "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/data"
output_dir <- "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/results/raquel_questions"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ──────────────────────────────────────────────────────────────────
counts_raw <- read.csv(file.path(data_dir, "gene_count.csv"),
                       sep = ";", header = TRUE, stringsAsFactors = FALSE)
samples    <- read.csv(file.path(data_dir, "samples.csv"),
                       sep = ";", header = TRUE, stringsAsFactors = FALSE)

ann_col_names <- c("gene_name","gene_chr","gene_start","gene_end",
                   "gene_strand","gene_length","gene_biotype","gene_description","Family")
gene_ann <- counts_raw[, c("gene_id", intersect(ann_col_names, colnames(counts_raw)))]
rownames(gene_ann) <- counts_raw$gene_id

sample_cols <- samples$sample_name
counts      <- as.matrix(counts_raw[, sample_cols])
mode(counts) <- "integer"
rownames(counts) <- counts_raw$gene_id

keep     <- rowSums(counts) >= 10
counts_f <- counts[keep, ]

rownames(samples)  <- samples$sample_name
samples$Treatment  <- factor(samples$Treatment, levels = c("CONT","NaCl"))
samples$Tol        <- factor(samples$Tol,       levels = c("NT","ST"))
samples$Time       <- factor(samples$Time,       levels = c("T0","T30"))
samples$Tolerance  <- factor(samples$Tolerance)

###############################################################
### Q3 FIRST: Outlier identification from PCA (all samples) ---
###############################################################
cat("\n=== Q3: Outlier identification ===\n")

dds_qc  <- DESeqDataSetFromMatrix(counts_f, colData = samples,
                                   design = ~ Tol + Treatment)
vst_all <- vst(dds_qc, blind = TRUE)

pca_data <- plotPCA(vst_all, intgroup = c("Tol","Treatment","Time"),
                    ntop = 1000, returnData = TRUE)
pct      <- round(100 * attr(pca_data, "percentVar"))

# Flag outliers: samples > 2 SD from their group centroid (Tol × Treatment × Time)
pca_data$group <- paste(pca_data$Tol, pca_data$Treatment, pca_data$Time, sep = "_")

outlier_flags <- pca_data %>%
  group_by(group) %>%
  mutate(
    mean_PC1 = mean(PC1), mean_PC2 = mean(PC2),
    sd_PC1   = sd(PC1),   sd_PC2   = sd(PC2),
    dist     = sqrt((PC1 - mean_PC1)^2 + (PC2 - mean_PC2)^2),
    group_sd = sqrt(sd_PC1^2 + sd_PC2^2),
    outlier  = dist > 2 * group_sd
  ) %>%
  ungroup()

outliers <- outlier_flags %>% filter(outlier) %>% select(name, group, dist)
cat("Potential outlier samples:\n")
print(outliers)
write.table(outliers, file.path(output_dir, "Q3_outlier_samples.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# PCA plot with outliers flagged
p_outlier <- ggplot(outlier_flags, aes(PC1, PC2, color = Tol, shape = Treatment)) +
  geom_point(aes(size = ifelse(outlier, 4, 2.5)), alpha = 0.9) +
  scale_size_identity() +
  geom_text_repel(aes(label = ifelse(outlier, name, "")),
                  size = 3, color = "black", max.overlaps = Inf,
                  box.padding = 0.4) +
  geom_point(data = filter(outlier_flags, outlier),
             shape = 21, size = 5, color = "black", fill = NA, stroke = 1.2) +
  facet_wrap(~Time, labeller = label_both) +
  scale_color_manual(values = c("NT" = "#0072B2", "ST" = "#D32F2F")) +
  xlab(paste0("PC1: ", pct[1], "% variance")) +
  ylab(paste0("PC2: ", pct[2], "% variance")) +
  theme_bw(base_size = 12) +
  ggtitle("PCA — potential outliers circled (>2 SD from group centroid)")

ggsave(file.path(output_dir, "Q3_PCA_outliers.png"),
       p_outlier, width = 12, height = 5, dpi = 300)
cat("Outlier PCA saved.\n")

###############################################################
### Q1: NaCl main effect for the three qPCR genes (T30) ------
###############################################################
cat("\n=== Q1: qPCR gene NaCl main effect ===\n")

qpcr_ids <- c("110697655", "110710750", "110703716")
qpcr_names <- c("PPR40", "AS-1", "trmH")

# Simple model on T30: ~ Tol + Treatment (no interaction) — extracts Treatment main effect
s_T30      <- samples[samples$Time == "T30", ]
counts_T30 <- counts_f[, rownames(s_T30)]

dds_main <- DESeqDataSetFromMatrix(counts_T30, colData = s_T30,
                                    design = ~ Tol + Treatment)
dds_main <- DESeq(dds_main)

res_main <- results(dds_main, name = "Treatment_NaCl_vs_CONT", alpha = 0.05)

qpcr_res <- as.data.frame(res_main[qpcr_ids, ]) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_name = qpcr_names) %>%
  select(gene_id, gene_name, baseMean, log2FoldChange, lfcSE, pvalue, padj)

cat("NaCl main effect (T30, controlling for Tol):\n")
print(qpcr_res)
write.table(qpcr_res, file.path(output_dir, "Q1_qPCR_genes_NaCl_main_effect.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

###############################################################
### Q2: Clustering analysis (degPatterns) — T30 DEGs ---------
###############################################################
cat("\n=== Q2: Clustering analysis ===\n")

# Re-run full interaction model on T30 to get VST matrix
dds_T30 <- DESeqDataSetFromMatrix(counts_T30, colData = s_T30,
                                   design = ~ Tol + Treatment + Tol:Treatment)
dds_T30  <- DESeq(dds_T30, test = "LRT", reduced = ~ Tol + Treatment)
vst_T30  <- vst(dds_T30, blind = FALSE)

# Load DEG list
deg_ids  <- trimws(readLines(
  "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/results/03_DEG_list_for_annotation.txt",
  warn = FALSE))
deg_ids  <- deg_ids[deg_ids %in% rownames(assay(vst_T30))]

# Use top 1000 most significant for clustering (degPatterns is slow on full list)
res_T30 <- results(dds_T30, name = "TolST.TreatmentNaCl")
top_genes <- as.data.frame(res_T30) %>%
  rownames_to_column("gene_id") %>%
  filter(gene_id %in% deg_ids, !is.na(padj)) %>%
  arrange(padj) %>%
  head(1000) %>%
  pull(gene_id)

expr_mat <- assay(vst_T30)[top_genes, ]

# degPatterns: Treatment as x-axis, Tol as color grouping
png(file.path(output_dir, "Q2_degPatterns_clusters.png"),
    width = 3000, height = 2400, res = 300)
clusters <- degPatterns(expr_mat,
                        metadata  = s_T30,
                        time      = "Treatment",
                        col       = "Tol",
                        plot      = TRUE,
                        reduce    = TRUE,
                        minc      = 15)
dev.off()

# Save cluster assignments
write.table(clusters$df,
            file.path(output_dir, "Q2_cluster_assignments.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

n_clusters <- length(unique(clusters$df$cluster))
cat("Clusters found:", n_clusters, "\n")
cat("Genes per cluster:\n")
print(table(clusters$df$cluster))

###############################################################
### Summary
###############################################################
cat("\n=== Summary ===\n")
cat("Q1 — qPCR gene NaCl main effect: Q1_qPCR_genes_NaCl_main_effect.tsv\n")
cat("Q2 — Cluster plot:               Q2_degPatterns_clusters.png\n")
cat("Q2 — Cluster assignments:        Q2_cluster_assignments.tsv\n")
cat("Q3 — Outlier samples:            Q3_outlier_samples.tsv\n")
cat("Q3 — Outlier PCA:                Q3_PCA_outliers.png\n")
cat("All results in:", output_dir, "\n")
