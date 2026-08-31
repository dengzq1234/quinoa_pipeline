###############################################################
### 04_cluster_analysis.R
### Answers Raquel's three follow-up questions:
###   Q1 — NaCl main effect for PPR40, AS-1, trmH (qPCR genes)
###   Q2 — Expression clustering (degPatterns)
###   Q3 — Outlier sample identification from PCA
###############################################################

library(DESeq2)
library(DEGreport)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)

# ── Configuration ──────────────────────────────────────────────────────────────
if (!exists("data_dir"))   data_dir   <- "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/data"
if (!exists("deg_file"))   deg_file   <- "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/results/03_DEG_list_for_annotation.txt"
if (!exists("output_dir")) output_dir <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/04_cluster"

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

sample_cols  <- samples$sample_name
counts       <- as.matrix(counts_raw[, sample_cols])
mode(counts) <- "integer"
rownames(counts) <- counts_raw$gene_id

keep     <- rowSums(counts) >= 10
counts_f <- counts[keep, ]

rownames(samples) <- samples$sample_name
samples$Treatment <- factor(samples$Treatment, levels = c("CONT","NaCl"))
samples$Tol       <- factor(samples$Tol,       levels = c("NT","ST"))
samples$Time      <- factor(samples$Time,       levels = c("T0","T30"))

###############################################################
### Q3 — Outlier sample identification ────────────────────────
###############################################################
message("\n[Q3] Identifying outlier samples from PCA...")

dds_qc  <- DESeqDataSetFromMatrix(counts_f, colData = samples,
                                   design = ~ Tol + Treatment)
vst_all <- vst(dds_qc, blind = TRUE)

pca_data <- plotPCA(vst_all, intgroup = c("Tol","Treatment","Time"),
                    ntop = 1000, returnData = TRUE)
pct      <- round(100 * attr(pca_data, "percentVar"))

# Flag samples > 2 SD from their group centroid (Tol × Treatment × Time)
pca_flagged <- pca_data %>%
  mutate(group = paste(Tol, Treatment, Time, sep = "_")) %>%
  group_by(group) %>%
  mutate(
    dist    = sqrt((PC1 - mean(PC1))^2 + (PC2 - mean(PC2))^2),
    grp_sd  = sqrt(sd(PC1)^2 + sd(PC2)^2),
    outlier = !is.na(grp_sd) & grp_sd > 0 & dist > 2 * grp_sd
  ) %>%
  ungroup()

outliers <- pca_flagged %>% filter(outlier) %>%
  select(name, group, PC1, PC2, dist) %>%
  arrange(desc(dist))

message("Potential outlier samples found: ", nrow(outliers))
print(outliers)
write.table(outliers, file.path(output_dir, "Q3_outlier_samples.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# PCA with outliers circled
p_pca <- ggplot(pca_flagged, aes(PC1, PC2, color = Tol, shape = Treatment)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_point(data = filter(pca_flagged, outlier),
             shape = 21, size = 5.5, color = "black", fill = NA, stroke = 1.3) +
  geom_text_repel(data = filter(pca_flagged, outlier),
                  aes(label = name), size = 3, color = "black",
                  box.padding = 0.5, max.overlaps = Inf) +
  facet_wrap(~Time, labeller = label_both) +
  scale_color_manual(values = c("NT" = "#0072B2", "ST" = "#D32F2F")) +
  xlab(paste0("PC1: ", pct[1], "% variance")) +
  ylab(paste0("PC2: ", pct[2], "% variance")) +
  theme_bw(base_size = 12) +
  ggtitle("PCA — outlier samples circled (> 2 SD from group centroid)")

ggsave(file.path(output_dir, "Q3_PCA_outliers.png"),
       p_pca, width = 12, height = 5, dpi = 300)
message("[Q3] Done.")

###############################################################
### Q1 — NaCl main effect for PPR40, AS-1, trmH (T30) ────────
###############################################################
message("\n[Q1] Extracting NaCl main effect for qPCR genes...")

qpcr_ids   <- c("110697655", "110710750", "110703716")
qpcr_names <- c("PPR40",     "AS-1",      "trmH")

s_T30      <- samples[samples$Time == "T30", ]
counts_T30 <- counts_f[, rownames(s_T30)]

# Simple model without interaction — isolates the NaCl main effect
dds_main <- DESeqDataSetFromMatrix(counts_T30, colData = s_T30,
                                    design = ~ Tol + Treatment)
dds_main <- DESeq(dds_main)

res_main  <- results(dds_main, name = "Treatment_NaCl_vs_CONT", alpha = 0.05)

qpcr_res <- as.data.frame(res_main[qpcr_ids, ]) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_name = qpcr_names,
         significant = ifelse(!is.na(padj) & padj < 0.05, "YES", "NO")) %>%
  select(gene_id, gene_name, baseMean, log2FoldChange, lfcSE, pvalue, padj, significant)

message("NaCl main effect (T30, model: ~ Tol + Treatment):")
print(qpcr_res)
write.table(qpcr_res, file.path(output_dir, "Q1_qPCR_NaCl_main_effect.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("[Q1] Done.")

###############################################################
### Q2 — Expression clustering with degPatterns ───────────────
###############################################################
message("\n[Q2] Running degPatterns clustering...")

# Re-run full model on T30 for VST matrix
dds_T30 <- DESeqDataSetFromMatrix(counts_T30, colData = s_T30,
                                   design = ~ Tol + Treatment + Tol:Treatment)
dds_T30  <- DESeq(dds_T30, test = "LRT", reduced = ~ Tol + Treatment)
vst_T30  <- vst(dds_T30, blind = FALSE)

# Load DEG list — use top 1000 by padj for clustering (manageable runtime)
deg_ids <- trimws(readLines(deg_file, warn = FALSE))
deg_ids <- deg_ids[nzchar(deg_ids)]

res_T30 <- results(dds_T30, name = "TolST.TreatmentNaCl")
top1000  <- as.data.frame(res_T30) %>%
  rownames_to_column("gene_id") %>%
  filter(gene_id %in% deg_ids, !is.na(padj)) %>%
  arrange(padj) %>%
  head(1000) %>%
  pull(gene_id)

expr_mat <- assay(vst_T30)[top1000, ]

# degPatterns: Treatment (CONT→NaCl) on x-axis, Tol (ST/NT) as color
png(file.path(output_dir, "Q2_degPatterns_clusters.png"),
    width = 3200, height = 2600, res = 300)
clusters <- degPatterns(expr_mat,
                        metadata = s_T30,
                        time     = "Treatment",
                        col      = "Tol",
                        plot     = TRUE,
                        reduce   = TRUE,
                        minc     = 15)
dev.off()

write.table(clusters$df,
            file.path(output_dir, "Q2_cluster_assignments.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

n_clusters <- length(unique(clusters$df$cluster))
message("Clusters found: ", n_clusters)
message("Genes per cluster:")
print(table(clusters$df$cluster))
message("[Q2] Done.")

###############################################################
### Final summary ─────────────────────────────────────────────
###############################################################
message("\n=== Results summary ===")
message("Q1  NaCl main effect (qPCR genes): Q1_qPCR_NaCl_main_effect.tsv")
message("Q2  Cluster plot:                  Q2_degPatterns_clusters.png")
message("Q2  Cluster assignments:           Q2_cluster_assignments.tsv")
message("Q3  Outlier samples:               Q3_outlier_samples.tsv")
message("Q3  Outlier PCA:                   Q3_PCA_outliers.png")
message("All results in: ", output_dir)
