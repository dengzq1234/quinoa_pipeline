###############################################################
### 05_cluster_go_enrichment.R
### GO enrichment analysis for each expression cluster
### Clusters from:    results/04_cluster/Q2_cluster_assignments.tsv
### Orthologs from:    results/degs_salt_tolerance/04_quinoa_arabidopsis_orthologs_2551.tsv
### GO universe from:  results/05_cluster_go/background_orthologs_all.tsv  (Step 5a)
###
### The universe is the Arabidopsis-ortholog set of ALL expressed quinoa genes
### (rowSums >= 10), not just the DEGs — this is the statistically correct
### background for enrichGO. Outputs are suffixed "_corrected".
###############################################################

library(clusterProfiler)
library(org.At.tair.db)
library(dplyr)
library(ggplot2)
library(stringr)

# ── Configuration ──────────────────────────────────────────────────────────────
if (!exists("output_dir")) output_dir <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/05_cluster_go"
if (!exists("cluster_file")) cluster_file <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/04_cluster/Q2_cluster_assignments.tsv"
if (!exists("ortholog_file")) ortholog_file <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/degs_salt_tolerance/04_quinoa_arabidopsis_orthologs_2551.tsv"
if (!exists("background_file")) background_file <- file.path(output_dir, "background_orthologs_all.tsv")

if (!file.exists(background_file))
  stop("Background ortholog table not found: ", background_file,
       "\nRun scripts/05a_build_background_orthologs.R first.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cluster descriptions based on degPatterns patterns
cluster_labels <- c(
  "1" = "Cluster 1 (n=257): Crossover — ST↑ / NT↓ under NaCl",
  "2" = "Cluster 2 (n=495): ST-specific strong induction under NaCl",
  "3" = "Cluster 3 (n=178): NT collapses under NaCl (ST stable)",
  "4" = "Cluster 4 (n=17):  ST shuts down under NaCl (NT stable)"
)

# ── Load data ──────────────────────────────────────────────────────────────────
message("Loading cluster assignments...")
clusters  <- read.table(cluster_file,  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
message("Loading ortholog table...")
orthologs <- read.table(ortholog_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Strip R-added "X" prefix and ensure string type for joining
clusters$gene_id <- sub("^X", "", clusters$genes)
orthologs$entrezgene_id <- as.character(orthologs$entrezgene_id)

# Filter for valid Arabidopsis genes
orthologs_valid <- orthologs %>%
  filter(!is.na(athaliana_eg_homolog_ensembl_gene),
         athaliana_eg_homolog_ensembl_gene != "")

# Join: quinoa cluster genes → Arabidopsis orthologs
merged <- clusters %>%
  left_join(orthologs_valid, by = c("gene_id" = "entrezgene_id"))

message("\nGenes with Arabidopsis orthologs per cluster:")
print(table(merged$cluster[!is.na(merged$athaliana_eg_homolog_ensembl_gene)]))

# Background (GO universe): Arabidopsis orthologs of ALL expressed quinoa genes,
# built by Step 5a. This is the correct enrichGO background — using only the DEG
# orthologs here would massively inflate enrichment.
message("Loading background ortholog table (Step 5a)...")
bg_tbl     <- read.table(background_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
background <- bg_tbl %>%
  filter(!is.na(athaliana_eg_homolog_ensembl_gene),
         athaliana_eg_homolog_ensembl_gene != "") %>%
  pull(athaliana_eg_homolog_ensembl_gene) %>%
  unique()
message("GO universe size (unique Arabidopsis genes): ", length(background))

# Sanity check: cluster genes should be a subset of the universe
cluster_at <- unique(na.omit(merged$athaliana_eg_homolog_ensembl_gene))
n_missing  <- sum(!cluster_at %in% background)
if (n_missing > 0)
  message("Note: ", n_missing, " / ", length(cluster_at),
          " cluster orthologs are not in the universe — adding them.")
background <- unique(c(background, cluster_at))

# ── GO enrichment per cluster ──────────────────────────────────────────────────
all_go_results <- list()

for (cl in 1:4) {
  message("\n", paste(rep("─", 60), collapse=""))
  message("Processing ", cluster_labels[as.character(cl)])

  genes_cl <- merged %>%
    filter(cluster == cl, !is.na(athaliana_eg_homolog_ensembl_gene)) %>%
    pull(athaliana_eg_homolog_ensembl_gene) %>%
    unique()

  message("Arabidopsis genes in cluster: ", length(genes_cl))

  if (length(genes_cl) < 5) {
    message("Too few genes (< 5) — skipping GO enrichment.")
    next
  }

  go <- tryCatch(
    enrichGO(
      gene          = genes_cl,
      universe      = background,
      OrgDb         = org.At.tair.db,
      keyType       = "TAIR",
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    ),
    error = function(e) {
      message("ERROR in enrichGO: ", e$message)
      NULL
    }
  )

  if (is.null(go)) next

  go_df <- as.data.frame(go)
  n_sig <- nrow(go_df)
  message("Significant GO terms: ", n_sig)

  if (n_sig == 0) {
    message("No significant GO terms found for cluster ", cl)
    next
  }

  # Save TSV
  tsv_out <- file.path(output_dir, paste0("cluster", cl, "_GO_corrected.tsv"))
  write.table(go_df, tsv_out, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Saved: ", basename(tsv_out))

  # ── Dotplot ────────────────────────────────────────────────────────────────
  # Wrap long descriptions to avoid label overlap
  go@result$Description <- str_wrap(go@result$Description, width = 38)

  n_show <- min(12, n_sig)

  p <- dotplot(go, showCategory = n_show, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free", space = "free_y") +
    theme(
      axis.text.y   = element_text(size = 7, lineheight = 0.85),
      axis.text.x   = element_text(size = 7),
      strip.text    = element_text(size = 8, face = "bold"),
      legend.text   = element_text(size = 7),
      legend.title  = element_text(size = 8),
      panel.spacing = unit(0.5, "lines"),
      plot.title    = element_text(size = 8, face = "bold", hjust = 0),
      plot.margin   = margin(10, 15, 10, 10)
    ) +
    ggtitle(cluster_labels[as.character(cl)])

  png_out <- file.path(output_dir, paste0("cluster", cl, "_GO_corrected_dotplot.png"))
  png(png_out, width = 2800, height = 3800, res = 300)
  print(p)
  dev.off()
  message("Saved: ", basename(png_out))

  # Store for combined table
  all_go_results[[cl]] <- go_df %>% mutate(cluster = cl)
}

# ── Combined summary table ─────────────────────────────────────────────────────
if (length(all_go_results) > 0) {
  combined <- bind_rows(all_go_results)
  combined_out <- file.path(output_dir, "all_clusters_GO_corrected_combined.tsv")
  write.table(combined, combined_out, sep = "\t", row.names = FALSE, quote = FALSE)

  message("\n", paste(rep("═", 60), collapse=""))
  message("Combined GO table: ", nrow(combined), " terms across ",
          length(all_go_results), " clusters")
  message("Saved: ", basename(combined_out))

  message("\nTop 3 GO terms per cluster (sorted by p.adjust):")
  combined %>%
    group_by(cluster, ONTOLOGY) %>%
    slice_min(p.adjust, n = 2, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(cluster, ONTOLOGY, p.adjust) %>%
    select(cluster, ONTOLOGY, Description, p.adjust, Count) %>%
    print(n = 50)
} else {
  message("\nNo significant GO enrichment found in any cluster.")
}

message("\nAll results saved to: ", output_dir)
