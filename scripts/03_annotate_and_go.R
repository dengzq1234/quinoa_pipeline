library(biomaRt)
library(clusterProfiler)
library(org.At.tair.db)
library(dplyr)

# ── Configuration ──────────────────────────────────────────────────────────────
if (!exists("input_file")) input_file <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/data/input/gene_list_1411.txt"
if (!exists("output_dir")) output_dir <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/gene_list_1411/"

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# ── Step 1: Load gene list ─────────────────────────────────────────────────────
ncbi_ids <- trimws(readLines(input_file, warn=FALSE))
ncbi_ids <- ncbi_ids[nzchar(ncbi_ids)]
message("Loaded ", length(ncbi_ids), " genes")

# ── Step 2: Connect to Ensembl Plants and check quinoa dataset ─────────────────
message("[+] Connecting to Ensembl Plants...")
mart <- useMart(biomart="plants_mart", host="https://plants.ensembl.org")

# List available datasets to find quinoa
datasets <- listDatasets(mart)
quinoa_ds <- datasets[grep("quinoa|chenopodium", datasets$description, ignore.case=TRUE), ]
message("Quinoa datasets found:")
print(quinoa_ds)

# ── Step 3: Connect to quinoa dataset and map to Arabidopsis orthologs ─────────
quinoa_mart <- useMart(
  biomart = "plants_mart",
  dataset = quinoa_ds$dataset[1],   # use first match
  host    = "https://plants.ensembl.org"
)

# Check available filters for NCBI gene IDs
filters <- listFilters(quinoa_mart)
message("\n[+] Filters matching 'entrez':")
print(filters[grep("entrez", filters$name, ignore.case=TRUE), ])

# Check available attributes for Arabidopsis orthology
attrs <- listAttributes(quinoa_mart)
message("\n[+] Attributes matching 'arabidopsis':")
print(attrs[grep("arabidopsis", attrs$name, ignore.case=TRUE), ])

# ── Step 4: Get Arabidopsis orthologs ─────────────────────────────────────────
message("\n[+] Fetching Arabidopsis orthologs...")

# Query 1: NCBI entrez ID → quinoa Ensembl gene ID
entrez_to_ensembl <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters    = "entrezgene_id",
  values     = ncbi_ids,
  mart       = quinoa_mart
)
message("Mapped ", nrow(entrez_to_ensembl), " entrez IDs to Ensembl gene IDs")

# Query 2: quinoa Ensembl gene ID → Arabidopsis ortholog
ensembl_to_orthologs <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "athaliana_eg_homolog_ensembl_gene",
    "athaliana_eg_homolog_orthology_type",
    "athaliana_eg_homolog_perc_id"
  ),
  filters    = "ensembl_gene_id",
  values     = entrez_to_ensembl$ensembl_gene_id,
  mart       = quinoa_mart
)
message("Found ", nrow(ensembl_to_orthologs), " ortholog mappings")

# Merge both queries
orthologs <- merge(entrez_to_ensembl, ensembl_to_orthologs, by="ensembl_gene_id")

message("Orthologs found: ", nrow(orthologs), " mappings for ",
        length(unique(orthologs$entrezgene_id)), " quinoa genes")

# Keep only one-to-one orthologs or best hit per quinoa gene
best_hits <- orthologs %>%
  filter(athaliana_eg_homolog_orthology_type == "ortholog_one2one" |
         athaliana_eg_homolog_perc_id > 0) %>%
  group_by(entrezgene_id) %>%
  slice_max(athaliana_eg_homolog_perc_id, n=1, with_ties=FALSE) %>%
  ungroup()

message("Best hits: ", nrow(best_hits), " quinoa → Arabidopsis mappings")

# Save ortholog table
write.table(best_hits,
            file.path(output_dir, "04_quinoa_arabidopsis_orthologs.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ── Step 5: GO enrichment with clusterProfiler ────────────────────────────────
arabidopsis_genes <- unique(best_hits$athaliana_eg_homolog_ensembl_gene)
arabidopsis_genes <- arabidopsis_genes[arabidopsis_genes != ""]

message("\n[+] Running GO enrichment on ", length(arabidopsis_genes), " Arabidopsis genes...")
message("    (note: small gene lists may yield few or no significant terms)")

go_results <- enrichGO(
  gene          = arabidopsis_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "ALL",          # BP, MF, CC combined
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

if (nrow(as.data.frame(go_results)) == 0) {
  message("[!] No significant GO terms found. Try relaxing pvalueCutoff.")
} else {
  message("Significant GO terms: ", nrow(as.data.frame(go_results)))

  # Save full GO results
  write.table(as.data.frame(go_results),
              file.path(output_dir, "04_GO_enrichment_results.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)

  # Dotplot — wrap long labels and give each ontology panel enough room
  go_plot <- go_results
  go_plot@result$Description <- stringr::str_wrap(go_plot@result$Description, width = 40)

  p_dot <- dotplot(go_plot, showCategory = 15, split = "ONTOLOGY") +
    ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free", space = "free_y") +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = 8, lineheight = 0.85),
      axis.text.x  = ggplot2::element_text(size = 8),
      strip.text   = ggplot2::element_text(size = 9, face = "bold"),
      legend.text  = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 8),
      panel.spacing = ggplot2::unit(0.6, "lines")
    )

  png(file.path(output_dir, "04_GO_dotplot.png"), width = 2800, height = 4200, res = 300)
  print(p_dot)
  dev.off()

  message("Done ", output_dir)
}
