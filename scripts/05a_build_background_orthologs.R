###############################################################
### 05a_build_background_orthologs.R
### Build a full background ortholog table from ALL expressed
### genes (rowSums >= 10), not just DEGs.
###
### Strategy (avoids biomaRt POST/405 issue):
###   Step 1: Download full quinoa gene table (no filter, GET)
###   Step 2: Filter to expressed genes locally
###   Step 3: Batch query Arabidopsis orthologs by Ensembl ID (GET)
###
### Output: results/05_cluster_go/background_orthologs_all.tsv
###############################################################

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(httr2))
suppressPackageStartupMessages(library(dplyr))

# ── Configuration ──────────────────────────────────────────────────────────────
if (!exists("data_dir"))   data_dir   <- "/home/ziqi/Projects/quinoa_raquel/rnaseq_analysis/data"
if (!exists("output_dir")) output_dir <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/05_cluster_go"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(output_dir, "background_orthologs_all.tsv")

if (file.exists(out_file)) {
  message("Background ortholog file already exists: ", out_file)
  message("Delete it to recompute. Exiting.")
  quit(save = "no")
}

# ── Patch biomaRt: replace POST with GET ──────────────────────────────────────
# biomaRt 2.66 uses httr2 POST, but Ensembl Plants returns 405.
# The BioMart REST endpoint accepts query as GET URL parameter.
.submitQueryXML_GET <- function(host, query, http_config) {
  req <- req_options(
    req_timeout(req_url_query(httr2::request(host), query = query),
                max(getOption("timeout", default = 300), 300)),
    !!!http_config
  )
  res <- httr2::req_perform(req)
  if (httr2::resp_is_error(res)) stop(httr2::resp_status(res))
  httr2::resp_body_string(res)
}
assignInNamespace(".submitQueryXML", .submitQueryXML_GET, ns = "biomaRt")
message("biomaRt patched: POST → GET")

# ── Load expressed gene IDs ────────────────────────────────────────────────────
message("Loading count matrix to get expressed gene IDs...")
counts_raw  <- read.csv(file.path(data_dir, "gene_count.csv"),
                        sep = ";", header = TRUE, stringsAsFactors = FALSE)
samples     <- read.csv(file.path(data_dir, "samples.csv"),
                        sep = ";", header = TRUE, stringsAsFactors = FALSE)
counts      <- as.matrix(counts_raw[, samples$sample_name])
mode(counts) <- "integer"
rownames(counts) <- counts_raw$gene_id
expressed   <- rownames(counts)[rowSums(counts) >= 10]
message("Expressed genes (rowSums >= 10): ", length(expressed))

# ── Connect to Ensembl Plants ──────────────────────────────────────────────────
message("Connecting to Ensembl Plants...")
plants_mart <- tryCatch(
  useEnsemblGenomes(biomart = "plants_mart", dataset = "cquinoa_eg_gene"),
  error = function(e) {
    message("Primary host failed: ", e$message)
    stop(e)
  }
)
message("Connected.")

# ── Step 1: Download full quinoa ID table (no filter) ─────────────────────────
# Two separate queries to avoid "multiple attribute pages" error.
message("\n[Step 1] Downloading full quinoa gene table (ensembl_gene_id + entrezgene_id)...")
q1_full <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  mart       = plants_mart
)
message("Total rows: ", nrow(q1_full))

# Filter to expressed genes
q1_full$entrezgene_id <- as.character(q1_full$entrezgene_id)
q1_expressed <- q1_full %>%
  filter(!is.na(entrezgene_id), entrezgene_id != "NA",
         entrezgene_id %in% expressed) %>%
  distinct()

n_expr_ens <- n_distinct(q1_expressed$ensembl_gene_id)
message("Expressed genes with Ensembl ID: ", n_expr_ens, " / ", length(expressed))

if (n_expr_ens == 0) stop("No expressed genes matched. Check gene ID format.")

# ── Step 2: Batch query Arabidopsis orthologs ──────────────────────────────────
message("\n[Step 2] Querying Arabidopsis orthologs (batch size 500)...")
ensembl_ids <- unique(q1_expressed$ensembl_gene_id)
batch_size  <- 500
n_batches   <- ceiling(length(ensembl_ids) / batch_size)
q2_results  <- vector("list", n_batches)

for (i in seq_len(n_batches)) {
  idx   <- ((i - 1) * batch_size + 1):min(i * batch_size, length(ensembl_ids))
  batch <- ensembl_ids[idx]

  if (i == 1 || i %% 5 == 0 || i == n_batches)
    message("  Batch ", i, "/", n_batches,
            " (", length(batch), " genes)")

  res <- tryCatch(
    getBM(
      attributes = c("ensembl_gene_id",
                     "athaliana_eg_homolog_ensembl_gene",
                     "athaliana_eg_homolog_orthology_type",
                     "athaliana_eg_homolog_perc_id"),
      filters    = "ensembl_gene_id",
      values     = batch,
      mart       = plants_mart
    ),
    error = function(e) {
      message("  Batch ", i, " failed: ", e$message, " — retrying in 5s")
      Sys.sleep(5)
      tryCatch(
        getBM(attributes = c("ensembl_gene_id",
                             "athaliana_eg_homolog_ensembl_gene",
                             "athaliana_eg_homolog_orthology_type",
                             "athaliana_eg_homolog_perc_id"),
              filters = "ensembl_gene_id", values = batch, mart = plants_mart),
        error = function(e2) { message("  Retry failed: ", e2$message); NULL }
      )
    }
  )

  q2_results[[i]] <- res
  Sys.sleep(0.3)
}

q2 <- bind_rows(q2_results) %>%
  filter(!is.na(athaliana_eg_homolog_ensembl_gene),
         athaliana_eg_homolog_ensembl_gene != "") %>%
  distinct()

message("Quinoa–Arabidopsis pairs: ", nrow(q2))

# ── Combine and save ───────────────────────────────────────────────────────────
background <- q1_expressed %>%
  left_join(q2, by = "ensembl_gene_id") %>%
  filter(!is.na(athaliana_eg_homolog_ensembl_gene),
         athaliana_eg_homolog_ensembl_gene != "") %>%
  distinct()

n_at  <- n_distinct(background$athaliana_eg_homolog_ensembl_gene)
n_qu  <- n_distinct(background$entrezgene_id)

message("\n=== Background set summary ===")
message("Expressed quinoa genes:             ", length(expressed))
message("  with Ensembl ID:                  ", n_expr_ens)
message("  with Arabidopsis ortholog:        ", n_qu,
        " (", round(100 * n_qu / length(expressed), 1), "%)")
message("Unique Arabidopsis background genes:", n_at)
message("Saving to: ", out_file)

write.table(background, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
message("Done.")
