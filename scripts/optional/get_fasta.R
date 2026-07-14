library(rentrez)

# ── Configuration ──────────────────────────────────────────────────────────────
input_file <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/demo_gene/01_demo_gene_filtered_genes.tsv"
output_dir <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/fastas/demo_gene/"

# ── Main ───────────────────────────────────────────────────────────────────────
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

gene_df <- read.table(input_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

for (i in seq_len(nrow(gene_df))) {
  row <- gene_df[i, ]

  if (is.na(row$refseq_peptide) || row$refseq_peptide == "") next

  accessions <- trimws(strsplit(row$refseq_peptide, ";")[[1]])

  for (acc in accessions) {
    filename <- file.path(output_dir, paste0(row$gene_name, ".", acc, ".fasta"))

    success <- FALSE
    for (attempt in 1:3) {
      tryCatch({
        fasta <- entrez_fetch(db="protein", id=acc, rettype="fasta", retmode="text")
        Sys.sleep(0.5)
        writeLines(fasta, filename)
        success <- TRUE
      }, error = function(e) {
        message("Attempt ", attempt, " failed for ", acc, ": ", e$message)
        Sys.sleep(1)
      })
      if (success) break
    }

    if (!success) message("Failed to fetch: ", acc)
  }

  message("[", i, "/", nrow(gene_df), "] ", row$gene_name, " done")
}

message("[✓] FASTAs saved to ", output_dir)
