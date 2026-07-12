library(rentrez)
library(xml2)

# ── Configuration ──────────────────────────────────────────────────────────────
input_file  <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/data/input/demo_gene.txt"
output_file <- "/home/ziqi/Projects/quinoa_raquel/quinoa_pipeline/results/demo_gene/01_demo_gene_filtered_genes.tsv"
only_unchar <- TRUE   # TRUE = keep only uncharacterized protein-coding genes

# ── Helpers ────────────────────────────────────────────────────────────────────

extract_refseq <- function(doc) {
  refseq_node <- xml_find_first(doc,
    "//Gene-commentary[Gene-commentary_heading[text()='NCBI Reference Sequences (RefSeq)']]")

  if (is.na(refseq_node)) return(list(genomic=character(), mrna=character(), peptide=character()))

  list(
    genomic = xml_text(xml_find_all(refseq_node,
      ".//Gene-commentary[Gene-commentary_type/@value='genomic']/Gene-commentary_accession")),
    mrna    = xml_text(xml_find_all(refseq_node,
      ".//Gene-commentary[Gene-commentary_type/@value='mRNA']/Gene-commentary_accession")),
    peptide = xml_text(xml_find_all(refseq_node,
      ".//Gene-commentary[Gene-commentary_type/@value='peptide']/Gene-commentary_accession"))
  )
}

extract_genbank <- function(doc) {
  related_node <- xml_find_first(doc,
    "//Gene-commentary[Gene-commentary_heading[text()='Related Sequences']]")
  if (is.na(related_node)) return(character())
  xml_text(xml_find_all(related_node, ".//Gene-commentary_accession"))
}

fetch_gene_info <- function(gene_id) {
  tryCatch({
    raw  <- entrez_fetch(db="gene", id=gene_id, rettype="xml", retmode="xml")
    Sys.sleep(0.4)
    doc  <- read_xml(raw)

    gene_name <- xml_text(xml_find_first(doc, "//Gene-ref_locus"))
    gene_desc <- xml_text(xml_find_first(doc, "//Prot-ref_desc"))
    gene_type <- xml_attr(xml_find_first(doc, "//Entrezgene_type"), "value")

    if (is.na(gene_name)) gene_name <- ""
    if (is.na(gene_desc)) gene_desc <- "uncharacterized"
    if (is.na(gene_type)) gene_type <- "unknown"

    refseq  <- extract_refseq(doc)
    genbank <- extract_genbank(doc)

    list(
      gene_id        = gene_id,
      gene_name      = gene_name,
      gene_desc      = gene_desc,
      gene_type      = gene_type,
      refseq_genomic = paste(refseq$genomic, collapse=";"),
      refseq_mrna    = paste(refseq$mrna,    collapse=";"),
      refseq_peptide = paste(refseq$peptide, collapse=";"),
      related_acc    = paste(genbank,         collapse=";")
    )
  }, error = function(e) {
    message("Error for gene ", gene_id, ": ", e$message)
    list(gene_id=gene_id, gene_name="", gene_desc="", gene_type="",
         refseq_genomic="", refseq_mrna="", refseq_peptide="", related_acc="")
  })
}

# ── Main ───────────────────────────────────────────────────────────────────────

gene_ids <- trimws(readLines(input_file, warn=FALSE))
gene_ids <- gene_ids[nzchar(gene_ids)]

results <- list()

for (i in seq_along(gene_ids)) {
  message("[", i, "/", length(gene_ids), "] Fetching gene ", gene_ids[i], "...")
  info <- fetch_gene_info(gene_ids[i])

  if (info$gene_type != "protein-coding") next
  if (only_unchar && !grepl("uncharacterize", tolower(info$gene_desc))) next

  results[[length(results) + 1]] <- info
}

out_df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors=FALSE))

dir.create(dirname(output_file), recursive=TRUE, showWarnings=FALSE)
write.table(out_df, output_file, sep="\t", row.names=FALSE, quote=FALSE)

message("[✓] Done. ", nrow(out_df), " genes written to ", output_file)
