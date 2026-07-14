# Uncharacterized Gene Annotation and Functional Inference Pipeline

This pipeline takes a list of quinoa NCBI gene IDs, filters for uncharacterized protein-coding genes, maps them to Arabidopsis orthologs via Ensembl Plants (biomaRt), and runs GO enrichment analysis using clusterProfiler — all within R.

## Requirements

Run once to install all required R packages:

```r
source("scripts/00_install_packages.R")
```

Packages installed: `rentrez`, `xml2`, `biomaRt`, `clusterProfiler`, `org.At.tair.db`, `dplyr`

## Project Structure

```
├── data/
│   ├── input/          # Input gene ID list (one NCBI gene ID per line)
│   └── reference/      # Arabidopsis DIAMOND database (used by optional BLAST step)
├── results/            # Output files per run
├── scripts/
│   ├── 00_install_packages.R       # one-time setup
│   ├── 01_filter_gene_list.R       # core: filter uncharacterized protein-coding genes
│   ├── 02_annotate_and_go.R        # core: biomaRt orthology + GO enrichment
│   ├── optional/
│   │   ├── get_fasta.R             # fetch protein FASTAs from NCBI (if needed)
│   │   └── blast_against_arabidopsis.sh  # DIAMOND BLAST vs Arabidopsis (if needed)
│   └── legacy/                     # original Python pipeline (archived)
│       ├── 01_filter_gene_list.py
│       ├── 02_get_fasta.py
│       ├── 04_parse_annotation.py
│       ├── 05_shinygo_analysis.md
│       └── run_pipeline.sh
└── readme.md
```

## Core Pipeline (two steps)

### Step 1 — Filter uncharacterized protein-coding genes

Edit the configuration block at the top of `scripts/01_filter_gene_list.R`:

```r
input_file  <- "data/input/your_gene_list.txt"   # one NCBI gene ID per line
output_file <- "results/your_prefix/01_filtered_genes.tsv"
only_unchar <- TRUE   # FALSE to keep all protein-coding genes
```

Run:
```r
source("scripts/01_filter_gene_list.R")
```

**Output:** TSV with columns `gene_id, gene_name, gene_desc, gene_type, refseq_genomic, refseq_mrna, refseq_peptide, related_acc`

---

### Step 2 — Arabidopsis ortholog mapping + GO enrichment

Edit the configuration block at the top of `scripts/02_annotate_and_go.R`:

```r
input_file <- "data/input/your_gene_list.txt"         # raw gene ID list
output_dir <- "results/your_prefix/"
```

Run:
```r
source("scripts/02_annotate_and_go.R")
```

**Outputs:**
- `04_quinoa_arabidopsis_orthologs.tsv` — quinoa gene → Arabidopsis ortholog mapping
- `04_GO_enrichment_results.tsv` — full GO enrichment table (BP, MF, CC)
- `04_GO_dotplot.png` — dotplot of top 20 enriched GO terms per ontology
- `04_GO_interpretation.html` — biological interpretation of results

## Optional Steps

These are only needed if protein FASTA sequences are required for downstream analyses (e.g. phylogenetics, structural prediction):

```r
source("scripts/optional/get_fasta.R")              # fetch FASTAs from NCBI
bash scripts/optional/blast_against_arabidopsis.sh   # DIAMOND BLAST
```

## Notes

- Input file should contain one NCBI Entrez gene ID per line, no header
- `02_annotate_and_go.R` can be run directly on the raw gene list — Step 1 is not required for GO enrichment
- Orthology is inferred using curated Ensembl Plants mappings (biomaRt), not BLAST
- GO enrichment uses Arabidopsis annotations (`org.At.tair.db`) since quinoa lacks a comprehensive GO database
- For best results, use a DEG list from DESeq2 as input rather than a background gene list
