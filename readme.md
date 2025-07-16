# 🔍 Uncharacterized Gene Annotation and Functional Inference Pipeline

This pipeline allows you to analyze a list of gene identifiers (e.g., from Novogene), filter for uncharacterized protein-coding genes, retrieve protein sequences, perform homology searches against the Arabidopsis proteome using DIAMOND, annotate hits with UniProt metadata, and prepare gene lists for functional enrichment analysis using g:Profiler.

## Install Environment

install conda
```
# Download the Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer
bash Miniconda3-latest-Linux-x86_64.sh

# Follow the prompts and restart your terminal or run:
source ~/.bashrc
```

Install environment
```
conda create -n quinoa python=3.10 -y
conda create -n quinoa_blast_env python=3.10 diamond pandas biopython requests  -y
conda activate quinoa_blast_env
```

## Project Structure
```
├── data/
│ ├── input/ # Input gene list
│ └── reference/ # for arabidosis sequence database/
├── fastas/ # Extracted FASTA files (protein) from target gene
├── results/ # Final annotation results and GO list
├── scripts/ # Python and Bash scripts used in pipeline
├── README.md # This file
```

## Pipeline Usage to get annotations
Place your gene list in the folder: data/input/, e.g., demo_gene.txt.

Run the pipeline using the provided bash script:

```
bash run_pipeline.sh data/input/demo_gene.txt
```

A folder will be automatically created under results/ using the input file name as the prefix (e.g., results/demo_gene/), containing the following files:

- `01_<prefix>_filtered_genes.tsv` — filtered table of uncharacterized protein-coding genes.

- `fastas/<prefix>/` — folder of FASTA files for each gene.

- `03_<prefix>_blast_results.tsv` — DIAMOND results against Arabidopsis proteome.

- `04_<prefix>_blast_with_annotations.tsv` — merged BLAST hits with UniProt annotation.

- `04_<prefix>_shinygo_gene_list.txt` — gene list for downstream GO analysis.

## Functional Enrichment Analysis
After running the pipeline, follow the step-by-step guide in:
```
scripts/05_shinygo_analysis.md
```
This guide explains how to use g:Profiler (https://biit.cs.ut.ee/gprofiler/gost) to analyze the gene list found in:

```
results/<prefix>/04_shinygo_gene_list.txt
```

Set the organism to Arabidopsis thaliana and export results including enriched GO terms and pathways.

Examples output
https://biit.cs.ut.ee/gplink/l/a2Uyjp24XTu

## Notes
Output folders are automatically created in `results/` using the prefix of your input file.

You can enable filtering for uncharacterized genes by adding the `--unchar` flag in Step 1.
