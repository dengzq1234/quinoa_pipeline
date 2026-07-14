#!/bin/bash

# Usage: bash run_pipeline.sh data/input/demo_gene.txt

set -e

INPUT_FILE="$1"

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "[ERROR] Input file not found: $INPUT_FILE"
    exit 1
fi

# Get prefix from full filename without extension (demo_gene.txt → demo_gene)
BASENAME=$(basename "$INPUT_FILE")
PREFIX="${BASENAME%.*}"  # keep full prefix before ".txt"

# Define folders
RESULT_DIR="results/$PREFIX"
FASTA_DIR="fastas/$PREFIX"
mkdir -p "$RESULT_DIR" "$FASTA_DIR"

# Define output files using step numbers and prefix
STEP1_OUTPUT="$RESULT_DIR/01_${PREFIX}_uncharacterized_protein_coding_genes.tsv"
STEP2_FASTAS="$FASTA_DIR/"
STEP3_BLAST_OUT="$RESULT_DIR/03_${PREFIX}_all_blast_results.tsv"
STEP4_ANNOT_OUT="$RESULT_DIR/04_${PREFIX}_blast_with_annotations.tsv"
STEP4_GO_OUT="$RESULT_DIR/04_${PREFIX}_shinygo_gene_list.txt"

echo "[1] Filtering uncharacterized protein-coding genes..."
python scripts/01_filter_gene_list.py "$INPUT_FILE" --unchar > "$STEP1_OUTPUT"

echo "[2] Extracting FASTA sequences..."
python scripts/02_get_fasta.py "$STEP1_OUTPUT" "$STEP2_FASTAS"

echo "[3] Running DIAMOND blastp..."
bash scripts/03_blast_against_arabidopsis.sh "$STEP2_FASTAS" "$STEP3_BLAST_OUT"

echo "[4] Annotating BLAST hits and extracting gene list for GO analysis..."
python scripts/04_parse_annotation.py "$STEP3_BLAST_OUT" --annot_out "$STEP4_ANNOT_OUT" --go_out "$STEP4_GO_OUT"

echo "[✓] Pipeline finished for $PREFIX"
echo "Results:"
echo " - Step 1: $STEP1_OUTPUT"
echo " - Step 2: $STEP2_FASTAS"
echo " - Step 3: $STEP3_BLAST_OUT"
echo " - Step 4: $STEP4_ANNOT_OUT"
echo " - Step 4: $STEP4_GO_OUT"
