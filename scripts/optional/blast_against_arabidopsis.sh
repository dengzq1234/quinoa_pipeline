#!/bin/bash

# scripts/03_blast_against_arabidopsis.sh

# Usage:
#   bash scripts/03_blast_against_arabidopsis.sh fastas/ results/output/all_blast_results.tsv


# Input arguments
FASTA_DIR="$1"
MERGED_OUT="$2"

# Configuration
REF_DIR="data/reference"
DB="$REF_DIR/arabidopsis_protein_db.dmnd"
EVALUE="1e-5"
THREADS="4"
COVERAGE="50"
IDENTITY="30"

# Create reference directory if needed
mkdir -p "$REF_DIR"

# Step 0: Download Arabidopsis protein FASTA from Ensembl Plants
if [ ! -f "$REF_DIR/Arabidopsis_thaliana.TAIR10.pep.all.fa" ]; then
    echo "[+] Downloading Arabidopsis protein FASTA..."
    wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz -P "$REF_DIR"
    gunzip -f "$REF_DIR/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"
else
    echo "[✓] Arabidopsis FASTA already exists."
fi

# Step 1: Build DIAMOND DB if not exists
if [ ! -f "$DB" ]; then
    echo "[+] Building DIAMOND database..."
    diamond makedb --in "$REF_DIR/Arabidopsis_thaliana.TAIR10.pep.all.fa" -d "$REF_DIR/arabidopsis_protein_db"
else
    echo "[✓] DIAMOND database already exists."
fi

# Step 2: Run DIAMOND BLASTP
echo "[+] Running DIAMOND BLASTP for all FASTA files in $FASTA_DIR..."
mkdir -p "$(dirname "$MERGED_OUT")"
echo -e "qseqid\tsseqid\tpident\tlength\tevalue\tbitscore\tqlen\tslen\tqcovhsp\tquery_file" > "$MERGED_OUT"

for fasta in "$FASTA_DIR"/*.fasta; do
    BASENAME=$(basename "$fasta" .fasta)
    echo "[+] Processing $BASENAME..."

    diamond blastp \
        --query "$fasta" \
        --db "$DB" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore qlen slen qcovhsp \
        --evalue "$EVALUE" \
        --threads "$THREADS" \
        --id "$IDENTITY" \
        --subject-cover "$COVERAGE" \
        --query-cover "$COVERAGE" \
        | awk -v file="$BASENAME" -v OFS='\t' '{print $0, file}' >> "$MERGED_OUT"

    echo "[✓] Results for $BASENAME appended."
done

echo "[✓] BLAST pipeline complete: $MERGED_OUT"