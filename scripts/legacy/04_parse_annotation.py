import pandas as pd
import requests
import time
import argparse
import os

def fetch_uniprot_annotations(ids):
    """
    Fetch annotations from UniProt for given gene IDs and print results per ID.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    annotations = {}

    for uniprot_id in ids:
        gene_symbol = uniprot_id.split('.')[0]  # Strip isoform if any
        query = f"gene:{gene_symbol} AND organism_id:3702"
        params = {
            "query": query,
            "format": "tsv",
            "fields": "accession,gene_names,protein_name,go,cc_domain,cc_function"
        }
        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            lines = response.text.strip().split("\n")
            if len(lines) > 1:
                headers = lines[0].split("\t")
                values = lines[1].split("\t")
                annotation = dict(zip(headers, values))
                annotations[uniprot_id] = annotation
                print(f"[✓] Annotated {uniprot_id}")
            else:
                print(f"[!] No results for {uniprot_id}")
        except Exception as err:
            print(f"[X] Error fetching {uniprot_id}: {err}")
        time.sleep(0.3)  # Respect UniProt's rate limit

    return annotations

def main():
    parser = argparse.ArgumentParser(description="Merge UniProt annotations with BLAST results and export ShinyGO gene list.")
    parser.add_argument("blast_input", help="Path to input BLAST result file (TSV format).")
    parser.add_argument("--annot_out", default="results/blast/blast_with_annotations.tsv", help="Path to output merged annotation file.")
    parser.add_argument("--go_out", default="results/blast/shinygo_gene_list.txt", help="Path to ShinyGO gene list output.")

    args = parser.parse_args()
    os.makedirs(os.path.dirname(args.annot_out), exist_ok=True)
    os.makedirs(os.path.dirname(args.go_out), exist_ok=True)

    # Load BLAST results
    blast_df = pd.read_csv(args.blast_input, sep="\t")

    # Extract unique base gene names from sseqid
    blast_df["gene_base"] = blast_df["sseqid"].str.split('.').str[0]
    unique_genes = blast_df["gene_base"].unique().tolist()

    # Fetch annotations
    annotations = fetch_uniprot_annotations(unique_genes)

    # Convert annotation dict to DataFrame
    ann_df = pd.DataFrame.from_dict(annotations, orient="index").reset_index()
    ann_df.rename(columns={"index": "query_gene"}, inplace=True)

    # Merge using query_gene as the key
    merged_df = pd.merge(blast_df, ann_df, left_on="gene_base", right_on="query_gene", how="left")
    merged_df.drop(columns=["gene_base", "query_gene"], inplace=True)

    # Save the full merged file
    merged_df.to_csv(args.annot_out, sep="\t", index=False)

    # Select best hit per quinoa gene based on highest bitscore
    best_hits = merged_df.sort_values("bitscore", ascending=False).drop_duplicates("query_file")

    # Extract Arabidopsis genes
    if "gene_names" in best_hits.columns and best_hits["gene_names"].notnull().any():
        arabidopsis_genes = best_hits["gene_names"].str.split().str[0]
    else:
        arabidopsis_genes = best_hits["sseqid"].str.split(".").str[0]

    # Save list for ShinyGO
    arabidopsis_genes.drop_duplicates().to_csv(args.go_out, index=False, header=False)

    print(f"[✓] Saved: {args.annot_out}")
    print(f"[✓] Saved: {args.go_out}")

if __name__ == "__main__":
    main()
