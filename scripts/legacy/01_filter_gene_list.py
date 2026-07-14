from Bio import Entrez
from Bio.Entrez import efetch, read
from time import sleep
import sys
import argparse

# Set your email
Entrez.email = "your_email@example.com"

# Function to extract RefSeq accessions
def extract_refseq_accessions(commentaries):
    refseq_info = {
        "genomic": [],
        "mRNA": [],
        "peptide": []
    }

    if isinstance(commentaries, list):
        for commentary in commentaries:
            extracted_info = extract_refseq_accessions(commentary)
            for key in refseq_info:
                refseq_info[key].extend(extracted_info[key])

    elif isinstance(commentaries, dict):
        if commentaries.get("Gene-commentary_heading") == "NCBI Reference Sequences (RefSeq)":
            if "Gene-commentary_comment" in commentaries:
                for sub_commentary in commentaries["Gene-commentary_comment"]:
                    if "Gene-commentary_products" in sub_commentary:
                        for product in sub_commentary["Gene-commentary_products"]:
                            if "Gene-commentary_accession" in product:
                                acc_type = str(product.get("Gene-commentary_type", ""))
                                accession = product["Gene-commentary_accession"]

                                if acc_type == "1":  # Genomic
                                    refseq_info["genomic"].append(accession)

                                if "Gene-commentary_products" in product:
                                    for sub_product in product["Gene-commentary_products"]:
                                        acc_type = str(sub_product.get("Gene-commentary_type", ""))
                                        accession = sub_product.get("Gene-commentary_accession", "")

                                        if acc_type == "3":  # mRNA
                                            refseq_info["mRNA"].append(accession)

                                        if "Gene-commentary_products" in sub_product:
                                            for sub_product2 in sub_product["Gene-commentary_products"]:
                                                acc_type = str(sub_product2.get("Gene-commentary_type", ""))
                                                accession = sub_product2.get("Gene-commentary_accession", "")
                                                if acc_type == "8":  # Peptide
                                                    refseq_info["peptide"].append(accession)
    return refseq_info

# Function to extract GenBank accessions
def extract_genbank_accessions(commentaries):
    genbank_accessions = []
    if isinstance(commentaries, list):
        for commentary in commentaries:
            genbank_accessions.extend(extract_genbank_accessions(commentary))
    elif isinstance(commentaries, dict):
        if commentaries.get("Gene-commentary_heading") == "Related Sequences":
            if "Gene-commentary_products" in commentaries:
                for product in commentaries["Gene-commentary_products"]:
                    if "Gene-commentary_accession" in product:
                        genbank_accessions.append(product["Gene-commentary_accession"])
    return genbank_accessions

# Function to fetch gene information
def fetch_gene_info(gene_ids, only_unchar=False):
    gene_info_list = []
    for gene_id in gene_ids:
        try:
            handle = efetch(db="gene", id=gene_id, rettype="xml", retmode="text")
            records = read(handle)
            handle.close()
            sleep(0.4)  # To respect NCBI's rate limit

            gene_info = {
                'gene_id': gene_id,
                'gene_name': "",
                'gene_desc': "uncharacterized",
                'gene_type': "unknown",
                'refseq_genomic': [],
                'refseq_mrna': [],
                'refseq_peptide': [],
                'genebank_accessions': [],
            }

            # Extract gene name
            if 'Entrezgene_gene' in records[0]:
                gene_ref = records[0]['Entrezgene_gene'].get('Gene-ref', {})
                gene_info['gene_name'] = gene_ref.get('Gene-ref_locus', "")

            # Extract gene description
            if 'Entrezgene_prot' in records[0]:
                prot_ref = records[0]['Entrezgene_prot'].get('Prot-ref', {})
                gene_info['gene_desc'] = prot_ref.get('Prot-ref_desc', "uncharacterized")

            # Extract gene type
            if 'Entrezgene_type' in records[0]:
                # gene_type_map = {
                #     1: "unknown",
                #     2: "protein-coding",
                #     3: "rRNA",
                #     4: "tRNA",
                #     5: "snRNA",
                #     6: "scRNA",
                #     7: "snoRNA",
                #     8: "pseudo",
                #     9: "non-coding",
                #     10: "other"
                # }
                gene_info['gene_type'] = records[0]['Entrezgene_type'].attributes.get("value", "unknown")

            # Extract GenBank and RefSeq accessions
            if 'Entrezgene_comments' in records[0]:
                comments = records[0]['Entrezgene_comments']
                gene_info['genebank_accessions'] = extract_genbank_accessions(comments)

                refseq_data = extract_refseq_accessions(comments)
                gene_info['refseq_genomic'] = refseq_data["genomic"]
                gene_info['refseq_mrna'] = refseq_data["mRNA"]
                gene_info['refseq_peptide'] = refseq_data["peptide"]

            gene_info_list.append(gene_info)

        except Exception as e:
            sys.stderr.write(f"Error fetching data for Gene ID {gene_id}: {e}\n")
            gene_info_list.append({
                'gene_id': gene_id,
                'gene_name': "",
                'gene_desc': "",
                'gene_type': "",
                'refseq_genomic': [],
                'refseq_mrna': [],
                'refseq_peptide': [],
                'genebank_accessions': [],
            })

        # Output data
        if gene_info["gene_type"] != "protein-coding":
            continue
        if only_unchar and "uncharacterize" not in gene_info["gene_desc"].lower():
            continue

        print(f"{gene_info['gene_id']}\t{gene_info['gene_name']}\t{gene_info['gene_desc']}\t{gene_info['gene_type']}\t"
                f"{';'.join(gene_info['refseq_genomic'])}\t{';'.join(gene_info['refseq_mrna'])}\t"
                f"{';'.join(gene_info['refseq_peptide'])}\t{';'.join(gene_info['genebank_accessions'])}")

    return gene_info_list

# Main function to process gene IDs
def main():
    parser = argparse.ArgumentParser(description="Fetch gene metadata from NCBI")
    parser.add_argument("input_file", help="File with gene IDs (one per line)")
    parser.add_argument("--unchar", action="store_true", help="Only return uncharacterized protein-coding genes")
    args = parser.parse_args()

    # Read Gene IDs from file
    with open(args.input_file, 'r') as f:
        gene_ids = f.read().splitlines()

    # Print header with new "gene_type" column
    print("gene_id\tgene_name\tgene_desc\tgene_type\trefseq_genomic\trefseq_mrna\trefseq_peptide\trelated_acc")

    # Fetch gene information
    fetch_gene_info(gene_ids, only_unchar=args.unchar)

if __name__ == '__main__':
    main()
