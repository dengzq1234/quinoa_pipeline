from Bio import Entrez, SeqIO
from Bio.Entrez import efetch, read
import pandas as pd

from time import sleep
import sys, os

# Set your email
Entrez.email = "your_email@example.com"



def fetch_protein_fasta(accession, gene_name, output_dir, max_retries=3):
    """Fetches protein sequence from NCBI, retries up to 3 times on failure, and saves as a FASTA file."""
    attempt = 0
    while attempt < max_retries:
        try:
            handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            sleep(0.5)  # Respect NCBI rate limit
            
            # Define filename
            filename = os.path.join(output_dir, f"{gene_name}.{accession}.fasta")

            # Save to FASTA file
            with open(filename, "w") as fasta_file:
                SeqIO.write(record, fasta_file, "fasta")
            
            # Print success message
            #print(f"Saved: {filename}")
            return  # Exit the function on success

        except Exception as e:
            #print(f" Attempt {attempt + 1} failed for {accession}: {e}")
            attempt += 1
            sleep(1)  # Wait before retrying
            
    print(f"{accession}")

# Main function to process gene IDs
def main():
    # Read Gene IDs from file
    gene_list = sys.argv[1]
    # Directory to save FASTA files
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)
    # get protein coding df
    protein_coding_df = pd.read_csv(gene_list, sep='\t')

    
    # Process each row in the DataFrame
    for index, row in protein_coding_df.iterrows():
        if pd.notna(row["refseq_peptide"]):  # Check if refseq_peptide is not NaN
            peptide_accessions = row["refseq_peptide"].split(";")  # Split if multiple accessions
            for accession in peptide_accessions:
                fetch_protein_fasta(accession.strip(), row["gene_name"], output_dir)  # Fetch and save

if __name__ == '__main__':
    main()