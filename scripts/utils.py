# Python utils

# Function to remotely run BLAST using `QBLAST`
# Based on: https://www.tutorialspoint.com/biopython/biopython_overview_of_blast.htm
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

def run_qblast(query_sequence, blast_program, database):
    result_handle = NCBIWWW.qblast(blast_program, database, query_sequence)

    blast_records = NCBIXML.parse(result_handle)

    table = []
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                accession = alignment.accession
                taxid = accession.split("|")[1]
                taxonomy = get_taxonomy_from_taxid(taxid)  # Function to get taxonomy based on NCBI TaxID
                table.append((blast_record.query_id, taxid, taxonomy))

    return table

# Function to get taxonomy based on NCBI taxid
def get_taxonomy_from_taxid(taxid):
    # Implement your logic to fetch taxonomy information based on the taxid.
    # Return a dummy taxonomy string for demonstration purposes.
    return f"Dummy Taxonomy for TaxID {taxid}"
