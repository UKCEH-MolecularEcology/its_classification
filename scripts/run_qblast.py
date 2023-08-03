#!/usr/bin/env python

###########
# MODULES #
###########
# importing relevant modules
import os, sys
import logging
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO


###########
# LOGGING #
###########
# logger
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG,
    format='[%(asctime)s] %(name)s %(levelname)s: %(message)s'
)
logger = logging.getLogger(__file__)


##########
# qBLAST #
#########
def run_qblast(query_sequence, blast_program, database, output_table):
    print(f"Running {blast_program} BLAST against {database}...")
    result_handle = NCBIWWW.qblast(blast_program, database, query_sequence)
    print("BLAST search completed!")

    blast_records = NCBIXML.parse(result_handle)

    table = []
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                accession = alignment.accession
                taxid = accession
                taxonomy = get_taxonomy_from_taxid(taxid)  # Function to get taxonomy based on NCBI TaxID
                table.append((blast_record.query, taxid, taxonomy))

    with open(output_table, "w") as table_file:
        table_file.write("Query_ID\tTaxID\tTaxonomy\n")
        for row in table:
            table_file.write("\t".join(row) + "\n")

    print(f"Taxonomy table written to {output_table}")

def get_taxonomy_from_taxid(taxid):
    # Implement your logic to fetch taxonomy information based on the taxid.
    # Return a dummy taxonomy string for demonstration purposes.
    return f"Dummy Taxonomy for TaxID {taxid}"

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python run_qblast.py <input_fasta> <blast_program> <blast_database> <output_table>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    blast_program = sys.argv[2]
    blast_database = sys.argv[3]
    output_table = sys.argv[4]

    print(f"Reading input FASTA file: {input_fasta}...")
    for record in SeqIO.parse(input_fasta, "fasta"):
        query_sequence = record.seq
        print(f"Processing query sequence: {record.id}...")
        run_qblast(query_sequence, blast_program, blast_database, output_table)

