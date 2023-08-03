# Rules for running qBLAST on the ASV fasta file
#
# Example call: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline --cores 1 -rpn
# Example call with conda-prefix-path: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline --cores $CORES -rpn


##############################
# default
rule blast:
    input:
        os.path.join(RESULTS_DIR, "blast/taxonomy_table.txt")
    output:
        touch("status/blast.done")


# Running BLAST
rule run_qblast_nt:
    input:
        query_fasta=os.path.join(FASTA)
    output:
        table=os.path.join(RESULTS_DIR, "blast/taxonomy_table.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/blast_nt.log")
    params:
        blast_program=config["blast"]["program"],
        blast_database=config["blast"]["db"],
        qblast=os.path.join(SRC_DIR, "run_qblast.py")
    conda:
        os.path.join(ENV_DIR, "biopython.yaml")
    message:
        "Running BLAST on ASV FASTA input"
    shell:
        "(date && "
        "{params.qblast} {input.query_fasta} {params.blast_program} {params.blast_database} {output.table} &&"
        "date) &> {log}"

