# Rules to run the RDP classifier for ITS sequences against the ASV FASTA sequences
#
# Example call: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline --cores 1 -rpn
# Example call with conda-prefix-path: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline --cores $CORES -rpn


##############################
# default
rule rdp_classifier_all:
    input:
        expand(os.path.join(RESULTS_DIR, "taxonomy/assigned_taxonomy_filtered_{confidence}.txt"), confidence=CONFIDENCE)
    output:
        touch("status/rdp_classifier.done")


# Downloading the database
rule download_db:
    output:
        expand(os.path.join(RESULTS_DIR, "db/Gweon-{amplicon}-20200325-rdp-trained/Gweon-{amplicon}.properties"), amplicon=AMPLICON_TYPE)
    log:
        os.path.join(RESULTS_DIR, "logs/download_db.log")
    message:
        "Downloading the database required for classification"
    params:
        url=config["database"][AMPLICON_TYPE]
    wildcard_constraints:
        amplicon="|".join(AMPLICON_TYPE)
    shell:
        "(date && "
        "wget -P $(dirname $(dirname {output})) {params.url} &&"
        "cd $(dirname $(dirname {output})) && tar -xzvf $(basename {params.url}) &&"
        "date) &> {log}"

# Running the classifier
rule rdp_classifier:
    input:
        fa=os.path.join(FASTA),
        db=rules.download_db.output[0]
    output:
        os.path.join(RESULTS_DIR, "taxonomy/assigned_taxonomy_prelim.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/rdp_classifier_prelim.log")
    conda:
        os.path.join(ENV_DIR, "rdp_classifier.yaml")
    params:
        mem=config["mem"]
    message:
        "Running the RDP classifier on the ASV input"
    shell:
        "(date && "
        "rdp_classifier -Xmx{params.mem} classify -t {input.db} -o {output} {input.fa} && "
        "date) &> {log}"

# Reformat the taxonomy
rule reformat:
    input:
        rules.rdp_classifier.output[0]
    output:
        os.path.join(RESULTS_DIR, "taxonomy/assigned_taxonomy_filtered_{confidence}.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/rdp_classifier_{confidence}.log")
    params:
        SRC=os.path.join(SRC_DIR, "honeypi_reformatAssignedTaxonomy.py")
    message:
        "Applying the confidence score of {wildcards.confidence} to the taxonomy output"
    shell:
        "(date && "
        "python3 {params.SRC} -i {input} -o {output} -c {wildcards.confidence} && "
        "date) &> {log}"

