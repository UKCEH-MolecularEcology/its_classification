# its_classification
- Snakemake workflow to run [RDP taxonomic classification](https://github.com/rdpstaff/classifier) on ITS FASTA sequences.
- based on [honeypi](https://github.com/hsgweon/honeypi/blob/master/bin/honeypi)

## Requirements
- ASV FASTA from DADA2 or Qiime2.
- Pipelines
  - [snakemake](https://snakemake.github.io/) (required)
- Software
  - Conda (can be used together w/ `snakemake`) (required)
    - [Conda](https://docs.conda.io/projects/conda/en/latest/index.html)
    - [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
    - [Anaconda](https://anaconda.org/)

## How to run
### Install conda
- Do the following to install `conda` on neohuxley.
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

### Clone the repository
- Download the GitHub repository.
```bash
git clone git@github.com:UKCEH-MolecularEcology/its_classification.git
cd its_classification
```

### Create the `snakemake` environment
```bash
# create the snakemake conda environment
conda env create -f requirements.yaml

# Activate the environment
conda activate its_classifier
```

### Set up the run
- Edit the `config/config.yaml`
- Following items in the `config.yaml` file need *USER input*
  - *amplicon_type:* choose from **ITS1** or **ITS2**
  - *fasta:* provide the path to the ASV FASTA file
  - `data_dir` and `results_dir` need to be provided
  - *confidence:* adjust the values as required; between 0 to 1

## Running the workflow
```bash
# dry-run to check if everything is in order
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 24 -rpn

# full run
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 24 -rp
```
