##############################
# MODULES
import os, re
import glob
import pandas as pd


##############################
# CONFIG
# can be overwritten by using --configfile <path to config> when calling snakemake
# configfile: "config/config.yaml"


##############################
# Paths
SRC_DIR = srcdir("../scripts")
ENV_DIR = srcdir("../envs")


##############################
# default executable for snakemake
shell.executable("bash")


##############################
# working directory
workdir:
    config["data_dir"]

##############################
# Relevant directories
DATA_DIR = config["data_dir"]
FASTA = config["fasta"]
RESULTS_DIR = config["results_dir"]
CONFIDENCE = config["confidence"]
AMPLICON_TYPE = config["amplicon_type"]

##############################
# Steps
STEPS = config["steps"]
