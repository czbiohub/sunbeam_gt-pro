#
# Sunbeam: an iridescent HTS pipeline
#
# Author: Erik Clarke <ecl@mail.med.upenn.edu>
# Created: 2016-04-28
#

import os
import re
import sys
import yaml
import configparser

from pprint import pprint
from pathlib import Path, PurePath

from snakemake.utils import update_config, listfiles
from snakemake.exceptions import WorkflowError

from sunbeamlib import load_sample_list, read_seq_ids
from sunbeamlib.config import *
from sunbeamlib.reports import *

import collections


## 20190502: we don't want to keep extra unnecessary files
def cz_load_sample_list(samplelist_fp):
  Samples = collections.defaultdict()
  with open(str(samplelist_fp)) as f:
    for line in f:
      kvs = line.rstrip().split(',')
      Samples[kvs[0]] = {'1': kvs[1], '2': kvs[2]}
  return Samples

# Disallow slashes in our sample names during Snakemake's wildcard evaluation.
# Slashes should always be interpreted as directory separators.
wildcard_constraints:
  sample="[^/]+"

sunbeam_dir = ""
try:
    sunbeam_dir = os.environ["SUNBEAM_DIR"]
except KeyError:
    raise SystemExit(
        "$SUNBEAM_DIR environment variable not defined. Are you sure you're "
        "running this from the Sunbeam conda env?")

# Load extensions
sbxs = list(listfiles(sunbeam_dir+"/extensions/{sbx_folder}/{sbx}.rules"))
for sbx in sbxs:
    sys.stderr.write("Found extension {sbx} in folder {sbx_folder}\n".format(**sbx[1]))

# Setting up config files and samples
Cfg = check_config(config)
Blastdbs = process_databases(Cfg['blastdbs'])

#Samples = load_sample_list(Cfg['all']['samplelist_fp'], Cfg['all']['paired_end'], Cfg['all']['download_reads'], Cfg["all"]['root']/Cfg['all']['output_fp'])
Samples = cz_load_sample_list(Cfg['all']['samplelist_fp'])
Pairs = ['1', '2'] if Cfg['all']['paired_end'] else ['1']
#print(Samples)

# Collect host (contaminant) genomes
sys.stderr.write("Collecting host/contaminant genomes... ")
if Cfg['qc']['host_fp'] == Cfg['all']['root']:
    HostGenomeFiles = []
else:
    HostGenomeFiles = [f for f in Cfg['qc']['host_fp'].glob('*.fna')]
    if not HostGenomeFiles:
        sys.stderr.write(
            "\n\nWARNING: No files detected in host genomes folder ({}). "
            "If this is not intentional, make sure all files end in "
            ".fasta and the folder is specified correctly.\n\n".format(
                Cfg['qc']['host_fp']
            ))
HostGenomes = {Path(g.name).stem: read_seq_ids(Cfg['qc']['host_fp'] / g) for g in HostGenomeFiles}
sys.stderr.write("done.\n")

sys.stderr.write("Collecting target genomes... ")
if Cfg['mapping']['genomes_fp'] == Cfg['all']['root']:
    GenomeFiles = []
    GenomeSegments = {}
else:
    GenomeFiles = [f for f in Cfg['mapping']['genomes_fp'].glob('*.fna')]
    GenomeSegments = {PurePath(g.name).stem: read_seq_ids(Cfg['mapping']['genomes_fp'] / g) for g in GenomeFiles}
sys.stderr.write("done.\n")

# ---- Change your workdir to output_fp
workdir: str(Cfg['all']['output_fp'])

# ---- Set up output paths for the various steps
DOWNLOAD_FP = output_subdir(Cfg, 'download')
QC_FP = output_subdir(Cfg, 'qc')
ASSEMBLY_FP = output_subdir(Cfg, 'assembly')
ANNOTATION_FP = output_subdir(Cfg, 'annotation')
CLASSIFY_FP = output_subdir(Cfg, 'classify')
MAPPING_FP = output_subdir(Cfg, 'mapping')


# ---- Download rules
if Cfg['all']['download_reads']:
	include: "rules/download/download.rules"


# ---- Targets rules
include: "rules/targets/targets.rules"


# ---- Quality control rules
#include: "rules/qc/qc.rules"
#include: "rules/qc/decontaminate.rules"


# ---- Assembly rules
include: "rules/assembly/assembly.rules"
include: "rules/assembly/coverage.rules"


# ---- Contig annotation rules
include: "rules/annotation/annotation.rules"
include: "rules/annotation/blast.rules"
include: "rules/annotation/orf.rules"


# ---- Classifier rules
include: "rules/classify/kraken.rules"


# ---- Mapping rules
include: "rules/mapping/mapping.rules"


# ---- Reports rules
include: "rules/reports/reports.rules"

#for sbx_path, wildcards in sbxs:
#    include: sbx_path
        

# ---- Rule all: run all targets
rule all:
    input: TARGET_ALL

rule samples:
    message: "Samples to be processed:"
    run:
        [print(sample) for sample in sorted(list(Samples.keys()))]

include: "rules/gtpro.rules"

#rule _all_filtered_cov:
#    input:
#       expand(str(MAPPING_FP/'filtered'/'{genome}'/'{sample}.depth'),genome=GenomeSegments.keys(), sample=Samples.keys())

#rule filtered_coverage:
#    input:
#        bam = str(MAPPING_FP/'filtered'/'{genome}'/'{sample}.bam'),
#        bai = str(MAPPING_FP/'filtered'/'{genome}'/'{sample}.bam.bai')
#    output:
#        str(MAPPING_FP/'filtered'/'{genome}'/'{sample}.depth')
#    shell:
#       """
#       samtools depth -aa {input.bam} > {output}
#       """
