#!/usr/bin/bash
set -e
set -x

## input
proj_fp=/mnt/chunyu6tb/sunbeam/Study9/sunbeam_output
data_fp=${proj_fp}/qc/cleaned

samples=`ls ${data_fp} | cut -d'_' -f1 | sort | uniq`

## output directories
midas_fp=${proj_fp}/classify/midas_results
mkdir -p ${midas_fp}

time_log=${proj_fp}/classify/log
mkdir -p ${time_log}


## tools
gnu_time_fp=/mnt/chunyu6tb/bin/time-1.9/time

## prepare the kmer database
midas_db=/mnt/chunyu6tb/dbs/IGGdb/v1.0.0/gtpro_repgenomes/midasdb_20190205

## run midas
for sample in ${samples}; do
  echo ${sample}
  mkdir -p ${midas_fp}/${sample}

  $gnu_time_fp -v \
   run_midas.py species ${midas_fp}/${sample} \
    -d $midas_db -t 1 -1 ${data_fp}/${sample}_1.fastq.gz \
    2>${time_log}/${sample}_species.log

  $gnu_time_fp -v \
   run_midas.py snps ${midas_fp}/${sample} \
    -d $midas_db -t 1 -1 ${data_fp}/${sample}_1.fastq.gz \
    --species_cov 3.0 \
    2>${time_log}/${sample}_snps.log
done
