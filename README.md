# Validation of gt-pro for microbial genotyping using shotgun metagenomics sequencing data

## CONDA environment

We need to have the proper channel order in the *.condarc* file, when downloading CONDA for the first time.

```bash
conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge
```

## Downalod SRA sequencing data

Use [grabseqs](https://github.com/louiejtaylor/grabseqs) to download sequencing data from NCBI SRA. 

### HMP dataset

Download the metadata file for HMP dataset.

There are 1576 paired-ends reads samples from 337 SubjectIDs. After filter by Stephen's metadata and Nandita’s SuppleN Table S1, there are **330** samples left.

```bash
grabseqs sra -t 32 -m -o hmp/ -r 3 SRP002163 -l
```

To better handle the potential incomplete downloads, we used *sunbeam init* command for secure download.

```bash
sunbeam init hmp --data_acc SRP002163
```

In this dataset, we have both the **technical replicates** (same sample, multipe sequencing runs) and **biological replicates** (samples from the same individual sequenced at different time points).

### Madagascar dataset

```bash
grabseqs sra -t 32 -m -o Madagascar/ -r 3 PRJNA485056 -l
sunbeam init Madagascar --data_acc PRJNA485056
```

## Install tools

### [Sunbeam](https://github.com/sunbeam-labs/sunbeam) pipeline

### [gt-pro2.0](https://github.com/zjshi/gt-pro2.0)

### [MIDAS](https://github.com/snayfach/MIDAS)

```bash
conda env update --name=midas --file ymlfiles/midas.yml

git clone https://github.com/snayfach/MIDAS
cd MIDAS
python setup.py install

cd ../dbs
wget http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.2.tar.gz
tar -zxvf midas_db_v1.2.tar.gz
dbs_fp=`pwd`
```

**MIDAS unit testing**

```bash
midasdb_fp=`pwd`
export PYTHONPATH=$PYTHONPATH:"${midasdb_fp}"
export PATH=$PATH:"${midasdb_fp}/scripts"
export MIDAS_DB="${dbs_fp}/midas_db_v1.2"

cd ../MIDAS/test
./test_midas.py -vf
```

**run MIDAS**

```bash
conda activate midas
ls $HOME/hmp/sunbeam_output/qc/cleaned/ | cut -d'_' -f1 | sort | uniq | \
      xargs -I{} -n1 -P 48 bash -c "bash run_midas_sample.sh {}"
```

## Bioinformatics processing for shotgun metagenomics samples

### Quality control

- qc/qc.rules: Trimmomatic + Komplexity

```bash
snakemake --configfile HMP/sunbeam_config.yml all_qc --cores 48
```

### Alignment-based SNP calling

- mapping/mapping.rules

```bash
snakemake --configfile HMP/sunbeam_config.yml _all_pileup --cores 48
```

### Run gt-pro2.0

- Split the samples into batches; and then we ran gt-pro2.p on R1 reads only (rules/gtpro.rules); rename and parse gt-pro output.

```
cd $HOME/hmp
split -l 48 <(cat samples.csv) split_R1/samples_

snakemake sunbeam_config.yml _all_raw_gtpro --cores 48
snakemake sunbeam_config.yml _all_rename_gtpro --cores 48
snakemake sunbeam_config.yml _all_parse_gtpro --cores 48
```

## Validation experiments

### Per sample validation using [xsnp](https://github.com/czbiohub/xsnp)

- correlation
- concordance
- output: *genome/contigs_stats* and *bi-plot* for each sample

```bash
cd $HOME/hmp/sunbeam_output/mapping/filtered/repgenomes/
$HOME/xsnp/scripts/go.py hmp_all.txt

snakemake sunbeam_config.yml _all_ps_gp_plot --cores 48
```

### Across samples core SNPs for pileup

```bash
for sample in `cat hmp_all.txt`; do for x in {0..47}; do 
  echo banded/${sample%.pileup}.pileup.dp2.gcb10.band${x}.tsv >> banded/band${x}.txt; done; 
done


for i in {0..47}; do echo band${i}.txt; $HOME/xsnp/scripts/go2.py "banded/band${i}".txt; done
```

### Across samples core SNPs using [xsnp-gtpro](https://github.com/zhaoc1/xsnp-gtpro)

```bash
$HOME/xsnp-gtpro/scripts/go.py hmp_all.txt

for sample in `cat hmp_all.txt`; do for x in {0..47}; do 
  echo banded/${sample%_1.tsv}.gtpro.dp2.gcb10.band${x}.tsv >> banded/band${x}.txt; done; 
done

for i in {0..47}; do echo band${i}.txt; $HOME/xsnp-gtpro/scripts/go2.py "banded/band${i}".txt; done
```

# Gtpro can detect low abundant species in the metagenomics samples

Dataset: The 29 C. difficile genome sequences and their episomes were deposited at NCBI
under PRJNA524299; metagenomic DNA sequence data is deposited under
PRJNA562600. Metabolomic data and computer code used in this study are available
at https://github.com/zhaoc1/nanoflow, https://github.com/zhaoc1/coreSNPs, and
https://github.com/reny1/cdiff.

```bash
grabseqs sra -t 32 -m -o cdiff/ -r 3 PRJNA562600 -l
sunbeam init cdiff --data_acc PRJNA562600
```


