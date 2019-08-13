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
