## Chunyu Zhao 20190326
args <- commandArgs(trailingOnly = TRUE)

library(tools)
library(tidyverse)
library(readr)
library(stringr)
library(reshape2)
library(readxl)
library(magrittr)
library(forcats)

## Input arguments
depth_file <- args[1]
genome_length_file <- args[2]
## Output argument
genome_stats_file <- args[3]


specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))

###################### GENOME_STATS: read in genome length and pileup depth
pileup.depth <- read_delim(depth_file, delim= "\t")

genome_stats <- pileup.depth %>%
  group_by(ref_id) %>%
  summarise(contig_total_depth = sum(depth), contig_covered_bases = n()) %>%
  ungroup()

genome_length <- read_delim(genome_length_file, delim="\t", col_names = c("contig_id", "contig_length")) %>%
  group_by(contig_id) %>%
  summarise(contig_length = sum(contig_length)) %>%
  ungroup()

genome_stats %<>% left_join(genome_length, by=c("ref_id" = "contig_id")) %>%
  dplyr::rename(contig_id = ref_id) %>%
  separate(contig_id, into=c("genome_id"), sep="\\|", extra="drop", remove = FALSE) %>%
  mutate(contig_id = gsub("\\|", "_", contig_id)) %>%
  mutate(frac_covered = specify_decimal(contig_covered_bases / contig_length, 3)) %>%
  mutate(mean_coverage = specify_decimal(contig_total_depth / contig_covered_bases, 2))
  
genome_stats %>% write.table(genome_stats_file, sep="\t", quote=F, row.names = F)
###################### GENOME_STATS: read in genome length and pileup depth
