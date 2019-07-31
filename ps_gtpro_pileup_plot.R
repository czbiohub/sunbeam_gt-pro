## Chunyu Zhao 20190712; 20190718
#' @pileup_file: {sample}_{band}.tsv
#' @gtpro: {sample}.tsv
#' 

args <- commandArgs(trailingOnly = TRUE)

library(tools)
library(tidyverse)
library(readr)
library(stringr)
library(reshape2)
library(readxl)
library(magrittr)
library(forcats)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(grid)
library(data.table)
library(qiimer)


## Input arguments
sample <- args[1]
project_dir <- args[2] #"/mnt/chunyu20TB/hmp/sunbeam_output"

## Output argument
plot_file <- args[3] ## str(CLASSIFY_FP/'gtpro-r1'/'ps_gp_plot'/'{genome}'/'{sample}_gp.pdf'),
output_file <- args[4] ## str(CLASSIFY_FP/'gtpro-r1'/'ps_gp_plot'/'{genome}'/'{sample}_gp.tsv'),
stats_file <- args[5] ## str(CLASSIFY_FP/'gtpro-r1'/'ps_gp_plot'/'{genome}'/'{sample}_stats.csv')
corr_by_genome_file <- args[6] ## str(CLASSIFY_FP/'gtpro-r1'/'ps_gp_plot'/'{genome}'/'{sample}_bygenome_corr.csv')
jc_by_genome_file <- args[7] ## str(CLASSIFY_FP/'gtpro-r1'/'ps_gp_plot'/'{genome}'/'{sample}_bygenome_jc.csv')

use_aws = TRUE
if (!use_aws) {
  sample <- "SRS012849"
  project_dir <- "~/Desktop/20190714_plot/"
  species_info <- file.path(project_dir, "species_taxonomy_ext.tsv")
  gtpro_dir <- file.path(project_dir, "parse")
  pileup_dir <- file.path(project_dir, "banded")
} else {
  species_info <- file.path("/mnt/chunyu20TB/midas-iggdb/IGGdb/v2.0.0/species_taxonomy_ext.tsv")
  gtpro_dir <- file.path(project_dir, "classify/gtpro-r1/parse")
  pileup_dir <- file.path(project_dir, "mapping/filtered/repgenomes/banded") 
}

###################### Common parameters and functions
min_dp <- 2
min_gcb <- 10
min_allele_freq <- 0.01
NUM_THREADS = 48

specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))

parse_gtpro_file <- function(gtpro_file, min_dp=2, min_gcb=10){
  gtpro <- read_delim(gtpro_file, delim="\t", col_names = F) %>%
    mutate(depth = X7 + X8) %>% 
    dplyr::rename(genome_id = X1, ref_id = X3, ref_pos = X4) %>%
    mutate(genome_id = as.character(genome_id)) %>%
    mutate(ref_pos = ref_pos + 1)  %>% ## 0-based
    mutate(allele = paste(X5, X6, sep=":")) %>% 
    mutate(counts = paste(X7, X8, sep=":")) %>% 
    select(-one_of(c("X2","X5", "X6", "X7", "X8"))) %>%
    ## filter by min depth
    filter(depth >= min_dp) %>%
    ## filter by min genome covered bases
    group_by(genome_id) %>%
    filter(n() >= min_gcb) %>%
    ungroup()
}

gtpro_wide2long <- function(gtpro) {
  gtpro %<>%
    separate_rows(allele, counts, sep=":") %>% 
    mutate(counts = as.numeric(counts)) %>%
    mutate(freq = specify_decimal(counts/ depth, 3)) %>%
    mutate(gtdb_type = ifelse(row_number() %% 2 == 0, "gtminor", "gtmajor" ))
  
  gtpro %<>% dplyr::rename(gtdepth = depth, gtallele = allele, gtcounts = counts, gtfreq = freq)
  
  gtpro
}

###################### GTPRO genotyped sites
gtpro_file <- file.path(gtpro_dir, paste(sample, "_1.tsv", sep=""))
gtpro <- parse_gtpro_file(gtpro_file, min_dp=2, min_gcb=10)
sites_gtproed <- gtpro %>% mutate(ID = paste(genome_id, ref_id, ref_pos, sep="|")) %>% .$ID

f5 <- gtpro %>%
  group_by(genome_id) %>%
  summarise(site_counts = n()) %>%
  ungroup() %>%
  ggplot(aes(x = site_counts)) +
  geom_histogram() +
  theme_bw() + 
  ggtitle("Distributuion of site counts genotyped by Gtpro")

###################### PILEUP SNPs
## Read banded pileup files
pileup <- list()
n.bi.pileup <- 0
for (band in seq(0, NUM_THREADS-1)){
  param <- paste(".pileup.dp", min_dp, ".gcb", min_gcb, ".band", band, sep="")
  banded_pileup_file <- file.path(pileup_dir, paste(sample, param, ".tsv", sep=""))
  pileup.band <- fread(banded_pileup_file, skip=0, header=T, fill = T, blank.lines.skip=T) %>%
    select(-one_of(c("nz_allele", "nz_allele_count")))
  # genome_id,ref_id,ref_pos,ref_allele,depth,A,C,G,T,number_alleles,nz_allele,nz_allele_count
  n.bi.pileup = n.bi.pileup + pileup.band %>% filter(number_alleles == 2) %>% nrow()
  
  pileup.band %<>%
    mutate(genome_id = as.character(genome_id)) %>%
    mutate(ID = paste(genome_id, ref_id, ref_pos, sep="|")) %>%
    filter(ID %in% sites_gtproed) %>%
    select(-ID) %>%
    gather(allele, counts, A:T) %>%
    mutate(freq = as.numeric(specify_decimal(counts/depth, 3))) %>%
    filter(freq >= min_allele_freq)
  pileup[[as.character(band)]] <- pileup.band
}
pileup <- bind_rows(pileup)

## this is no longer making sense, since we didn't read in the whole pileup results
f4 <- pileup %>%
  group_by(genome_id) %>%
  summarise(site_counts = n()) %>%
  ungroup() %>%
  ggplot(aes(x = site_counts)) +
  geom_histogram() +
  theme_bw() + 
  ggtitle("(nah) Distributuion of site counts genotyped by Pileup")


###################### GP time
n.gtpro <- nrow(gtpro)
gtpro <- gtpro_wide2long(gtpro)

gp <- left_join(gtpro, pileup, by=c("genome_id", "ref_id", "ref_pos", "gtallele"="allele")) %>%
  mutate(isPileuped = ifelse(is.na(ref_allele), 0, 1)) %>% # agreed sites
  mutate(isBiallelicAgreed = ifelse(isPileuped==1 & number_alleles == 2, TRUE, FALSE)) # agreed bi-allelic SNPs defined by pileup

## these sites were genotyped by GTPRO but were never seen by Pileup
#gp %>% filter(isPileuped == 0) %>% View()
gp %<>% replace_na(list(depth = 0, number_alleles = 0, counts = 0, freq = 0))

gp %<>%
  group_by(genome_id, ref_id, ref_pos) %>%
  arrange(desc(isPileuped)) %>%
  filter(row_number() == 1) %>%
  ungroup()

## following should only be recovered fixed alleles
#gp %>% filter(isBiallelicAgreed == FALSE & isPileuped == 1)

###################### FIGURE time
gp.overlap <- gp %>% filter(isBiallelicAgreed)
## i need to see an example... and then re-run this part... fine...
corr <- specify_decimal(cor(gp.overlap$gtfreq, gp.overlap$freq),3)

corr_by_genome <- gp %>% 
  filter(isBiallelicAgreed) %>%
  group_by(genome_id) %>%
  summarise(correlation = specify_decimal(cor(gtfreq, freq, method="spearman"), 3), site_counts = n()) %>%
  ungroup()

## since gtpro.freq == 0 doesn't mean absence, then `jc` variable is completely wrong.
gp.pileuped <- gp %>% filter(isPileuped == 1) 
vec <- gp.pileuped %>% select(gtfreq, freq) %>% as.matrix()
jc <- specify_decimal(dist(t(vec), method = "binary"), 3)

cal_jaccard <- function(vec) {
  jc = specify_decimal(dist(t(as.matrix(vec)), method = "binary"), 3)
  return(jc)
}

## attention: 0 from gtpro doesn't mean absence, but 0 from pileup does. Need to confirm this with Katie.
## mutate(gtfreq = ifelse(gtfreq == 0, 0.001, as.numeric(gtfreq))) %>% 
## however I asked Jason and he said gtpro = 0 (for alt allele) indeed represent absence 
vec <- gp %>% 
  select(gtfreq, freq) %>% as.matrix()
jc.all <- specify_decimal(dist(t(vec), method = "binary"), 3)

jaccard_by_genome <- gp %>%
  #mutate(gtfreq = ifelse(gtfreq == 0, 0.001, as.numeric(gtfreq))) %>%  #<- this seems wrong 
  group_by(genome_id) %>%
  do(jaccard = cal_jaccard(.[, c("gtfreq", "freq")])) %>%
  mutate(jaccard = jaccard[[1]]) %>%
  ungroup()

n.agreed.bi <- gp %>% filter(isBiallelicAgreed) %>% nrow()
n.agreed <- gp %>% filter(isPileuped == 1) %>% nrow()
n.yaxis <- gp %>% filter(freq > 0)  %>% filter(freq < 1) %>% filter(gtfreq == 0) %>% nrow()
n.xaxis <- gp %>% filter(gtfreq > 0 ) %>% filter(gtfreq < 1) %>% filter(freq == 0) %>% nrow()
n.gtpro.fn <- gp %>% filter(isBiallelicAgreed) %>% filter(gtcounts == 0 & counts > 0) %>% nrow()

statsToshow <- data.frame(n.bi.pileup=n.bi.pileup, n.gtpro=n.gtpro, n.agreed.bi=n.agreed.bi, n.agreed=n.agreed, 
                          n.yaxis = n.yaxis, n.xaxis = n.xaxis, n.gtpro.fn=n.gtpro.fn, corr = corr, jaccard = jc,
                          jaccard.all=jc.all)

colorsList<- c("#064A77", "#2d83bc","#A6CEE3","#6a3d9a","#8e5ec1","#cfb1ef","#076302", "#33A02C", "#B2DF8A","#06ad97", "#63e2d1","#c6f2ec","#ffe100", "#f78002", "#f9a804","#fcc071",
               "#8c0101","#ff0000","#ffc4c4", "#af0166","#ff0094", "#fc9fd5","#dbdad0")
colourCount <- 8
mycolor <- colorRampPalette(colorsList)(colourCount+1)

## This among all the agreed bi-allelic gtpro SNPs
f1 <- gp %>%
  filter(isPileuped == 1) %>% #<- i think we should include these...
  select(genome_id, ref_id, ref_pos, gtfreq, freq, isBiallelicAgreed) %>%
  count(isBiallelicAgreed, gtfreq, freq) %>% 
  dplyr::rename(site_counts = n) %>% 
  ggplot(aes(x = gtfreq, y = freq, size = site_counts, color = isBiallelicAgreed)) +
  geom_point() + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "gtpro allele frequency", y = "correponding alleles's frequency by pileup") +
  theme(axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~isBiallelicAgreed) +
  ggtitle(paste("Bi-allelic sites agreed by GTPRO and Pileup for", sample))

f2 <- tableGrob(statsToshow, rows = NULL)

###################### ABUNDANT SPECIES
iggdb_species <- read_delim(species_info, delim="\t", col_names = c("orig_genome_id", "genome_id", "classification", "ext")) %>%
  mutate(classification = ifelse(is.na(ext), classification, paste(classification, ext, sep = " ")))
iggdb_species <- cbind(iggdb_species, qiimer::split_assignments(iggdb_species$classification, split = ";")) %>%
  mutate(species_name = ifelse(Species == "s__", paste(Genus, Species, sep=""), as.character(Species))) %>%
  select(genome_id, species_name) %>%
  mutate(genome_id = as.character(genome_id))

top.genomes <- gp %>% 
  group_by(genome_id) %>% 
  summarise(n_sites = n(), totalCount = sum(gtcounts), avgCov = totalCount/n_sites) %>%
  arrange(desc(avgCov)) %>% 
  head(n = 8) %>% 
  .$genome_id

f3 <- gp %>%
  select(genome_id, ref_id, ref_pos, gtfreq, freq, isPileuped) %>%
  count(genome_id, isPileuped, gtfreq, freq) %>%
  dplyr::rename(site_counts = n) %>% 
  left_join(iggdb_species, by=c("genome_id")) %>%
  mutate(species_name = ifelse(genome_id %in% top.genomes, paste(genome_id, species_name, sep="\n"), "Other")) %>%
  ggplot(aes(x = gtfreq, y = freq, color = species_name, size=site_counts)) + 
  geom_point() +
  scale_color_manual(values = mycolor, na.value="grey", name="species_name") +
  #geom_abline(intercept = 0, slope = 1, linetype="dashed", size=1, color = "red") + 
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.5)) + ## min maf is 0.1
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.5)) +
  theme_bw() + 
  guides(colour = FALSE) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
  labs(x = "gtpro allele frequency", y = "correponding pileup allele's frequency") +
  theme(axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste("common allele frequency for", sample)) +
  facet_grid(isPileuped ~ species_name ) 

###################### PLOT
set.seed(123)
grid.newpage()
pdf(plot_file,width=18,height=12,paper='special') 
grid.arrange(f1, f4, f5, f3, f2, ncol = 3,
             heights=c(1.4,3,0.8),widths = c(1.4, 1,1),
             layout_matrix = rbind(c(1,2, 3), c(4,4,4), c(5, 5, 5)))
dev.off()
fwrite(gp, output_file, row.names=F, quote=F, sep=",")
fwrite(statsToshow, stats_file, row.names=F, quote=F, sep=",")
fwrite(corr_by_genome, corr_by_genome_file, row.names=F, quote=F, sep=",")
fwrite(jaccard_by_genome, jc_by_genome_file, row.names=F, quote=F, sep=",")
