
rm(list = ls())
# neutrality_analysis.R
#
# Author: Mouhamadou F. DIOP
# Date: 2021-11-20
#
# Purpose:
# Produce Minor Allele Frequency for each species.
#
# ------------------------------------------------------------------

## Load libraries
library(tidyverse)
library(data.table)
library(PopGenome)

## LOAD VCF
VCF <- "data/references/MISO_maf_filtered.vcf.gz"

## LOAD METADATA
metadata <- readxl::read_xlsx("data/metadata/global_Metadata.xlsx")

genename <- gsub("_maf_filtered.vcf.gz", "", basename(VCF))

Species <- metadata %>% 
  select(species) %>% unique() %>% pull() %>% sort()

for(j in 1:length(Species))
{
  cat('===============================\n')
  cat("  Processing ", Species[j], " ...\n")
  cat('===============================\n')
  # EXTRACT POPULATION SAMPLE IDs
  samplesTokeep <- "data/references/samplesTokeep.txt"
  targetIDs <- metadata[which(metadata$species == Species[j]),]$sample_ID
  number_samples <- length(targetIDs)
  
  # SAVE SAMPLEIDs
  write.table(targetIDs, samplesTokeep, col.names = FALSE, 
              row.names = FALSE, quote = FALSE)
  
  #=============================================
  #---- Extract Isolates from filtered vcf file
  #=============================================
  
  OUT <- paste0("data/references/", Species[j])
  system(paste0("vcftools --gzvcf ", VCF,
                " --keep ", samplesTokeep, 
                " --recode --recode-INFO-all --out ", OUT))
  
  file.remove(samplesTokeep)
  
  system(paste0("mv ", OUT, ".recode.vcf ", OUT, ".vcf"))
  system(paste0("bgzip ", OUT, ".vcf"))
  system(paste0("tabix ", OUT, ".vcf.gz"))
  
  ####################################
  ### Get readVCF parameters from file
  ####################################
  cat('==========================================\n')
  cat("  Running Neutrality Stats on ", Species[j], " ...\n")
  cat('==========================================\n')
  vcfname <- paste0(OUT, ".vcf.gz")
  skipNum <- as.numeric(system(paste0("zgrep  -P \"##\" ", vcfname, " | wc -l"), TRUE))
  vcf <- read.table(vcfname, skip = skipNum, header = TRUE, comment.char = "", 
                    stringsAsFactors = FALSE, check.names = FALSE)
	
  names(vcf)[1] <- "CHROM"
  
  chroms <- vcf %>% 
	select(CHROM) %>% 
	unique() %>% 
	pull()
  frompos <- vcf$POS[1]
  topos <- vcf$POS[nrow(vcf)]
  samplenames <- names(vcf)[-c(1:9)]
  
  ##########################
  ## Set list of populations
  ##########################
  Population <- metadata %>% 
    filter(species == Species[j]) %>% 
    select(countries, sample_ID)
    
  population <- data.frame( value = Population, stringsAsFactors = T) %>%
    group_by(value.countries) %>%
    summarize(list_result = list(value.sample_ID)) %>% .$list_result
  
  #################
  ### Read vcf file
  #################
  Genome <- readVCF(vcfname, numcols = 10000, tid = chroms, from = frompos, 
                    to = topos, samplenames = samplenames, include.unknown = TRUE)
  
  ## Set the populations
  Genome <- set.populations(Genome, population, diploid = TRUE)
  
  ### Compute Neutrality Tests
  GENOME.class <- neutrality.stats(Genome, FAST = TRUE)
  
  seg_sites <- GENOME.class@n.segregating.sites %>% as.vector()
  Tajima <- round(GENOME.class@Tajima.D, 3) %>% as.vector()
  FuLi.F <- round(GENOME.class@Fu.Li.F, 3) %>% as.vector()
  FuLi.D <- round(GENOME.class@Fu.Li.D, 3) %>% as.vector()
  Theta <- round(GENOME.class@theta_Watterson, 3) %>% as.vector()
  
  ###########################################
  #---- Nucleotide and haplotype diversity
  #---- Based on Nielson Theory
  ##########################################
  Genome.diversity <- diversity.stats(Genome)
  
  Nucleotide_Diversity <- round(Genome.diversity@nuc.diversity.within, 3) %>% as.vector()
  Haplotype_Diversity <- round(Genome.diversity@hap.diversity.within, 3) %>% as.vector()
  
  
  number_samples <- table(Population$countries) %>% as.vector()
  Populations <- unique(Population$countries) %>% sort()
  Neutrality.stats <- tibble(Species = Species[j], Populations = Populations,
                             N = number_samples, NSS = seg_sites, Pi = Nucleotide_Diversity,
                             Hd = Haplotype_Diversity, Theta = Theta,
                             Tajima = Tajima, FuLi_F = FuLi.F, FuLi_D = FuLi.D)
  
	per_species_stats <- paste0("results/tables/", genename, "_PerSpecies_population_neutrality_stats.tsv")
  ## Save file
  if(!file.exists(per_species_stats)){
    file.create(per_species_stats)
    write.table(Neutrality.stats, per_species_stats, col.names = TRUE,
                row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  else 
    write.table(Neutrality.stats, per_species_stats, append = TRUE, 
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Remove VCF file
  file.remove(vcfname, paste0(vcfname, ".tbi"))
  
}
