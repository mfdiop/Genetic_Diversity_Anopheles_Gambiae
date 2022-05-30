
rm(list = ls())
# generate_maf_plot.R
#
# Author: Mouhamadou F. DIOP
# Date: 2021-11-20
#
# Purpose:
# Pi and Tajima's D.
#
# ------------------------------------------------------------------

## Load libraries
library(tidyverse)
library(data.table)
library(PopGenome)

## LOAD VCF
VCF <- "data/references/HPX14_maf_filtered.vcf.gz"
genename <- gsub("_maf_filtered.vcf.gz", "", basename(VCF))

OUTPUT <- "results/tables/sliding_windows/"
if(!dir.exists(OUTPUT)) dir.create(OUTPUT)

## LOAD METADATA
metadata <- readxl::read_xlsx("data/metadata/global_Metadata.xlsx")

Species <- metadata %>% 
  select(species) %>% unique() %>% pull() %>% sort()

for(j in 2:length(Species))
{
  cat('===============================\n')
  cat("  Processing ", Species[j], " ...\n")
  cat('===============================\n')
  # EXTRACT POPULATION SAMPLE IDs
  samplesTokeep <- "data/references/samplesTokeep.txt"
  target <- metadata[which(metadata$species == Species[j]),]
  
  ## Extract Populations and their sample IDs
  Population <- target %>% 
      select(countries, sample_ID) %>% 
      group_by(countries) %>% # Group by Population
      filter(n() >= 10) %>% # Remove populations with low number of isolates
      ungroup()
  
  ## Create a list for each Pop/sampleIDs
  population <- data.frame( value = Population, stringsAsFactors = T) %>%
      group_by(value.countries) %>%
      summarize(list_result = list(value.sample_ID)) %>% .$list_result
  
  targetIDs <- Population %>% pull(sample_ID)
  
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
  filtered_vcf <- paste0(OUT, ".vcf.gz")
  skipNum <- as.numeric(system(paste0("zgrep  -P \"##\" ", filtered_vcf, " | wc -l"), TRUE))
  vcf <- read.table(filtered_vcf, skip = skipNum, header = TRUE, comment.char = "", 
                    stringsAsFactors = FALSE, check.names = FALSE)
  
  names(vcf)[1] <- "CHROM"
  
  chroms <- vcf %>% pull(CHROM) %>% unique(); 
  frompos <- vcf$POS[1]; 
  topos <- vcf$POS[nrow(vcf)];
  samplenames <- names(vcf)[-c(1:9)]
  
  ##########################
  ## Set list of populations
  ##########################
  # ## Extract Populations and their sample IDs
  # Population <- metadata %>% 
  #   filter(species == Species[j]) %>% 
  #   select(countries, sample_ID) %>% 
  #     group_by(countries) %>% 
  #     filter(n() >= 10)
  # 
  # ## Create a list for each Pop/sampleIDs
  # population <- data.frame( value = Population, stringsAsFactors = T) %>%
  #   group_by(value.countries) %>%
  #   summarize(list_result = list(value.sample_ID)) %>% .$list_result
  
  #################
  ### Read vcf file
  #################
  Genome <- readVCF(filtered_vcf, numcols = 10000, tid = chroms, from = frompos, 
                    to = topos, samplenames = samplenames, include.unknown = TRUE)
  
  Genome@n.biallelic.sites
  ## Set the populations
  Genome <- set.populations(Genome, population, diploid = TRUE)
  
  ## Sliding windows
  # split the data in 10kb consecutive windows
  slide <- sliding.window.transform(Genome,100,100, type = 2)
  
  # total number of windows
  length(slide@region.names)
  # Statistics
  slide <- diversity.stats(slide)
  nucdiv <- slide@nuc.diversity.within
  
  # the values have to be normalized by the number of nucleotides in each window
  nucdiv <- nucdiv/length(slide@region.names)
  head(nucdiv)
  
  nuc.div <- as_tibble(nucdiv); 
  names(nuc.div) <- Population %>% 
    select(countries) %>% 
    unique() %>% 
    pull() %>% 
    sort()
  
  write.table(nuc.div, paste0(OUTPUT, genename, "_", Species[j], "_pi.tsv"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
  
  slide <- neutrality.stats(slide, FAST = TRUE)
  tajima <- round(slide@Tajima.D, 3)
  
  # the values have to be normalized by the number of nucleotides in each window
  # tajima <- tajima/length(slide@region.names)

  tajima <- as_tibble(tajima)
  names(tajima) <- Population %>% 
    select(countries) %>% 
    unique() %>% 
    pull() %>% 
    sort()
  
  write.table(tajima, paste0(OUTPUT, genename, "_", Species[j], "_tajima.tsv"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
  
  # Remove VCF file
  file.remove(paste0(OUT, ".vcf.gz"), paste0(OUT, ".vcf.gz.tbi"))
  
}
