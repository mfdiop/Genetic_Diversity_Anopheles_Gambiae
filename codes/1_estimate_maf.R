#===================================================
# 1_estimate_maf.R
#
# Author: Mouhamadou F. DIOP
# Date: 2022-07-10
#
# Purpose:
# Estimate Minor Allele Frequency for each species.
#
# https://speciationgenomics.github.io/filtering_vcfs/
#===================================================

rm(list = ls())

VCF <- commandArgs(trailingOnly = TRUE)

if(length(VCF) == 0){
    print("===================================================", quote = FALSE)
    print("        Estimate Minor Allele Frequency            ", quote = FALSE)
    print("===================================================", quote = FALSE)
    print("", quote = FALSE)
    print(" Usage: Rscript 1_estimate_maf.R file.vcf.gz", quote = FALSE)
    stop(" Not enough arguments supplied ...")
    
}

## Uncomment these lines if you want to install the packages
# packages <- c("tidyverse", "data.table")
# install.packages(setdiff(packages, rownames(installed.packages())))

# Load  packages
library(tidyverse)
library(data.table)

#==============
#   FUNCTIONS
#=============
# MINOR ALLELE FREQUENCY ACROSS SPECIES
computeMAF <- function(FILTERED_VCF, OUT)
{
  system(sprintf("vcftools --gzvcf %s --freq2 --out %s ", FILTERED_VCF, OUT))

  ## LOAD ALLELE FREQUENCY FILES
  var_freq <- fread(paste0(OUT, ".frq"))
  names(var_freq) = c("chr", "pos", "nalleles", "nchr", "a1", "a2", "a3", "a4")   
  
  var_freq <- var_freq %>% 
    group_by(id = row_number()) %>%   # for each row
    nest(data = c(a2, a3, a4)) %>%    # nest selected columns
    mutate(SUM = map_dbl(data, sum))  # calculate the sum of those columns
  
  var_freq <- var_freq %>% 
    select(-c(id, data)) %>% 
    rename(a2 = SUM)
  
  # FIND MINOR ALLELE FREQUENCY
  var_freq$maf <- var_freq %>% 
    select(a1, a2) %>% 
    apply(1, function(z) min(z))
  
  ## SAVE MAF FILE
  maf_file <- file.path(OUT, "MAF.tsv")
  if(!file.exists(maf_file)) file.create(maf_file)

  write.table(var_freq, maf_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
  file.remove(paste0(OUT, ".frq"))
}

#================================================
#================================================
# MINOR ALLELE FREQUENCY PER SPECIES
computeMAF_perSpecies <- function(VCF_OUT, targetSamples, OUT, MAF_OUTPUT)
{
  system(sprintf("vcftools --gzvcf %s --keep %s --freq2 --out %s ", VCF_OUT, targetSamples, OUT))
  file.remove(targetSamples)
  
  ## LOAD ALLELE FREQUENCY FILES
  var_freq <- fread(paste0(OUT, ".frq"))
  names(var_freq) = c("chr", "pos", "nalleles", "nchr", "a1", "a2", "a3", "a4")   
  
  var_freq <- var_freq %>% 
    group_by(id = row_number()) %>%   # for each row
    nest(data = c(a2, a3, a4)) %>%    # nest selected columns
    mutate(SUM = map_dbl(data, sum))  # calculate the sum of those columns
  
  var_freq <- var_freq %>% 
    select(-c(id, data)) %>% 
    rename(a2 = SUM) %>% 
    mutate(Species = basename(OUT))
  
  # FIND MINOR ALLELE FREQUENCY
  var_freq$maf <- var_freq %>% 
    select(a1, a2) %>% 
    apply(1, function(z) min(z))
  
  ## SAVE MAF FILE
  maf_file <- paste0(MAF_OUTPUT, "/PerSpecies_MAF.tsv")
  if (!file.exists(maf_file)) {
    file.create(maf_file)
    write.table(var_freq, maf_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
  }
  else 
    write.table(var_freq, maf_file, col.names = FALSE, append = TRUE,
                row.names = FALSE, quote = FALSE, sep = '\t')
  
  file.remove(paste0(OUT, ".frq"))
}

#=======================================================================
## LOAD METADATA
metadata <- readxl::read_xlsx("data/metadata/global_Metadata.xlsx")

working_dir <- "data/process/"
output_dir <- "results/tables/Allel_Freq/"
system(sprintf("mkdir -p %s", output_dir))
system(sprintf("mkdir -p %s", working_dir))


## FILTER VCF BASED ON 
# set filters
MISS = 0.9
QUAL = 30 # Phred Site Quality: how much confidence we have in our variant calls. 
# Phred score of 30 represents a 1 in 1000 chance that our SNP call is erroneous.
# Variant mean depth
MIN_DEPTH = 10
MAX_DEPTH = 50
FILTERED_VCF <- paste0(working_dir, gsub("\\.[aA-zZ]*", "", basename(VCF)))

system(sprintf("vcftools --gzvcf %s --remove-indels --max-missing %s --minQ %s --min-meanDP %s --max-meanDP %s --minDP %s --maxDP %s --recode --recode-INFO-all --out %s",
               VCF, MISS, QUAL, MIN_DEPTH, MAX_DEPTH, MIN_DEPTH, MAX_DEPTH, FILTERED_VCF))

## Rename vcf file
system(paste0("mv ", FILTERED_VCF, ".recode.vcf ", FILTERED_VCF, ".vcf"))

## Compress vcf
system(sprintf("bgzip -f %s", paste0(FILTERED_VCF, ".vcf")))

FILTERED_VCF <- paste0(FILTERED_VCF, ".vcf.gz")

#====================================
# MINOR ALLELE FREQUENCY ACROSS SPECIES
#====================================
genename <- gsub(".vcf.gz", "", basename(VCF))
MAF_OUTPUT <- paste0("results/tables/", genename)

if(!dir.exists(MAF_OUTPUT)) dir.create(MAF_OUTPUT)

computeMAF(FILTERED_VCF, MAF_OUTPUT)

#====================================
# MINOR ALLELE FREQUENCY PER SPECIES
#====================================
Species <- metadata %>% 
  select(species) %>% unique() %>% pull() %>% sort()

for(j in 1:length(Species))
{
  cat('===============================\n')
  cat("Processing ", Species[j], " ...\n")
  cat('===============================\n')
  
  ## EXTRACT POPULATION SAMPLE IDs
  targetSamples <- paste0(output_dir, Species[j], ".txt")
  targetIDs <- metadata[which(metadata$species == Species[j]),]$sample_ID
  
  ## SAVE SAMPLEIDs
  write.table(targetIDs, targetSamples, col.names = FALSE, 
              row.names = FALSE, quote = FALSE)
  
  OUT <- paste0(output_dir, Species[j])
  cat('===============================\n')
  cat("Computing MAF on ", Species[j], " data ...\n")
  cat('===============================\n')
  
  ## COMPUTE MINOR ALLELE FREQUENCY
  computeMAF_perSpecies(FILTERED_VCF, targetSamples, OUT, MAF_OUTPUT)
  # computeMAF_perSpecies(VCF, targetSamples, OUT, genename)
}

system(paste0("rm -rf ", output_dir)) # will delete directory called 'output_dir'

