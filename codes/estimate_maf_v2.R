
# estimate_maf.R
#
# Author: Fatoumata SECK HPC
# Date: 2022-01-19
#
# Purpose:
# Estimate Minor Allele Frequency for each species.
#
# ------------------------------------------------------------------

arguments <- commandArgs(trailingOnly = TRUE)

if(length(arguments) == 0){
    cat("=============================================\n")
    cat("SCRIPT SHOULD BE RUN LIKE: \n")
    cat("Rscript estimate_maf.R file.vcf.gz")
    cat("\n")
    stop("NOT ENOUGH ARGUMENTS SUPPLIED .....")
    cat("=============================================\n")
}

## Uncomment these lines if you want to install the packages
# packages <- c("tidyverse", "data.table")
# install.packages(setdiff(packages, rownames(installed.packages())))

# load  packages
library(tidyverse)
library(data.table)

working_dir <- "../Data/process/"
output_dir <- "../results/tables/Allel_Freq/"
system(sprintf("mkdir -p %s", output_dir))
system(sprintf("mkdir -p %s", working_dir))

#==============
#   FUNCTION
#=============
computeMAF <- function(VCF_OUT, targetSamples, OUT, genename)
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
  maf_file <- paste0("../results/tables/", genename, "_PerSpecies_MAF.tsv")
  if (!file.exists(maf_file)) {
    file.create(maf_file)
    write.table(var_freq, maf_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
  }
  else 
    write.table(var_freq, maf_file, col.names = FALSE, append = TRUE,
                row.names = FALSE, quote = FALSE, sep = '\t')
  
  file.remove(paste0(OUT, ".frq"))
}

## LOAD METADATA
metadata <- readxl::read_xlsx("../Data/metadata/global_Metadata.xlsx")

## LOAD VCF
VCF <- arguments[1]

# ## FILTER VCF BASED ON 
# # set filters
# MISS = 0.9
# QUAL = 30 
# MIN_DEPTH = 10
# MAX_DEPTH = 50
# VCF_OUT <- paste0(working_dir, gsub(".vcf.gz", "", basename(VCF)))
# 
# system(sprintf("vcftools --gzvcf %s --remove-indels --max-missing %s --minQ %s --min-meanDP %s --max-meanDP %s --minDP %s --maxDP %s --recode --recode-INFO-all --out %s",
#                VCF, MISS, QUAL, MIN_DEPTH, MAX_DEPTH, MIN_DEPTH, MAX_DEPTH, VCF_OUT))
# 
# ## Rename vcf file
# system(paste0("mv ", VCF_OUT, ".recode.vcf ", VCF_OUT, ".vcf"))
# 
# ## Compress vcf
# system(sprintf("bgzip %s", paste0(VCF_OUT, ".vcf")))
# 
# VCF_OUT <- paste0(VCF_OUT, ".vcf.gz")

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
  targetSamples <- paste0(working_dir, Species[j], ".txt")
  targetIDs <- metadata[which(metadata$species == Species[j]),]$sample_ID
  
  ## SAVE SAMPLEIDs
  write.table(targetIDs, targetSamples, col.names = FALSE, 
              row.names = FALSE, quote = FALSE)
  
  OUT <- paste0(output_dir, Species[j])
  cat('===============================\n')
  cat("Computing MAF on ", Species[j], " data ...\n")
  cat('===============================\n')
  
  genename <- gsub(".vcf.gz", "", basename(VCF))
  ## COMPUTE MINOR ALLELE FREQUENCY
  computeMAF(VCF_OUT, targetSamples, OUT, genename)
}

system(paste0("rm -rf ", output_dir)) # will delete directory called 'output_dir'

