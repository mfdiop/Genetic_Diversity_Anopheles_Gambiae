
##################################
## This script performs sliding windows
## analysis of Nucleotide diversity, 
## Tajima's D for each population
##################################
library(tidyverse)
library(data.table)
library(readxl)

######################
## Check if zip file
######################
# # 1. Load library 'tools'
# library("tools")
# 
# # 2. Get extension for file
# if(file_ext(vcf_file) != "gz")
# {
#   system(sprintf("bgzip %s", vcf_file))
#   system(sprintf("tabix %s", paste0(vcf_file, ".gz")))
#   
#   vcf_file <- paste0(vcf_file, ".gz")
# }
# 
# system(sprintf("tabix %s", paste0(vcf_file)))

## LOAD VCF
VCF <- "../data/proceed/all_HPX15.vcf.gz"

## LOAD METADATA
metadata <- readxl::read_xlsx("../data/metadata/metadata.xlsx")

## Minor allele freq file
maf_file <- fread("../results/tables/per_species_minor_allel_freq.tsv", header = TRUE)

Species <- metadata %>% 
  select(species) %>% unique() %>% pull() %>% sort()

working_dir <- "../data/proceed/"
output_dir <- "../data/references/"
tajima <- "../results/tables/Tajima/"
nucleotide <- "../results/tables/Nucleotide/"
dir.create(tajima)
dir.create(nucleotide)

for(j in 1:length(Species))
{
  cat('===============================\n')
  cat("  Processing ", toupper(Species[j]), " ...\n")
  cat('===============================\n')
  # EXTRACT POPULATION SAMPLE IDs
  samplesTokeep <- paste0(working_dir, "samplesTokeep.txt")
  targetIDs <- metadata[which(metadata$species == Species[j]),]$sample_ID
  
  # SAVE SAMPLEIDs
  write.table(targetIDs, samplesTokeep, col.names = FALSE, 
              row.names = FALSE, quote = FALSE)
  
  # GET POSITIONS BASED ON MAF > 0.01
  snpsTokeep <- paste0(working_dir, "snpsTokeep.txt")
  
  maf_file %>% 
    filter(Species == Species[j] & maf > 0.01) %>% 
    select(chr, pos) %>% 
    # mutate(chr = paste0(chr, "L")) %>% 
    write.table(snpsTokeep, col.names = FALSE,
                row.names = FALSE, quote = FALSE, sep = '\t')
  
  #=========================================================
  #---- Extract SNPs and Isolates from filtered vcf file
  #========================================================
  
  OUT <- paste0(output_dir, Species[j])
  system(paste0("vcftools --gzvcf ", VCF,
                " --keep ", samplesTokeep, 
                " --positions ", snpsTokeep, 
                " --recode --recode-INFO-all --out ", OUT))
  
  file.remove(samplesTokeep, snpsTokeep)
  
  system(paste0("mv ", OUT, ".recode.vcf ", OUT, ".vcf"))
  system(paste0("bgzip ", OUT, ".vcf"))
  system(paste0("tabix ", OUT, ".vcf.gz"))
  
  cat('==========================================\n')
  cat("  Running Sliding Windows Analysis \n\ton ", toupper(Species[j]), " ...\n")
  cat('==========================================\n')
  vcfname <- paste0(OUT, ".vcf.gz")
  
  Populations <- metadata %>% filter(species == Species[j]) %>% 
    select(countries) %>% unique() %>% pull() %>% sort()
  
  pathToPiOutput <- paste0(nucleotide, Species[j]); dir.create(pathToPiOutput)
  pathToTjdOutput <- paste0(tajima, Species[j]); dir.create(pathToTjdOutput)
  
  for (p in 1:length(Populations)) {
    
    samplesTokeep <- paste0(working_dir, "samplesTokeep.txt")
    metadata %>% filter(species == Species[j] & countries == Populations[p]) %>% 
      select(sample_ID) %>% unique() %>% pull() %>% 
      write.table(samplesTokeep, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # cat('==========================================\n')
    # cat("  Running Pi on ", toupper(Populations[p]), " ...\n")
    # cat('==========================================\n')
    # Output <- file.path(pathToPiOutput, Populations[p])
    # # system(sprintf("vcftools --gzvcf %s --keep %s --window-pi 100 --window-pi-step 25 --out %s", 
    # #                vcfname, samplesTokeep, Output))
    # system(sprintf("vcftools --gzvcf %s --keep %s --window-pi 100 --out %s",
    #                vcfname, samplesTokeep, Output))
    
    cat('==========================================\n')
    cat("  Running TajimaD on ", toupper(Populations[p]), " ...\n")
    cat('==========================================\n')
    Output <- file.path(pathToTjdOutput, Populations[p])
    system(sprintf("vcftools --gzvcf %s --keep %s --TajimaD 50 --out %s",
                   vcfname, samplesTokeep, Output))
    
    file.remove(samplesTokeep)
  }
  
  system(sprintf("rm -r %s", vcfname))
  system(sprintf("rm -r %s", paste0(vcfname, ".tbi")))
}
