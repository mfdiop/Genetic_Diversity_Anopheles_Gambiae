
rm(list = ls())

# fasta_from_vcf.R
#
# Author: Mouhamadou F. DIOP
# Date: 2021-11-20
#
# Purpose:
# Estimate .
#
# ------------------------------------------------------------------
library(tidyverse)

VCF <- "../data/references/HPX15_maf_filtered.vcf.gz"

## LOAD METADATA
metadata <- readxl::read_xlsx("../data/metadata/metadata.xlsx")

## Minor allele freq file
genename <- gsub("_maf_filtered.vcf.gz", "", basename(VCF))
maf_file <- data.table::fread(paste0("../results/tables/", genename, "_MAF.tsv"), 
                              header = TRUE)

Genotypes = "../data/references/Genotypes.txt"
AllelicDeph = "../data/references/ReadDepth.txt"
GTexpression = '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
ADexpression = '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n'
samples = "../data/references/isolates.txt"

system(sprintf("bcftools query -f'%s' %s > %s", GTexpression, VCF, Genotypes))  
system(sprintf("bcftools query -f'%s' %s > %s", ADexpression, VCF, AllelicDeph))
system(sprintf("bcftools query -l %s > %s", VCF, samples))

# ---------- Phased the final vcf file: FinalSelectionData.vcf.gz

genotypeData = fread(Genotypes, header = F)
allelicDepthData = fread(AllelicDeph, header = F)
isolates <- scan(samples, what = 'character')

first4Columns = subset(genotypeData, select=c(1:4))
genotypeData = as.data.frame(subset(genotypeData, select=-c(1:4)))
allelicDepthData = as.data.frame(subset(allelicDepthData, select=-c(1:4)))

write.table(first4Columns, "../data/references/first4Columns.txt", 
            quote = FALSE, row.names = FALSE, col.names = F)  
FileName ='../results/PhasedData.txt'

# genotypeData <- as.matrix(genotypeData)
# genotypeData[genotypeData=="0/0"]=0
# genotypeData[genotypeData=="1/1"]=1
# genotypeData[genotypeData=="2/2"]=2
# genotypeData[genotypeData=="3/3"]=3
# genotypeData <- as.data.frame(genotypeData)

phasedData = matrix(9, nrow=dim(genotypeData)[1], ncol=dim(genotypeData)[2])
het <- c("0/1", "0/2", "0/3", "1/2", "1/3", "2/3")

for(j in 1:nrow(genotypeData))
{
  k=1
  while(k<=ncol(genotypeData))
  {
    if(genotypeData[j,k]=='0/0')
      phasedData[j,k]=0   
    else if(genotypeData[j,k]=='./.')
      phasedData[j,k]='.'
    else if(genotypeData[j,k]=='1/1')
      phasedData[j,k]=1
    else if(genotypeData[j,k]=='2/2')
      phasedData[j,k]=2
    else if(genotypeData[j,k]=='3/3')
      phasedData[j,k]=3
    else if(genotypeData[j,k] %in% het)
    {
      target = as.integer(unlist(strsplit(allelicDepthData[j,k], ',')))
      index <- which(target == max(target))
      if(index == 1) phasedData[j,k] = 0
      else if(index == 2) phasedData[j,k] = 1
      else if(index == 3) phasedData[j,k] = 2
      else phasedData[j,k] = 3
    }
    k = k+1
  }
}

phasedData = as.data.frame(cbind(first4Columns, phasedData))
names(phasedData) <- c("CHROM", "POS", "REF", "ALT", isolates)

phasedData <- as.matrix(phasedData)
phasedData[phasedData =='.'] <- NA

phasedData <- as.tibble(phasedData)

phasedData <- phasedData[ , colSums(is.na(phasedData)) == 0]

write.table(phasedData, FileName, quote = FALSE, row.names = FALSE, col.names = TRUE)

#===============================
# CONVERT GENOTYPE TO HAPLOTYPE
#==============================
ref <- phasedData %>% 
  select(REF) %>% 
  pull()

alt <- phasedData %>% 
  select(ALT) %>% 
  pull()

phasedData <- phasedData %>% 
  select(-c(1:4))
phasedData <- type_convert(phasedData)

subset.data <- phasedData[, sample(ncol(phasedData), 150)]

haplotype <- matrix(NA, nrow=dim(phasedData)[1], ncol=dim(phasedData)[2])

for (i in 1:nrow(phasedData)) {
  cat("Computing SNPs ", i, "\n")
  alts <- unlist(strsplit(alt[i], ','))
  for (j in 1:ncol(phasedData)) {
    if(phasedData[i,j] == 0) haplotype[i,j] <- ref[i]
    else if(phasedData[i,j] == 1) haplotype[i,j] <- alts[1]
    else if(phasedData[i,j] == 2) haplotype[i,j] <- alts[2]
    else if(phasedData[i,j] == 3) haplotype[i,j] <- alts[3]
  }
}

haplotype <- as.tibble(haplotype); names(haplotype) <- names(phasedData)
write.table(haplotype, "../results/hpx15.haplotype.txt", sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = TRUE)

#=============================
# CONVERT HAPLOTYPE TO FASTA
#============================
outputFasta <- "../results/hpx15.fasta"
fasta_from_haplotype <- function(haplotype, sample_names, outputFasta)
{
  vecteur<- vector(mode="character", length(sample_names))
  for (i in 1: length(sample_names))
  {
    vecteur[i] = paste(">",sample_names[i])
  }
  
  transposer = t(haplotype)
  
  fichier = file(outputFasta,open="w")
  
  for (i in 1:nrow(transposer))
  {
    cat(vecteur[i],file=fichier,sep = '\n')
    sample_seqs <- as.character(transposer[i,])
    sample_seqs = paste0(sample_seqs,collapse ='')
    if(nchar(sample_seqs)>500)
    {
      
      cat(substring(sample_seqs, 1, 500), file=fichier, sep = '\n')
      RestOfLine <- nchar(substring(sample_seqs, 501, nchar(sample_seqs)))
      if(RestOfLine >500)
      {
        cat(substring(sample_seqs, 501, 1000), file=fichier, sep = '\n')
        RestOfLine2 <- nchar(substring(sample_seqs, 1001, nchar(sample_seqs)))
        if(RestOfLine2 >500)
        {
          cat(substring(sample_seqs, 1001, 1500), file=fichier, sep = '\n')
          RestOfLine3 <- nchar(substring(sample_seqs, 1501, nchar(sample_seqs)))
          if(RestOfLine3 >500)
          {
            cat(substring(sample_seqs, 1501, 2000), file=fichier, sep = '\n')
            cat(substring(sample_seqs, 2001, nchar(sample_seqs)), file=fichier, sep = '\n')
          }
          else
          {
            cat(substring(sample_seqs, 1501, nchar(sample_seqs)), file=fichier, sep = '\n')
          }
        }
        else
        {
          cat(substring(sample_seqs, 1001, nchar(sample_seqs)), file=fichier, sep = '\n')
        }
        
        
      }
      else
      {
        cat(substring(sample_seqs, 501, nchar(sample_seqs)), file=fichier, sep = '\n')
      }
      
    }
    else
      cat(sample_seqs, file = fichier, sep = '\n')
  }
  
  close(fichier)
}

fasta_from_haplotype(haplotype, names(haplotype), outputFasta)

# Use this link to convert the Fasta to Nexus Format
# http://phylogeny.lirmm.fr/phylo_cgi/data_converter.cgi
#===================
# CREATE TRAIT FILE
common <- metadata[metadata$sample_ID %in% names(haplotype),]

trait.labels <- common %>% 
  select(species) %>% 
  unique() %>% 
  pull() %>% 
  sort()

Matrix <- data.frame(INDV = common$sample_ID,
                     arabiensis = 0,
                     coluzzii = 0,
                     gambiae = 0, 
                     S_M = 0)

for (t in 1:nrow(Matrix)) {
  cat("Computing SNPs ", t, "\n")
  for (c in 2:ncol(Matrix)) {
    if(names(Matrix)[c] == common$species[t]) Matrix[t, c] <- 1
  }
}

trait.file <- file("../results/hpx15_trait.txt", open = "w")
cat("Begin Traits;", file = trait.file, sep = "\n")
cat(paste0("Dimensions NTraits=", ncol(Matrix)-1,";"), file = trait.file, append = TRUE, sep = "\n")
cat("Format labels=yes missing=? separator=Comma;", file = trait.file, append = TRUE, sep = "\n")
cat("TraitLabels ", names(Matrix)[-1], ";", file = trait.file, append = TRUE)
cat("\nMatrix ", file = trait.file, append = TRUE, sep = "\n")

Matrix$New <- paste(Matrix$arabiensis, Matrix$coluzzii, Matrix$gambiae, Matrix$S_M, sep = ',')

write.table(Matrix[c(1,6)], file = trait.file, append = TRUE, sep = "\t",
            col.names = FALSE, quote = FALSE, row.names = FALSE)
cat(";", file = trait.file, append = TRUE, sep = "\n")
cat("END", file = trait.file, append = TRUE, sep = "\n")
close(trait.file)
