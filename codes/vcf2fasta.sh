#!/bin/sh

conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n vcf-kit \
  danielecook::vcf-kit=0.2.6 \
  "bwa>=0.7.17" \
  "samtools>=1.10" \
  "bcftools>=1.10" \
  "blast>=2.2.31" \
  "muscle>=3.8.31" \
  "primer3>=2.5.0"

conda activate vcf-kit

VCF=[1]
FILENAME=$(basename $VCF _maf_filtered.vcf.gz)
PATH="/media/Data/Data/Documents_Karim/Fadel/MRC_projects/Fatima/programs"
python=/usr/bin/python3

#vk phylo fasta $VCF > $FILENAME.fasta


$python ${PATH}/vcf2phylip/vcf2phylip.py 
-i HPX15_maf_filtered.vcf.gz 
--output-prefix HPX15 
-n --min-samples-locus 2000 
--resolve-IUPAC 
--phylip-disable