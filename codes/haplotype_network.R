
# Clear workspace
rm(list = ls())

# Load library
library(tidyverse)

plink <- "/home/karim/Documents/Mes_Programmes/plink_linux_x86_64_20181202/plink"
vcf_file <-  "../Fatima/data/references/HPX14_maf_filtered.vcf.gz"
annotation_file <-  "../Fatima/data/references/HPX14_annotation.tab"

system(paste0("cat < \\(zgrep '#CHROM' ", vcf_file," | cut -f 1-3\\) < \\(paste < \\(zgrep -v '#' ", vcf_file,
              "| cut -f 1,2\\) < \\(zgrep -v '#' ", vcf_file, " | cut -f 1,2 | sed 's/\t/_/g'\\)\\) > ", annotation_file))

system(paste0(bgzip annotation_file.tab))
system(paste0(tabix -p vcf annotation_file.tab.gz))
system(paste0(bcftools annotate -c CHROM,POS,ID -a annotation_file.tab.gz $VCF_in -o $VCF_out))


## Set Variables
vcf_file <-  "../Fatima/results/HPX14_annot.vcf"
Output_dir <- "../Fatima/results"

mind <- 0.1; geno <- 0.1; hwe <- 0.0000001
maf <- 0.05; ld_window <- 100; ld_window_kb <- 10


Output <- file.path(Output_dir, gsub(".vcf", "_afterQC", basename(vcf_file)))

# quality control
# select one chromosome
system(paste0(plink, " --vcf ", vcf_file, " --mind ", mind, " --geno ", geno, 
              " --maf ", maf, " --hwe ", hwe, " --make-bed --allow-extra-chr --out ", Output))

# LD pruning -  remove SNPs with high LD with each other (removes one from each pair)

# replace --nonfounders with --make-founders!
system(paste0(plink, " --bfile ", Output, " --allow-extra-chr --indep-pairwise 100 10 0.5 --out ", Output))

# keep only selected markers for your data
system(paste0(plink, " --bfile ", Output," --allow-extra-chr --exclude ", Output, 
              ".prune.out --make-bed --out ", Output, "_Prunned"))

# IDB and IBS calculations
system(paste0(plink, " --bfile ", Output," --allow-extra-chr --exclude ", Output, 
              ".prune.out --genome --recode --out ", Output))

ibd <- read_table(paste0(Output, ".genome"))       

ibd <- ibd %>% 
    filter(!is.na(RATIO))

id1 <- ibd %>% select(FID1) %>% unique() %>% pull()
id2 <- ibd %>% select(FID2) %>% unique() %>% pull()

c(id1, id2) %>% unique() %>% length()

common <- intersect(id1, id2)
write.table(common, "../Fatima/results/HPX14_unrelated_isolates.txt", 
            col.names = F, row.names = F, quote = F)

# keep only selected markers and isolates for the data
system(paste0(plink, " --bfile ", Output," --allow-extra-chr --exclude ", Output, 
              ".prune.out --recode vcf --out ", Output, "_forHN"))
