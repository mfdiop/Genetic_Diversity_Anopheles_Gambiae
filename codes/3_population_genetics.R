
# ----------------------------
# 3_population_genetics.R
#
# Author: Mouhamadou F. DIOP
# Date: 2022-07-29
#
# Purpose:
# Estimate Population Genetic.
# ----------------------------
rm(list = ls())

FILTERED_VCF <- commandArgs(trailingOnly = TRUE)

if(length(FILTERED_VCF) == 0){
    print("===================================================", quote = FALSE)
    print("        Running Population Genetic Analysis        ", quote = FALSE)
    print("===================================================", quote = FALSE)
    print("", quote = FALSE)
    print(" Usage: Rscript 3_population_genetics.R file.vcf.gz", quote = FALSE)
    stop(" Not enough arguments supplied ...")
    
}

# Uncomment these lines to install packages
# install.packages("devtools")
# devtools::install_github("jgx65/hierfstat")

library(data.table)
library(hierfstat)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(tidyverse)
library(adegenet)
library(reshape2)

## LOAD VCF
# FILTERED_VCF <- "data/process/HPX15.vcf.gz"

## LOAD METADATA
metadata <- readxl::read_xlsx("data/metadata/global_Metadata.xlsx")

## Minor allele freq file
genename <- gsub(".vcf.gz", "", basename(FILTERED_VCF))

output_dir <- paste0("results/figures/", genename)
if(!dir.exists(output_dir)) dir.create(output_dir)

# Modify ID column in VCF file
system(paste0("bcftools annotate --set-id +'%CHROM:%POS' ", 
              FILTERED_VCF, " > ", gsub(".gz", "", FILTERED_VCF)))

# Index output file
system(paste0("bgzip -f ", gsub(".gz", "", FILTERED_VCF)))

#==================
# Quality Control
#==================
geno = 0.01; maf = 0.01
mind = 0.01; hwe = 0.05
LD.window = 50; LD.step = 5; LD.correlation = 0.25
plink <- "/home/karim/Documents/Mes_Programmes/plink_linux_x86_64_20181202/plink"

data.out <- paste0("results/tables/", genename, "/", genename)

system(paste0(plink , " --vcf  ", FILTERED_VCF,
              " --hwe ", hwe,
              " --mind ", mind,
              " --geno ", geno,
              " --maf ", maf,
              " --make-bed --allow-extra-chr",
              " --out ", data.out, "_qc"))

system(paste0(plink , " --bfile  ", data.out,
              "_qc -indep-pairwise ", LD.window, " ",
              LD.step, " ", LD.correlation,
              " --allow-extra-chr",
              " --out ",data.out))

system(paste0(plink , " --bfile ", data.out,
              "_qc --extract  ", data.out, ".prune.in",
              " --recode vcf --allow-extra-chr ", 
              "--out ", data.out))

#===================
# Read filtered VCF
#===================
gene.VCF <- read.vcfR(paste0(data.out, ".vcf"), verbose = FALSE)

#===================
# Convert VCF 
# to a genind object
#===================
gen.gene <- vcfR2genind(gene.VCF, return.alleles = TRUE, NA.char = "./.") 

# Extract IDs from genind object
indiv <- indNames(gen.gene)

# Remove duplicates
indiv <- gsub('_.*', '', indiv)

# Extract Species
pop <- metadata[which(metadata$sample_ID %in% indiv), "species"] %>% pull()
pop(gen.gene) <- as.factor(pop)

#===========================
# Calculate heterozygosity 
# per site
#==========================

# Calculate basic stats using hierfstat
basic_gene = basic.stats(gen.gene, diploid = TRUE)

# Mean observed heterozygosity
Ho_gene = apply(basic_gene$Ho, MARGIN = 2, FUN = mean, na.rm = FALSE) %>%
  round(digits = 3)

# Mean expected heterozygosity
He_gene = apply(basic_gene$Hs, MARGIN = 2, FUN = mean, na.rm = FALSE) %>%
  round(digits = 2)

#===========================
# Visualize heterozygosity
#===========================

# Create a data.frame of site names, Ho and He and then convert to long format
Het_gene_df = data.frame(Site = names(Ho_gene), Ho = Ho_gene, He = He_gene) %>%
  melt(id.vars = "Site")

# Custom theme for ggplot2
custom_theme = theme(
  axis.text.x = element_text(size = 10, angle = 35, vjust = 0.5, face = "bold"),
  axis.text.y = element_text(size = 10, face = "bold"),
  axis.title.y = element_text(size = 12, face = "bold"),
  axis.title.x = element_blank(),
  axis.line.y = element_line(size = 0.5),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, size = 15, face="bold")
)

# Italic label
hetlab.o <- expression(italic("H")[o])
hetlab.e <- expression(italic("H")[e])

ylimits <- c(0, max(Het_gene_df$value)+0.01)

ggplot(data = Het_gene_df, aes(x = Site, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  scale_y_continuous(expand = c(0,0), limits = ylimits) + 
  scale_fill_manual(values = c("royalblue", "#bdbdbd"), labels = c(hetlab.o, hetlab.e)) +
  ylab("Heterozygosity") +
  ggtitle(genename) +
  custom_theme

ggsave(paste0(output_dir, "/heterozygosity.pdf"),
       width = 12, height = 8, dpi = 600)

#==========================
# Compute pairwise FST 
# (Weir & Cockerham 1984)
#==========================
gene_fst <- genet.dist(gen.gene, method = "WC84")

#=========================
# Visualize pairwise FST.
#=========================

# Convert dist object to data.frame
fst.matrix <- as.matrix(gene_fst)

# Save Fst data
# write.table(fst.matrix, paste0("results/tables/", genename, "/Fst.xlsx"),
#             col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')

## Sort column names
fst.matrix <- fst.matrix[order(rownames(fst.matrix)), order(colnames(fst.matrix))]
ind <- which( upper.tri(fst.matrix), arr.ind = TRUE)
fst.df <- data.frame(Site1 = dimnames(fst.matrix)[[2]][ind[,2]],
                    Site2 = dimnames(fst.matrix)[[1]][ind[,1]],
                    Fst = fst.matrix[ ind ] %>% round(digits = 3))

# Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] <- 0

# Fst italic label
fst.label <- expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid <- max(fst.df$Fst) / 2

# Plot heatmap
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red",
                       midpoint = mid, name = fst.label, # , breaks = c(0, 0.15, 0.25, 0.35, 0.45)
                       limits = c(0, max(fst.df$Fst)))+ 
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)
  )

ggsave(paste0(output_dir, "/Fst.pdf"),
       width = 12, height = 8, dpi = 600)

#================
# Perform DAPC
#================
set.seed(123)

# Replace missing data with the mean allele frequencies
x <- tab(gen.gene, NA.method = "mean")

#==============
# Perform PCA
#==============

pca <- dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)

# Analyse how much percent of genetic variance is explained by each PC
percent <- pca$eig/sum(pca$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)",
        ylim = c(0,12), names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates
pca_coords <- as.data.frame(pca$li)

#=======================
# Visualize PCA results
#=======================

# Rename columns of dataframe
colnames(pca_coords) <- c("PC1","PC2","PC3")

# Add a column containing individuals
pca_coords$Ind <- indiv

# Add a column with the site IDs
pca_coords$Site <- gen.gene$pop

# Calculate centroid (average) position for each population
centroid <- aggregate(cbind(PC1, PC2, PC3) ~ Site,
                      data = pca_coords,
                      FUN = mean)

# Add centroid coordinates to pca_coords dataframe
pca_coords <- left_join(pca_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define colour palette
# cols <- brewer.pal(nPop(gen.gene), "Set1")
cols <- brewer.pal(nPop(gen.gene), "Set2")

# Custom x and y labels
xlab <- paste("PC 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab <- paste("PC 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom theme for ggplot2
ggtheme <- theme(axis.text = element_text(colour = "black", size = 12, face = "bold"),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour = "black", size = 12, face = "bold"),
                axis.line = element_line(size = 0.5),
                panel.border = element_rect(colour = "black", fill = NA, size = 1),
                panel.background = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 15, face = "bold") 
)

# Scatter plot PC 1 vs. 2
ggplot(data = pca_coords, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    
    # spider segments
    geom_segment(aes(xend = PC1.cen, yend = PC2.cen, colour = Site), show.legend = FALSE) +
    
    # points
    geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE) +
    
    # centroids
    geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE) +
    
    # colouring
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    
    # custom labels
    labs(x = xlab, y = ylab) +
    ggtitle(paste0(genename, " PCA")) +
    
    # custom theme
    ggtheme

# Export plot
ggsave(paste0(output_dir, "/PCA.pdf"),
       width = 12, height = 8, dpi = 600)


# Perform cross validation to find the optimal number of PCs to retain in DAPC
crossval <- xvalDapc(x, gen.gene$pop, result = "groupMean", xval.plot = TRUE)

# Number of PCs with best stats (lower score = better)
numPCs <- as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)

# Run a DAPC using population IDs as priors
dapc1 <- dapc(gen.gene, gen.gene$pop, n.pca = numPCs, n.da = 3)

scatter(dapc1, col = cols) # plot of the group

# Analyse how much percent of genetic variance is explained by each PC
peig <- dapc1$eig/sum(dapc1$eig)*100
barplot(peig, ylab = "Percent of genetic variance explained by eigenvectors", 
        names.arg = round(peig, 2))

#========================
# Visualize DAPC results.
#========================
# Create a dataframe containing individual coordinates
dapc_coords <- as.data.frame(dapc1$ind.coord)

# Rename columns of dataframe
colnames(dapc_coords) <- c("PC1","PC2","PC3")

# Add a column containing individuals
dapc_coords$Ind <- indiv

# Add a column with the site IDs
dapc_coords$Site <- gen.gene$pop

# Calculate centroid (average) position for each population
centroid <- aggregate(cbind(PC1, PC2, PC3) ~ Site,
                      data = dapc_coords,
                      FUN = mean)

# Add centroid coordinates to ind_coords dataframe
dapc_coords <- left_join(dapc_coords, centroid, by = "Site", suffix = c("",".cen"))

# Custom x and y labels
xlab <- paste("PC 1 (", format(round(peig[1], 1), nsmall=1)," %)", sep="")
ylab <- paste("PC 2 (", format(round(peig[2], 1), nsmall=1)," %)", sep="")

# Scatter plot PC1 vs. PC2
ggplot(data = dapc_coords, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    
    # spider segments
    geom_segment(aes(xend = PC1.cen, yend = PC2.cen, colour = Site), show.legend = FALSE) +
    
    # points
    geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE) +
    
    # centroids
    geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE) +
    
    # colouring
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    
    # custom labels
    labs(x = xlab, y = ylab) +
    ggtitle(paste0(genename, " DAPC")) +
    
    # custom theme
    ggtheme

# Save plot
ggsave(paste0(output_dir, "/DAPC.pdf"),
       width = 12, height = 8, dpi = 600)

system(paste0("rm -rf ", data.out, "*"))
