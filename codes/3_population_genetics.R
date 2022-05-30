
rm(list = ls())
# population_genetics.R
#
# Author: Mouhamadou F. DIOP
# Date: 2021-11-20
#
# Purpose:
# Estimate Population Genetic.
# ----------------------------
# packages <- c("hierfstat", "vcfR", "poppr", "ape", "adegenet")
# install.packages(setdiff(packages, rownames(installed.packages())))
# install.packages("devtools")
# library(devtools)
# install_github("jgx65/hierfstat")
# library("hierfstat")

library(hierfstat)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(tidyverse)
library(adegenet)
library(reshape2)

## LOAD VCF
FILTERED_VCF <- "data/process/HPX15.vcf.gz"

## LOAD METADATA
metadata <- readxl::read_xlsx("data/metadata/global_Metadata.xlsx")

## Minor allele freq file
genename <- gsub(".vcf.gz", "", basename(FILTERED_VCF))
maf_file <- data.table::fread(paste0("results/tables/", genename, "_MAF.tsv"), 
                              header = TRUE)

# GET POSITIONS BASED ON MAF > 0.01
snpsTokeep <- "data/process/snpsTokeep.txt"

maf_file %>% 
  filter(maf > 0.01) %>% 
  distinct(pos, .keep_all = TRUE) %>% 
  select(chr, pos) %>%
  write.table(snpsTokeep, col.names = FALSE,
              row.names = FALSE, quote = FALSE, sep = '\t')

#==========================================
#---- Extract SNPs from filtered vcf file
#==========================================
MAF_FILTERED_VCF <- paste0("data/references/", gsub(".vcf.gz", "_maf_filtered", basename(FILTERED_VCF)))
system(paste0("vcftools --gzvcf ", FILTERED_VCF,
              " --positions ", snpsTokeep, 
              " --recode --recode-INFO-all --out ", MAF_FILTERED_VCF))

file.remove(snpsTokeep)

system(paste0("mv ", MAF_FILTERED_VCF, ".recode.vcf ", MAF_FILTERED_VCF, ".vcf")) # Rename MAF_FILTERED_VCFput VCF file
system(paste0("bgzip ", MAF_FILTERED_VCF, ".vcf")) # compress VCF 
system(paste0("tabix ", MAF_FILTERED_VCF, ".vcf.gz")) # Index VCF

vcfname <- paste0(MAF_FILTERED_VCF, ".vcf.gz")

gene.VCF <- read.vcfR(vcfname, verbose = FALSE)

# Convert VCF to a genind object
gen.gene <- vcfR2genind(gene.VCF, return.alleles = TRUE, NA.char = "./.") 

indiv <- indNames(gen.gene)
pop <- metadata[which(metadata$sample_ID %in% indiv), "species"] %>% pull()

pop(gen.gene) <- as.factor(pop)

# Calculate heterozygosity per site
# Calculate basic stats using hierfstat
basic_gene = basic.stats(gen.gene, diploid = TRUE)

# Mean observed heterozygosity per site
Ho_gene = apply(basic_gene$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)

# Mean expected heterozygosity per site
He_gene = apply(basic_gene$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)

# Visualise heterozygosity per site
# Create a data.frame of site names, Ho and He and then convert to long format
Het_gene_df = data.frame(Site = names(Ho_gene), Ho = Ho_gene, He = He_gene) %>%
  melt(id.vars = "Site")

# Custom theme for ggplot2
custom_theme = theme(
  axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, face = "bold"),
  axis.text.y = element_text(size = 10),
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

ylimits <- c(0, max(Het_gene_df$value)+0.02)
ggplot(data = Het_gene_df, aes(x = Site, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  scale_y_continuous(expand = c(0,0), limits = ylimits) +                                 #c(0,0.06)
  scale_fill_manual(values = c("royalblue", "#bdbdbd"), labels = c(hetlab.o, hetlab.e)) +
  ylab("Heterozygosity") +
  ggtitle(genename) +
  custom_theme

ggsave(paste0("../results/figures/", genename, "_heterozygosity.pdf"),
       width = 12, height = 8, dpi = 600)

# ## Inbreeding coefficient (FIS)
# ## Calculate mean FIS per site.
# apply(basic_hpx15$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
#   round(digits = 3)

## FST Compute pairwise FST (Weir & Cockerham 1984).
gene_fst <- genet.dist(gen.gene, method = "WC84")

## Visualise pairwise FST for hpx15.
# Convert dist object to data.frame
fst.matrix <- as.matrix(gene_fst)
write.table(fst.matrix, paste0("../results/tables/", genename, "_Fst.xlsx"),
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')

## Sort column names
fst.matrix <- fst.matrix[order(rownames(fst.matrix)), order(colnames(fst.matrix))]
ind <- which( upper.tri(fst.matrix), arr.ind = TRUE)
fst.df <- data.frame(Site1 = dimnames(fst.matrix)[[2]][ind[,2]],
                    Site2 = dimnames(fst.matrix)[[1]][ind[,1]],
                    Fst = fst.matrix[ ind ] %>% round(digits = 3))

# Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] <- 0

# Print data.frame summary
fst.df %>% str

# Fst italic label
fst.label <- expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid = max(fst.df$Fst) / 2

# Plot heatmap
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, 
                       limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.15, 0.25, 0.35, 0.45))+ 
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

ggsave(paste0("../results/figures/", genename, "_Fst.pdf"),
       width = 12, height = 8, dpi = 600)

# Perform a PCA (Principle Components Analysis)
# Replace missing data with the mean allele frequencies
x <- tab(gen.gene, NA.method = "mean")

# Perform PCA
pca1 <- dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)

# Analyse how much percent of genetic variance is explained by each axis
percent <- pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))

## Visualise PCA results.
# Create a data.frame containing individual coordinates
ind_coords <- as.data.frame(pca1$li)

# Rename columns of dataframe
colnames(ind_coords) <- c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind <- indNames(gen.gene)

# Add a column with the site IDs
ind_coords$Site <- gen.gene$pop

# Calculate centroid (average) position for each population
centroid <- aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords <- left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define colour palette
cols <- brewer.pal(nPop(gen.gene), "Set1")

# Custom x and y labels
xlab <- paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab <- paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) 
)

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE) +
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
ggsave(paste0("../results/figures/", genename, "_PCA.pdf"),
       width = 12, height = 8, dpi = 600)

file.remove("./Rplots.pdf")
