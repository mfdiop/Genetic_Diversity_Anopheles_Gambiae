

rm(list = ls())
## Create phylogeny tree for HPX15 by converting the VCF file to a genind object
## and run hclust on the distance matrix. Save the tree as a Newick format
## Load the tree into ITOL to draw it
## https://grunwaldlab.github.io/Population_Genetics_in_R/Pop_Structure.html
## https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
## Let’s first load the libraries needed for analysis:

library(dendsort)
library(dendextend)
library(ctc)
library(tidyverse)
library(data.table)
library(vcfR)
library(poppr)
library(ape) # To visualize the tree using the "nj" function
library(RColorBrewer)
library(magrittr)
# library(mmod)

## LOAD VCF
VCF <- "../data/references/MISO_maf_filtered.vcf.gz"
genename <- gsub("_maf_filtered.vcf.gz", "", basename(VCF))

## LOAD METADATA
metadata <- readxl::read_xlsx("../data/metadata/global_Metadata.xlsx")

cat('===============================\n')
cat("  Load filtered VCF file \n")
cat('===============================\n')

gene.VCF <- read.vcfR(VCF, verbose = FALSE)

# Convert VCF to a genind object
gen.gene <- vcfR2genind(gene.VCF, return.alleles = TRUE, NA.char = "./.") 

indiv <- indNames(gen.gene)
pop <- metadata[which(metadata$sample_ID %in% indiv), "species"] %>% pull()

pop(gen.gene) <- as.factor(pop)

## Note
# prevosti.dist() and diss.dist() are exactly the same, 
# but diss.dist() scales better for large numbers of individuals (n > 125) 
# at the cost of required memory. So we are using diss.dist

# However, if the user supplies a [genind] or
# [genclone]object, [prevosti.dist()] will be used for
# calculation.

hpx15.dist <- prevosti.dist(gen.gene)
my_hclust <- hclust(hpx15.dist, method = "complete")
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

my_hclust_sort <- sort_hclust(my_hclust)

plot(my_hclust_sort, main = "Sorted Dendrogram", xlab = "", sub = "")

write(hc2Newick(my_hclust_sort), paste0("../results/", genename, ".newick"))

#=======================================
## CREATE THE ANNOTATIONS FILE FOR ITOL
#=======================================
## SPECIES ANNOTATION
colors <- c("#00FF00", "#0000FF", "#66A61E", "#E6AB02")

annotation_species <- data_frame(NODE_ID = metadata$sample_ID,
                              TYPE = "range",
                              COLOR = metadata$species,
                              LABEL = toupper(metadata$species))

annotation_species$COLOR <- colors[match(annotation_species$COLOR, unique(metadata$species))]
write.table(annotation_species, paste0("../results/", genename, ".ann.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

# ## POPULATION ANNOTATION
# colors <- c(grey(0.8), "#00FF00", "#0000FF", "#66A61E", "#E6AB02",
#             "#A6761D", "#666666", RColorBrewer::brewer.pal(12, "Paired"),
#             RColorBrewer::brewer.pal(8, 'Spectral'))
# annotation_pop <- data_frame(NODE_ID = metadata$sample_ID,
#                              TYPE = "branch",
#                              COLORS = metadata$country_date,
#                              LABEL = "normal",
#                              SIZE= 2 )
# 
# annotation_pop$COLORS <- colors[match(annotation_pop$COLORS, unique(metadata$country_date))]
# 
# write.table(annotation_pop, paste0("../results/", genename, ".ann_pop.txt"),
#             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

my_hclust <- hclust(hpx15.dist, method = "complete")

my_dendogram <- as.dendrogram(my_hclust)
cols <- brewer.pal(n = nPop(gen.gene), name = "Dark2")

dendogram <- my_dendogram %>%
  # set("labels_col", value = c("skyblue", "orange", "grey", "red")) %>%
  set("branches_k_color", value = c("skyblue", "orange", "grey", "red")) %>%
  set("labels_cex", 0.2)

circlize_dendrogram(dendogram, facing = "outside",
                    labels = TRUE,
                    labels_track_height = 0.1,
                    dend_track_height = 0.7)
