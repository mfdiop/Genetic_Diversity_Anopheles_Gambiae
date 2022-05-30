
arguments <- commandArgs(trailingOnly = TRUE)

if(length(arguments) == 0){
    cat("==========================================================\n")
    cat("SCRIPT SHOULD BE RUN LIKE: \n")
    cat("Rscript pop_struct.R file.vcf  output metadata.txt ")
    cat("\n")
    stop("NOT ENOUGH ARGUMENTS SUPPLIED .....")
    cat("==========================================================\n")
}

devtools::install_github("bobverity/bobfunctions2", ref = "version1.0.0")

load.lib <- c("data.table", "tidyverse", "RColorBrewer", "gridExtra")


install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib, require, character=TRUE)

require(data.table)
require(tidyverse)
require(RColorBrewer)
require(bobfunctions2)

plink <- "/home/karim/Documents/Mes_Programmes/plink_linux_x86_64_20181202/plink"

# data.in = arguments[1];
# data.out = arguments[2];
data.in = "../data/raw/all_HPX15.vcf.gz";
data.out = "../data/proceed/hpx15";
geno = 0.01;
mind = 0.01;
LD.window = 50;
LD.step = 5;
LD.correlation = 0.2

# metadata <- read_tsv(arguments[3])
metadata <- readxl::read_xlsx("../data/metadata/metadata.xlsx")

metadata$sample_ID <- gsub("-.*", "", metadata$sample_ID)

system(paste0(plink, " --vcf  ", data.in, " --allow-extra-chr  --geno ", geno,
              " --mind ", mind, " --make-bed --out ", data.out,"_qc"))

system(paste0(plink, " --bfile ", data.out, "_qc --allow-extra-chr -indep-pairwise ",
              LD.window, " ", LD.step, " ", LD.correlation, " --out ", data.out))

system(paste0(plink, " --bfile ", data.out, "_qc --allow-extra-chr  --extract ",
              data.out, ".prune.in ", " --make-bed --out ", data.out, "_qc2"))

### PCA USING PLINK WITH --pca 
system(paste0(plink, " --bfile ", data.out, "_qc2 --allow-extra-chr ",
              "--pca 10 header tabs --out plink_pca"))

## Load Plink output
pca <- read_tsv('plink_pca.eigenvec')
pca$FID <- gsub("-.*", "", pca$FID)

metadata <- metadata[which(pca$FID %in% metadata$sample_ID), ]

pca$IID <- metadata$countries
names(pca)[2] <- 'Country'

pca_percentage <- fread('plink_pca.eigenval', header = FALSE)

# make axis labels
x_lab <- sprintf("PC1 (%s%%)", round(pca_percentage$V1[1], digits = 2))
y_lab <- sprintf("PC2 (%s%%)", round(pca_percentage$V1[2], digits = 2))
z_lab <- sprintf("PC3 (%s%%)", round(pca_percentage$V1[3], digits = 2))

# plotting arguments
pointsize_2D <- 1
pointsize_3D <- 0.3
panel_letter_size <- 8

# change order of countries. Changes order in which points are plotted so that DRC is behind
w <- match(pca$Country, c( "Angola", "Burkina", "Cameroon", "Centr_Africa", "Equatorial_Guinea",      
                           "Gabon", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Ivory_Coast", "Kenya",    
                           "Malawi", "Mali", "Mayotte", "Mozambique", "RD_Congo", "Tanzania", "Uganda"))
pca <- pca[order(w),]

# define colour legend
# col_pal <- c(grey(0.8), RColorBrewer::brewer.pal(8, "Dark2")[c(2,1,3,4)])    # "#0000FF", "#FF0000"
col_pal <- c(grey(0.8), "#00FF00", "#0000FF", "#66A61E", "#E6AB02",
             "#A6761D", "#666666", RColorBrewer::brewer.pal(12, "Paired"))

# produce 2D scatterplot
plot1 <- ggplot() + theme_bw(base_size = 8) +
    geom_point(aes(x = PC1, y = PC3, col = Country), data = pca, size = pointsize_2D) + 
    xlab(x_lab) + ylab(y_lab) + ggtitle("a)") + 
    theme(plot.title = element_text(size = panel_letter_size)) +
    scale_color_manual(values = col_pal, name = "Country")

# produce 3D scatterplot
plot2 <- gg3d_scatterplot(pca$PC1, pca$PC2, pca$PC3, colour = pca$Country, size = pointsize_3D, theta = 135,
                          d = 2.0, x_lab = x_lab, y_lab = y_lab, z_lab = z_lab,
                          axis_lab_dist = 1.5, axis_lab_size = 1.8, grid_size = 0.15, zero_line_size = 0.4) +
    scale_color_manual(values = col_pal, name = "Country") +
    ggtitle("b)") + theme(plot.title = element_text(size = panel_letter_size),
                          legend.text = element_text(size = 8),
                          legend.title = element_text(size = 8))

# plot side-by-side
plot3 <- gridExtra::grid.arrange(plot1, plot2, ncol = 2, widths = c(2.1,3))

plot(plot3)

# save to file
file_ext <- c("jpeg", "pdf", "png")
for (i in seq_along(file_ext)) {
    ggsave(sprintf("figure1_PCA/figure1_PCA.%s", file_ext[i]),
           plot = plot3, device = file_ext[i], width = 179, height = 80, units = "mm")
}