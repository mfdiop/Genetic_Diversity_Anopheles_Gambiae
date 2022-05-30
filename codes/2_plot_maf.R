
rm(list = ls())
# generate_maf_plot.R
#
# Author: Mouhamadou F. DIOP
# Date: 2021-11-20
#
# Purpose:
# Produce Minor Allele Frequency for each species.
#
# ------------------------------------------------------------------
arguments <- commandArgs(trailingOnly = TRUE)

if(length(arguments) == 0){
  cat("=============================================\n")
  cat("SCRIPT SHOULD BE RUN LIKE: \n")
  cat("Rscript 2_plot_maf.R gene_maf.tsv PerSpecies_maf.tsv")
  cat("\n")
  stop("NOT ENOUGH ARGUMENTS SUPPLIED .....")
  cat("=============================================\n")
}

# NOTE - uncomment these lines to install packages as needed
devtools::install_github("thomasp85/patchwork")

# Load packages
library(tidyverse)
library(grid)
library(patchwork)

# all_maf <- read_tsv("../results/tables/HPX15_MAF.tsv")
all_maf <- arguments[1]

genename <- gsub("_MAF.tsv", "", basename(all_maf))
all_maf <- read_tsv(all_maf)

## CREATE HISTOGRAM FOR EACH SPECIES
barfill <- "#4271AE" # "#56B4E9" "gold1"
barlines <- "#1F3552"

All <- ggplot(all_maf, aes(x = maf)) +
  geom_histogram(aes(y = ..count..), colour = barlines, fill = "#666666") +
  scale_x_continuous(name = "Minor Allele Frequency") +
  scale_y_continuous(name = "Frequency") +
  ggtitle("All") +
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 12, family = "Tahoma", face = "bold", hjust = 0.5),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 7),
        axis.text.y=element_text(colour="black", size = 7))

## Load maf file
# maf <- read_tsv("../results/tables/HPX15_PerSpecies_MAF.tsv")
maf <- read_tsv(arguments[2])

## EXTRACT MAF FOR EACH SPECIES
arabiensis <- maf %>% filter(Species == "Arabiensis")
coluzzi <- maf %>% filter(Species == "Coluzzii")
gambiae <- maf %>% filter(Species == "Gambiae")
sm <- maf %>% filter(Species == "S_M")


p1 <- ggplot(arabiensis, aes(x = maf)) +
  geom_histogram(aes(y = ..count..), colour = barlines, fill = "#00FF00") +
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "Frequency") +
  ggtitle("Arabiensis") +
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 9, family = "Tahoma", face = "bold", hjust = 0.5),
        text = element_text(family="Tahoma"),
        axis.text.x = element_text(colour="black", size = 5),
        axis.text.y = element_text(colour="black", size = 5),
        axis.title = element_text(size = 7))

p2 <- ggplot(coluzzi, aes(x = maf)) +
  geom_histogram(aes(y = ..count..), colour = barlines, fill = "#0000FF") +
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "") +
  ggtitle("Coluzzi") +
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 9, family = "Tahoma", face = "bold", hjust = 0.5),
        text = element_text(family="Tahoma"),
        axis.text.x = element_text(colour="black", size = 5),
        axis.text.y = element_text(colour="black", size = 5),
        axis.title = element_text(size = 7))

p3 <- ggplot(gambiae, aes(x = maf)) +
  geom_histogram(aes(y = ..count..), colour = barlines, fill = "#FF0000") +
  scale_x_continuous(name = "Minor Allele Frequency") +
  scale_y_continuous(name = "Frequency") +
  ggtitle("Gambiae") +
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 9, family = "Tahoma", face = "bold", hjust = 0.5),
        text = element_text(family="Tahoma"),
        axis.text.x = element_text(colour="black", size = 5),
        axis.text.y = element_text(colour="black", size = 5),
        axis.title = element_text(size = 7))


p4 <- ggplot(sm, aes(x = maf)) +
  geom_histogram(aes(y = ..count..), colour = barlines, fill = "#4271AE") +
  scale_x_continuous(name = "Minor Allele Frequency") +
  scale_y_continuous(name = "") +
  ggtitle("S_M") +
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 9, family = "Tahoma", face = "bold", hjust = 0.5),
        text = element_text(family="Tahoma"),
        axis.text.x = element_text(colour="black", size = 5),
        axis.text.y = element_text(colour="black", size = 5),
        axis.title = element_text(size = 7))


All + wrap_plots(p1, p2, p3, p4)

# save to file
ggsave(paste0("../results/figures/", genename, "_maf.png"))

# #========================
# ## GROUPED BARCHART PLOT
# #========================
# 
# # 1. Create a vector of break points:
# b <- c(-Inf, 0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
# 
# # 2. Create a vector of names for break points:
# names <- c("0", "0 - 0.01", "0.01 - 0.025", "0.025 - 0.05", 
#            "0.05 - 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.3 - 0.4", "0.4 - 0.5")
# 
# # 3. Create a vector of names for break points:
# maf$maf_bin <- cut(maf$maf, breaks = b, labels = names)
# 
# # Grouped
# ggplot(maf, aes(fill = Species, x = maf_bin)) + 
#   geom_bar(position = "dodge", stat = "count") +
#   scale_fill_manual(values = c("#00FF00", "#0000FF", "#FF0000", "#4271AE")) +
#   xlab("Minor Allele Frequency") + ylab("Frequency") + # ggtitle("") + "green", "blue", gold1, steelblue,  grey
#   theme_bw() +
#   theme(axis.line = element_line(size = 0.5, colour = "black"),
#         panel.grid.major = element_line(colour = "#d3d3d3"),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(), panel.background = element_blank(),
#         # plot.title = element_text(size = 9, family = "Tahoma", face = "bold", hjust = 0.5),
#         text = element_text(family="Tahoma", face = "bold"),
#         axis.text.x = element_text(colour="black", size = 5),
#         axis.text.y = element_text(colour="black", size = 5),
#         axis.title = element_text(size = 7))
# 
# # save to file
# ggsave("../results/figures/minor_allele_freq2.png")

