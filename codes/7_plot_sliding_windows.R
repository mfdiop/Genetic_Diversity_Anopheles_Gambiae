
rm(list = ls())
# plot_sliding_windows.R
#
# Author: Mouhamadou F. DIOP
# Date: 2021-11-20
#
# Purpose:
# Produce Minor Allele Frequency for each species.
#
# ------------------------------------------------------------------

# library
library(tidyverse)
library(ggplot2)
library(hrbrthemes)

cols <- c(grey(0.8), "#00FF00", "#0000FF", "#66A61E", "#E6AB02",
          "#A6761D", "#666666", RColorBrewer::brewer.pal(12, "Paired"),
          RColorBrewer::brewer.pal(11, 'Spectral'), RColorBrewer::brewer.pal(8, 'Dark2'))

#==========================
#  NUCLEOTIDE DIVERSITY
#==========================
Filename1 <- "../results/tables/sliding_windows/Vg_Arabiensis_pi.tsv"; 
# species1 <- gsub("\\_[.Aa-zZ]*", "", basename(Filename1))
Filename2 <- "../results/tables/sliding_windows/Vg_Coluzzii_pi.tsv";
# species2 <- gsub("\\_[.Aa-zZ]*", "", basename(Filename2))
Filename3 <- "../results/tables/sliding_windows/Vg_Gambiae_pi.tsv";
# species3 <- gsub("\\_[.Aa-zZ]*", "", basename(Filename3))
Filename4 <- "../results/tables/sliding_windows/Vg_S_M_pi.tsv";
# species4 <- gsub("\\_[.Aa-zZ]*", "", basename(Filename4))

species1 <- "Arabiensis"; species2 <- "Coluzzii"; 
species3 <- "Gambiae"; species4 <- "SM"

pi_wide1 <- read_tsv(Filename1)
names(pi_wide1) <- paste(colnames(pi_wide1), species1, sep = " ")
pi_long1 <- gather(pi_wide1, Populations, Nuc_Div, factor_key=TRUE)

pi_wide2 <- read_tsv(Filename2)
names(pi_wide2) <- paste(colnames(pi_wide2), species2, sep = " ")
pi_long2 <- gather(pi_wide2, Populations, Nuc_Div, factor_key=TRUE)

pi_wide3 <- read_tsv(Filename3)
names(pi_wide3) <- paste(colnames(pi_wide3), species3, sep = " ")
pi_long3 <- gather(pi_wide3, Populations, Nuc_Div, factor_key=TRUE)

pi_wide4 <- read_tsv(Filename4)
names(pi_wide4) <- paste(colnames(pi_wide4), species4, sep = " ")
pi_long4 <- gather(pi_wide4, Populations, Nuc_Div, factor_key=TRUE)

pi_long <- rbind(pi_long1, pi_long2, pi_long3, pi_long4)

pi_long %>% 
  filter(Nuc_Div != 0) %>% 
  write.table("../results/tables/Vg_pi.xlsx", sep = '\t',
              col.names = TRUE, row.names = FALSE, quote = FALSE)

png(paste0('../results/figures/Vg_ND.png'),width = 12, height = 6, 
    family = 'serif', units = 'in', res = 600)

p <- pi_long %>%
  mutate(Populations = fct_reorder(Populations, Nuc_Div)) %>% # Reorder data
  ggplot( aes(x = Populations, y = Nuc_Div)) + # , fill=Populations, color=Populations
  geom_boxplot(lwd = 0.3, fill = cols) +
  coord_flip() +
  theme_bw() +
  theme_ipsum() +
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8, face = 'bold')) +
  xlab("") + ylab("Nucleotide Diversity")

print(p)
dev.off()

file.remove(Filename1, Filename2, Filename3, Filename4)

#==============
#  TAJIMA'S D
#==============
cols <- c(grey(0.8), "#00FF00", "#0000FF", "#66A61E", "#E6AB02",
          "#A6761D", "#666666", RColorBrewer::brewer.pal(12, "Paired"),
          RColorBrewer::brewer.pal(11, 'Spectral'), RColorBrewer::brewer.pal(7, 'Dark2'))

Filename1 <- "../results/tables/sliding_windows/Vg_Arabiensis_tajima.tsv"; 
Filename2 <- "../results/tables/sliding_windows/Vg_Coluzzii_tajima.tsv";
Filename3 <- "../results/tables/sliding_windows/Vg_Gambiae_tajima.tsv";
Filename4 <- "../results/tables/sliding_windows/Vg_S_M_tajima.tsv"

pi_wide1 <- read_tsv(Filename1)
names(pi_wide1) <- paste(colnames(pi_wide1), species1, sep = " ")
pi_long1 <- gather(pi_wide1, Populations, Tajima, factor_key=TRUE)

pi_wide2 <- read_tsv(Filename2)
names(pi_wide2) <- paste(colnames(pi_wide2), species2, sep = " ")
pi_long2 <- gather(pi_wide2, Populations, Tajima, factor_key=TRUE)

pi_wide3 <- read_tsv(Filename3)
names(pi_wide3) <- paste(colnames(pi_wide3), species3, sep = " ")
pi_long3 <- gather(pi_wide3, Populations, Tajima, factor_key=TRUE)

pi_wide4 <- read_tsv(Filename4)
names(pi_wide4) <- paste(colnames(pi_wide4), species4, sep = " ")
pi_long4 <- gather(pi_wide4, Populations, Tajima, factor_key=TRUE)

pi_long <- rbind(pi_long1, pi_long2, pi_long3, pi_long4)

pi_long %>% 
  drop_na() %>% 
  write.table("../results/tables/Vg_tajima.xlsx", sep = '\t',
              col.names = TRUE, row.names = FALSE, quote = FALSE)

png(paste0('../results/figures/Vg_tajima.png'), width = 12, height = 6, 
    family = 'serif', units = 'in', res = 600)

p <- pi_long %>%
  mutate(Populations = fct_reorder(Populations, Tajima)) %>% # Reorder data
  drop_na() %>% 
  ggplot( aes(x = Populations, y = Tajima)) +
  geom_boxplot(lwd = 0.3, fill = cols) +
  geom_hline(yintercept = 0, linetype="dashed", color = grey(0.8)) +
  
  # This switch X and Y axis and allows to get the horizontal version
  coord_flip() +
  theme_bw() +
  theme_ipsum() +
  theme(legend.position="none", axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8, face = 'bold')) +
  xlab("") +
  ylab("Tajima's D")

print(p)
dev.off()

file.remove(Filename1, Filename2, Filename3, Filename4)