
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

#==========================
#  NUCLEOTIDE DIVERSITY
#==========================
Filename1 <- "results/tables/sliding_windows/HPX14_Arabiensis_pi.tsv"; 
Filename2 <- "results/tables/sliding_windows/HPX14_Coluzzii_pi.tsv";
Filename3 <- "results/tables/sliding_windows/HPX14_Gambiae_pi.tsv";
Filename4 <- "results/tables/sliding_windows/HPX14_S_M_pi.tsv";

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
  write.table("results/tables/HPX14_pi.xlsx", sep = '\t',
              col.names = TRUE, row.names = FALSE, quote = FALSE)


# define color legend
cols <- c(grey(0.8), "#00FF00", "#0000FF", "#66A61E", "#E6AB02",
          "#A6761D", "#666666", RColorBrewer::brewer.pal(12, "Paired"),
          RColorBrewer::brewer.pal(11, 'Spectral'), RColorBrewer::brewer.pal(8, 'Dark2'))

mycol <- sample(cols, unique(pi_long$Populations) %>% length())

png(paste0('results/figures/HPX14_ND.png'),width = 12, height = 6, 
    family = 'serif', units = 'in', res = 600)

p <- pi_long %>%
  mutate(Populations = fct_reorder(Populations, Nuc_Div)) %>% # Reorder data
  ggplot( aes(x = Populations, y = Nuc_Div)) + # , fill=Populations, color=Populations
  geom_boxplot(lwd = 0.3, fill = mycol) +
  coord_flip() +
  theme_bw() +
  theme_ipsum() +
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8, face = 'bold')) +
  xlab("") + ylab("Nucleotide Diversity")

print(p)
dev.off()

# RIDGELINE PLOT
for (i in seq_len(nrow(pi_long))) {
    # pi_long$Species[i] <- unlist(strsplit(as.character(pi_long$Populations)[i], " ", fixed = TRUE))[2]
    pi_long$Population[i] <- unlist(strsplit(as.character(pi_long$Populations)[i], " ", fixed = TRUE))[1]
}

# basic example
pi_long %>%
    mutate(Populations = fct_rev(as_factor(Populations))) %>%
    ggplot(aes(x = Nuc_Div, y = Populations, fill = Species)) +
    # geom_density_ridges_gradient(scale = 5, size = 0.3, rel_min_height = 0.005) +
    stat_density_ridges(quantile_lines = TRUE, quantiles = 0.5) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#0000FFA0")) +
    theme_ridges() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0, vjust = 0),
          plot.title.position = "plot",
          panel.spacing = unit(0.1, "lines"),
          axis.text = element_text(size = 8, color = 'black', face = "bold.italic"),
          axis.title = element_text(size = 7, color = 'black', face = "bold.italic"))


file.remove(Filename1, Filename2, Filename3, Filename4)

#==============
#  TAJIMA'S D
#==============

Filename1 <- "results/tables/sliding_windows/HPX14_Arabiensis_tajima.tsv"; 
Filename2 <- "results/tables/sliding_windows/HPX14_Coluzzii_tajima.tsv";
Filename3 <- "results/tables/sliding_windows/HPX14_Gambiae_tajima.tsv";
Filename4 <- "results/tables/sliding_windows/HPX14_S_M_tajima.tsv"

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
  write.table("results/tables/HPX14_tajima.xlsx", sep = '\t',
              col.names = TRUE, row.names = FALSE, quote = FALSE)


mycol <- sample(cols, unique(pi_long$Populations) %>% length())

png(paste0('results/figures/HPX14_tajima.png'), width = 12, height = 6, 
    family = 'serif', units = 'in', res = 600)

p <- pi_long %>%
  mutate(Populations = fct_reorder(Populations, Tajima)) %>% # Reorder data
  drop_na() %>% 
  ggplot( aes(x = Populations, y = Tajima)) +
  geom_boxplot(lwd = 0.3, fill = mycol) +
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