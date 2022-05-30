
# https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html

# DAPC requires the adegenet package. Let's load this package:
library(adegenet)
data(H3N2) # load the H3N2 influenza data. Type ?H3N2 for more info.
pop(H3N2) <- H3N2$other$epid
dapc.H3N2 <- dapc(H3N2, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(H3N2) - 1)
scatter(dapc.H3N2, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)

## LOAD VCF
VCF <- "../data/proceed/all_HPX15.vcf.gz"

## LOAD METADATA
metadata <- readxl::read_xlsx("../data/metadata/metadata.xlsx")
hpx15.VCF <- read.vcfR(VCF)

pop <- metadata %>% select(species) %>% pull()
ind <- metadata %>% select(sample_ID) %>% pull()

gen.hpx15 <- vcfR2genind(hpx15.VCF, return.alleles = TRUE, ind.names = ind, pop = pop) 

dapc.hpx15 <- dapc(gen.hpx15, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(gen.hpx15) - 1)
scatter(dapc.hpx15, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)

# Next, let’s assess if there are alleles that most differentiate the 2006 cluster from those in other years.
set.seed(4)
contrib <- loadingplot(dapc.hpx15$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)

temp    <- seploc(gen.hpx15)       # seploc {adegenet} creates a list of individual loci.
snp1  <- tab(temp[["3L_10786285"]]) # tab {adegenet} returns a matrix of genotypes
snp2  <- tab(temp[["3L_10786185"]])

# The following two commands find the average allele frequencies per population
(freq1 <- apply(snp1, 2, function(e) tapply(e, pop(gen.hpx15), mean, na.rm = TRUE)))
