
# https://yulab-smu.top/treedata-book/index.html
# https://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#installation-of-the-package-snprelate
# https://rpubs.com/adel922/560260
rm(list = ls())

devtools::install_version('rvcheck',version='0.1.8')
BiocManager::install("ggtree")

## Loading the following R libraries
library(gdsfmt)
library(SNPRelate) 
library(ggplot2)
library(ggtree)
library(ape)


### Save the path to the GATK vcfs to  variables  /media/Data/Data/Documents_Karim/Fadel/Aoua/My_vcf/Bamako.filtered.vcf.gz

gene.genotypes <- "data/HPX15_maf_filtered.vcf"
### turn the VCF file into a less data intensive form (GDS) for easier computing

snpgdsVCF2GDS(gene.genotypes,"HPX15_genotype.gds",method ="biallelic.only")
### Preparing the data so it is formatted correctly to create a dissimilarity matrix.

gene_genofile <- snpgdsOpen("HPX15_genotype.gds")

set.seed(100) ## making the code reproducible

ibs_gene <- snpgdsHCluster(snpgdsIBS(gene_genofile,
                                     num.thread=2, 
                                     autosome.only=FALSE,
                                     remove.monosnp=FALSE))

rvGene <- snpgdsCutTree(ibs_gene)

### Saving the dendrograms to new Variables

treeGene = rvGene$dendrogram
plot(rvGene$dendrogram,horiz=T, main ="HPX15 SNP Tree" )

tree2 = ggtree(as.phylo(as.hclust(treeGene)), layout="circular",color='darkgreen', branch.length="branch.length") + 
    geom_tiplab(size=2.5, aes(angle=angle)) + 
    ggtitle("CEU.exon.2010_03.genotypes.vcf SNP Tree")
tree2


showfile.gds(closeall=TRUE)
