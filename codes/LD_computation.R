
# https://pastebin.com/FMdqpiEJ
# https://www.cog-genomics.org/plink/2.0/ld
# https://wlz0726.github.io/2017/05/27/LD-prune-with-plink/
# https://avikarn.com/2019-07-30-prunning/
# https://www.biostars.org/p/300381/

# Clear workspace
rm(list = ls())

# Load library
library(tidyverse)

plink <- "/home/karim/Documents/Mes_Programmes/plink_linux_x86_64_20181202/plink"

## Set Variables
vcf_file <-  "../Fatima/data/references/HPX14_maf_filtered.vcf.gz"
Output_dir <- "../Fatima/results/LD"
Gene <- gsub("_maf_filtered.vcf.gz", "", basename(vcf_file))

# Output_dir <- paste0(Output_dir, Gene)
if(!dir.exists(Output_dir)) dir.create(Output_dir)

### Load Metadata file
metadata  <- readxl::read_xlsx("../Fatima/data/metadata/global_Metadata.xlsx")

## Arrange the list of countries by name
num_isolates <- metadata %>% group_by(species) %>% count(countries)

Species <- metadata %>% 
    select(species) %>% 
    unique() %>% 
    pull()

mind <- 0.1; geno <- 0.1; hwe <- 0.0000001
maf <- 0.05; ld_window <- 100; ld_window_kb <- 10

for(j in 1:length(Species))
{
    ## Extract specific population sample IDs
    targetID <- metadata %>% 
        filter(species == Species[j]) %>% 
        pull(sample_ID)
    
    pathToTargetID <- file.path(Output_dir, paste0(Gene, "_", tolower(Species[j]), '.txt'))
    
    ## Save the file
    write.table(targetID, pathToTargetID, row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    ## Extract isolates from the raw VCF file
    pathToTargetVCF <- file.path(Output_dir, paste0(Gene, "_", tolower(Species[j]), '.vcf.gz'))
    
    system(sprintf("bcftools view -S %s -o %s -O z %s", pathToTargetID, pathToTargetVCF, vcf_file))

    ########################
    ## LD within chromosome
    #######################
    Output <- file.path(Output_dir, gsub(".vcf.gz", "_afterQC", basename(pathToTargetVCF)))
    
    # quality control
    # select one chromosome
    system(paste0(plink, " --vcf ", pathToTargetVCF, " --mind ", mind, " --geno ", geno, 
                  " --maf ", maf, " --hwe ", hwe, " --make-bed --allow-extra-chr --out ", Output))
    
    # LD output file
    ld_output <- file.path(Output_dir, gsub(".vcf.gz", "", basename(pathToTargetVCF)))
    system(paste0(plink, " --bfile ", Output, " --r2 --ld-window ", ld_window," --ld-window-kb ", ld_window_kb,
                  " --ld-window-r2 0 --allow-extra-chr --out ", ld_output))
    
    # LD block length distribution,
    system(paste0(plink, " --bfile ", Output, " --blocks no-pheno-req --blocks-max-kb 5 --allow-extra-chr --out ", ld_output))
    
    system(sprintf("rm -r %s %s %s %s %s %s %s %s %s %s", pathToTargetID, pathToTargetVCF, 
                   paste0(ld_output, ".nosex"), paste0(ld_output, ".log"),
                   paste0(Output, ".nosex"), paste0(Output, ".log"),
                   paste0(Output, ".bed"),paste0(Output, ".bim"),
                   paste0(Output, ".fam"), paste0(Output, ".irem")))
}

#=====================
# LD DECAY LINE PLOT
#=====================
Path <- "../Fatima/results/LD"
files <- list.files(Path, pattern = ".ld", full.names = TRUE)

names(files) <- str_replace(string = files,
                            pattern = paste0(Path, "/HPX14_(.*).ld"),
                            replacement = "\\1")

## Combine files by adding filenames
ld <- map_dfr(.x = files, .f = read_table, .id = "Species")

ld <- ld %>%
    mutate(markerDistance = abs(BP_B - BP_A))

ld <- ld %>% 
    group_by(Species) %>%
    arrange(Species, markerDistance) %>% 
    # mutate(Loess = predict(loess(R2 ~ markerDistance, span = .5, data=., degree=1, control = loess.control(surface = "direct")))) # remove data=., 
    mutate(Loess = predict(loess(R2 ~ markerDistance, span = .5, degree=1, control = loess.control(surface = "direct"))))

# PLOT LD 
ggplot(ld) +
    geom_point(aes(x=markerDistance, y=R2, color = Species))



# # read in LD results file
# LdValues <- read_table(paste0(Output, "_1.ld"))
# 
# # calculate LD in 20 kb bins to display the trendline
# averageLD <- LdValues %>%
#     mutate(markerDistance = abs(BP_B - BP_A)) %>%
#     # dplyr::filter(markerDistance < 5000) %>%
#     mutate(intervals = cut_width(markerDistance, 100, boundary = 0)) %>%
#     group_by(intervals) %>%
#     summarise_at(vars(R2),funs(mean(., na.rm=TRUE))) %>%
#     rename(averageR2=R2)
# 
# # calculate inter marker distances
# fullLD <- LdValues %>%
#     mutate(markerDistance = abs(BP_B - BP_A)) %>%
#     # dplyr::filter(markerDistance < 5000) %>%
#     mutate(intervals = cut_width(markerDistance, 100, boundary = 0))
# 
# #merge the two data sets (full LD info and average per bin)
# mergedLD <- full_join(fullLD,averageLD, by = "intervals")
# 
# # visualize LD decay
# ggplot(mergedLD) +
#     geom_point(aes(x=markerDistance, y=R2)) +
#     geom_line(aes(x=markerDistance, y=averageR2), color="red", size=2)
# 
# 
# ########################
# ## Scatter LD Plot
# ########################
# 
# LD$distancekb <- with (LD, LD$BP_B-LD$BP_A)/1000 ## the distance between snp1 and snp2 in kb
# LD$grp <- cut(LD$distancekb, 0:5) ## bin 70kb
# r2means <- with (LD, tapply(LD$R2, LD$grp, FUN = mean)) ##r2 mean every 1kb
# 
# ggplot(LD, aes(x = distancekb,y = R2))+
#     geom_point(shape = 20, size = 1, alpha = 0.7) +
#     geom_smooth() +
#     labs(x = "Distance (Bases)", y = expression(italic(LD~(r^{2})))) +
#     theme_bw(base_size = 14) +
#     theme(panel.border = element_blank(),
#           axis.ticks = element_blank()) %>%
#     
#     return()
# 
# ------------------------------------------------------------------------------------------------
# cat HPX14_coluzzii.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > HPX14_coluzzii.ld.summary
# dfr <- read.delim(paste0(Output, ".ld.summary"),sep="",header=F,check.names=F,stringsAsFactors=F)
# colnames(dfr) <- c("dist","rsq")
# 
# dfr$distc <- cut(dfr$dist, breaks = seq(from = min(dfr$dist)-1, to = max(dfr$dist)+1, by=500))
# dfr1 <- dfr %>% group_by(distc) %>% summarise(mean = mean(rsq), median = median(rsq))
# 
# dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
#                         end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
#                         mid=start+((end-start)/2))
# 
# ggplot()+
#     geom_point(data=dfr1,aes(x=start,y=mean),size=0.4,colour="grey20")+
#     geom_line(data=dfr1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
#     labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
#     # scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8"))+
#     theme_bw()
# 
# 
# 
# 
# dfr <- read.delim("../Fatima/results/HPX14_coluzzii_block.blocks.det",sep="",header=T,check.names=F,stringsAsFactors=F)
# colnames(dfr) <- tolower(colnames(dfr))
# 
# # ld block density
# p <- ggplot(dfr,aes(x=kb))+
#     geom_density(size=0.5,colour="grey40")+
#     labs(x="LD block length (Kb)",y="Density")+
#     theme_bw()
# 
# p
# 
# ggsave("snp-thin-ld-blocks.png",p,height=8,width=8,units="cm",dpi=250)
