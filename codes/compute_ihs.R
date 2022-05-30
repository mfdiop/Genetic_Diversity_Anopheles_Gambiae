
library(rehh)

### create haps files per chromosome
hpx14 <- data2haplohh(hap_file="../Fatima/data/references/HPX14_maf_filtered.vcf.gz", 
                                     min_perc_geno.mrk = 50,min_perc_geno.hap = 50, polarize_vcf = F)
###run scan
hpx14_eehscan <- scan_hh(hpx14, lower_ehh_y_bound = 0, phased = F, discard_integration_at_border = F,
                         polarized = TRUE)


####generate ihs
ihs_hpx14 <- ihh2ihs(hpx14_eehscan, standardize = T, min_maf = 0.05, freqbin = 0.15)

###extract ihs table
ihs_hpx14_dat <- ihs_hpx14$ihs

###create df for cmplot dfs

ggplot(ihs_hpx14_dat) +
    # geom_line(aes(x=POSITION, y=-log10(LOGPVALUE)))
    geom_line(aes(x=POSITION, y=IHS))

# Extended Haplotype Homozygosity (EHH) at a given SNP
# res <- calc_ehh(hpx14, mrk = "3L", include_nhaplo = TRUE)
