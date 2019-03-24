rm(list=ls())

setwd("C:/Users/cavin/Desktop/methyl_course/")

library(data.table)
library(TwoSampleMR)
library(MRInstruments)

basic <- fread(file="GSE99788_basic_results.csv")
full <- fread(file="GSE99788_adj_results.csv")

basic <- subset(basic, P.value < 1E-3)
full <- subset(full, P.value < 1E-3)

data("aries_mqtl")

adult_mqtl_basic <- subset(aries_mqtl, cpg%in%basic$CPG.Labels)
aries_exp_basic <- format_aries_mqtl(subset(adult_mqtl_basic, age=="Middle age"))
basic <- subset(basic, CPG.Labels%in%adult_mqtl_basic$cpg)
### clump SNPs that are in LD
aries_exp_basic <- clump_data(aries_exp_basic)

adult_mqtl_full <- subset(aries_mqtl, cpg%in%full$CPG.Labels)
aries_exp_full <- format_aries_mqtl(subset(adult_mqtl_full, age=="Middle age"))
full <- subset(full, CPG.Labels%in%adult_mqtl_full$cpg)

ao <- available_outcomes()

chrons <- subset(ao, grepl("Crohn's", ao$trait) & mr==1)
chrons <- chrons[order(chrons$sample_size, decreasing=TRUE),]

### extract relevant instruments
chrons_instruments <- extract_outcome_data(outcomes=chrons$id, snps=c(aries_exp_basic$SNP, aries_exp_full$SNP))

### clump SNPs that are in LD
dat_basic <- harmonise_data(exposure_dat = aries_exp_basic,
                            outcome_dat = chrons_instruments)
dat_basic <- power.prune(dat_basic,method.size=T)

dat_full <- harmonise_data(exposure_dat = aries_exp_full,
                            outcome_dat = chrons_instruments)
dat_full <- power.prune(dat_full,method.size=T)

mr_basic <- mr(dat_basic)
mr_basic <- mr_basic[order(mr_basic$pval),]

mr_full <- mr(dat_full)
mr_full <- mr_full[order(mr_full$pval),]

mr_basic[mr_basic$pval < 0.05,]
