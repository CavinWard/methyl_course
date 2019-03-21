######
## calculation of PhenoAge DNAmAge

rm(list=ls())

library(DMwR)

setwd("Y:/1. Raw Data Files/DNA methylation aging/PhenoAge")

load("Y:/1. Raw Data Files/DNA methylation_450k data/DNHS192_meta_beta_postComBat_PTSDpm_Gender.Rdata")

PhenoAge <- read.csv("PhenoAge coefficients.csv", stringsAsFactors=FALSE, fill=TRUE, header=TRUE)

PhenoAge.nomiss <- subset(PhenoAge, PhenoAge$CpG %in% rownames(reversbeta)) ### all present

h.beta <- reversbeta[PhenoAge.nomiss$CpG,]
h.beta.0impute <- h.beta
h.beta.0impute[is.na(h.beta.0impute)] <- 0 ### removes those CpGs which are missing in individuals

h.ages <- t(t(PhenoAge.nomiss$Weight) %*% h.beta) + PhenoAge[PhenoAge$CpG=="Intercept","Weight"]
h.ages.0impute <- t(t(PhenoAge.nomiss$Weight) %*% h.beta.0impute) + PhenoAge[PhenoAge$CpG=="Intercept","Weight"]
h.ages.miss <- t(t(colSums(is.na(h.beta))))

knn.beta <- knnImputation(h.beta, k=10, scale=TRUE, meth="weighAvg") ### uses a weighted averge (as opposed to median), to compute imputation based on 10 knn
rownames(knn.beta) <- rownames(h.beta)
h.ages.knn <- t(t(PhenoAge.nomiss$Weight) %*% knn.beta) + PhenoAge[PhenoAge$CpG=="Intercept","Weight"]

samp <- read.csv("../DNHS_EpigeneticAging_Sample.csv")

samp <- merge(samp, h.ages, by.x="rownames", by.y=0)

names(samp)[ncol(samp)] <- "DNAmPhenoAge"

samp <- merge(samp, h.ages.0impute, by.x="rownames", by.y=0)
names(samp)[ncol(samp)] <- "DNAmPhenoAge_0impute"

samp <- merge(samp, h.ages.knn, by.x="rownames", by.y=0)
names(samp)[ncol(samp)] <- "DNAmPhenoAge_knn_impute"

samp <- merge(samp, h.ages.miss, by.x="rownames", by.y=0)
names(samp)[ncol(samp)] <- "DNAmPhenoAge_NmissCpGs"

save(samp, file="DNAm PhenoAge Estimation.RData")

write.table(samp, file="DNAm PhenoAge Estimation.csv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")
