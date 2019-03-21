###### script for analyzing methylation data

rm(list=ls())

setwd("C:/Users/cavin/Desktop/methyl_course/")
library(minfi)
load("GSE99788_Data/Processed Data/WB.noob.RData") # phenotype data
dim(WB.noob)
cellprop<-read.csv("GSE99788_Data/GSE99788 Pheno/GSE99788_cellprops.csv") # cell type composition
load("GSE99788_Data/Processed Data/betas.rcp.RData") # processed betas
load("GSE99788_Data/Processed Data/Gbeta.RData") # annotation file

#' load packages
suppressPackageStartupMessages({
  library(CpGassoc) # for running association analysis between methylation levels values and phenotype of interest
  library(data.table) # for fast aggregation of large data 
  library(qqman) # for visualization of data
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) # for annotation for Illumina's EPIC methylation arrays
  library(bumphunter) # for regional analysis (try to get working for ChAMP)
  library(MASS) # for basic statistics
  library(sandwich) # for linear regression (robust sandwich variance estimator)
  library(lmtest) # for testing Linear Regression Models
  library(stringi) # string manipulation
  library(knitr) # prints prettily
})

source("R_code/rmSNPandCH_DMRcate.R")

#### set up phenotypes
pheno = data.frame(pData(WB.noob))
pheno = pheno[,c("geo_accession","disease_state","gender","age","chip")]
pheno = merge(pheno, cellprop, by.x="geo_accession", by.y="samp")
pheno$disease_state <- factor(pheno$disease_state)
pheno$chip <- factor(pheno$chip)
pheno$gender <- factor(pheno$gender)
rownames(pheno) <- pheno$geo_accession
pheno <- pheno[,-1] #### remove geo_accession column, is now the rownames

rm(WB.noob)

#### quick look at the balance of disease_sate across the chips
table(pheno[,c("chip","disease_state")])

#### remove snp and non CpG probes
betas.clean = rmSNPandCH(betas.rcp,  mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY= TRUE)
nCpG = dim(betas.clean)[1]
nCpG

##### first let's look at one random snp
CpG.name = "cg05575921"
CpG.level <- betas.clean[CpG.name,]

knitr::kable(cbind(Min   = round( tapply(CpG.level,pheno$disease_state,min   ),3),
                   Mean  = round( tapply(CpG.level,pheno$disease_state,mean  ),3), 
                   Median= round( tapply(CpG.level,pheno$disease_state,median),3),
                   Max   = round( tapply(CpG.level,pheno$disease_state,max   ),3),
                   SD    = round( tapply(CpG.level,pheno$disease_state,sd    ),3),
                   N     = table( pheno$disease_state )))

#' boxplot by disease state
par(mfrow=c(1,1))
boxplot(CpG.level ~ pheno$disease_state, main=paste0("Beta-values\n", CpG.name), col=c("blue","red"),ylim=c(.5,1))

#' linear regression on betas
summary(lm(CpG.level~pheno$disease_state))$coefficients[2,c("Estimate", "Pr(>|t|)","Std. Error")]

#### now take a look at M-values
CpG.mlevel = log2(CpG.level/(1-CpG.level))

knitr::kable(cbind(Min    = round( tapply(CpG.mlevel, pheno$disease_state,min   ),3),
                   Mean   = round( tapply(CpG.mlevel, pheno$disease_state,mean  ),3), 
                   Median = round( tapply(CpG.mlevel, pheno$disease_state,median),3),
                   Max    = round( tapply(CpG.mlevel, pheno$disease_state,max   ),3),
                   SD     = round( tapply(CpG.mlevel, pheno$disease_state,sd    ),3),
                   N      = table(pheno$disease_state)))

par(mfrow=c(1,2))
boxplot(CpG.level  ~ pheno$disease_state, main=paste0("Beta-values\n",CpG.name), col=c("blue","red"))
boxplot(CpG.mlevel ~ pheno$disease_state, main=paste0("M-values\n"   ,CpG.name), col=c("blue","red"))

### now do regression on M-values
summary(lm(CpG.mlevel~pheno$disease_state))$coefficients[2,c("Estimate", "Pr(>|t|)","Std. Error")]

#' we can always extract measures of the relative quality of statistical models - e.g. adjusted R2 - to look at model performance  
#' model on betas
summary(lm(CpG.level~pheno$disease_state))$adj.r.squared

#' model on mvalues
summary(lm(CpG.mlevel~pheno$disease_state))$adj.r.squared

#'## EWAS and results using CpGassoc
#'see [Barfield et al. Bioinformatics 2012](http://www.ncbi.nlm.nih.gov/pubmed/22451269)  

#' Smoking as predictor  
#' note that CpGassoc is quite fast for running almost a million regressions!

pheno$Chrons = ifelse(pheno$disease_state=="Crohn's Disease",1,0)
system.time(results.basic <- cpg.assoc(betas.clean, pheno$Chrons)) ### very very fast

#' Bonferroni significant hits
table(results.basic$results[,3] < 0.05/(nCpG))
#' FDR significant hits
table(results.basic$results[,5] < 0.05)

#### now adjusted
results.adj = cpg.assoc(
  betas.clean
  ,pheno$Chrons
  ,covariates=pheno[,c("gender","age","chip","CD8T","CD4T","NK","Bcell","Mono","Neu")]
)

print(results.adj) ### now there are no significant results

###################################################################
#### compare bonferonni CpGs before and after adjustment


#' Lambda - this is a summary measure of genomic inflation  
#' ratio of observed vs expected median p-value - is there early departure of the qqline
#' estimated at -log10(median=0.5) ~ 0.3 on the x-axis of a qqplot  
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)

par(mfrow=c(1,2))
plot(results.basic)
plot(results.adj)

#' Lambda before cell type adjustment
lambda(results.basic$results[,3])
#' Lambda after cell type adjustment
lambda(results.adj$results[,3])

####################################################################
##### map results to genome
#' Map the results to the epigenetic annotation
IlluminaAnnot<-as.data.frame(getAnnotation(Gbeta))

#' Restrict to good quality probes and order data frames
IlluminaAnnot <- IlluminaAnnot[IlluminaAnnot$Name %in% results.adj$results$CPG.Labels,]
#### 791 probes in results but not annotation, check as to why. possibly do double filtering
IlluminaAnnot <- IlluminaAnnot[match(results.adj$results$CPG.Labels, IlluminaAnnot$Name),]

#' Check that CpGs are align
identical(IlluminaAnnot$Name,results.adj$results$CPG.Labels)

#### switch this to using a merge if CpGs don't align
datamanhat <- data.frame(CpG=results.adj$results[,1],Chr=as.character(IlluminaAnnot$chr),
                         Mapinfo=IlluminaAnnot$pos, UCSC_RefGene_Name=IlluminaAnnot$UCSC_RefGene_Name, 
                         Pval=results.adj$results[,3], Eff.Size = results.adj$coefficients[,4], Std.Error = results.adj$coefficients[,5])

#' see where the top hits are
head(datamanhat[order(datamanhat$Pval), ],n=7)

#' Volcano Plot-results2
#' with Bonferroni threshold and current FDR
plot(results2$coefficients[,4],-log10(results2$results[,3]), 
     xlab="Estimate", ylab="-log10(Pvalue)", main="Volcano Plot\nadjusted for cell proportions",ylim=c(0,8))
#Bonferroni threshold & FDR threshold
abline(h = -log10(0.05/(nCpG)), lty=1, col="red", lwd=2)

#'## Manhattan plot for cell-type adjusted EWAS  
#' Reformat the variable Chr (so we can simplify and use a numeric x-axis)
datamanhat$Chr <- as.numeric(sub("chr","",datamanhat$Chr))

#' the function manhattan needs data.frame including CpG, Chr, MapInfo and Pvalues
#' ##### fix this
manhattan(datamanhat,"Chr","Mapinfo", "Pval", "CpG", 
          genomewideline = -log10(0.05/(nCpG)), suggestiveline = FALSE,
          main = "Manhattan Plot \n adjusted for cell proportions",ylim=c(0,8))

