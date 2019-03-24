###### script for analyzing methylation data
#' data is from: https://www.ncbi.nlm.nih.gov/pubmed/30589872

rm(list=ls())

home=TRUE

if(home)
{
  setwd("C:/Users/cavin/Desktop/methyl_course/")
} else
{
  setwd("M:/Methyl Course/methyl_course/")
}

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
  library(bumphunter) # for regional analysis 
  library(DMRcate)
  library(MASS) # for basic statistics
  library(sandwich) # for linear regression (robust sandwich variance estimator)
  library(lmtest) # for testing Linear Regression Models
  library(stringi) # string manipulation
  library(knitr) # prints prettily
})

### extracted code from DMRcate as its install doesn't work for all R versions
#source("R_code/rmSNPandCH_DMRcate.R")

#### set up phenotypes
pheno = data.frame(pData(WB.noob))
pheno = pheno[,c("geo_accession","disease_state","gender","age","chip")]
pheno = merge(pheno, cellprop, by.x="geo_accession", by.y="samp")
pheno$disease_state <- factor(pheno$disease_state)
pheno$chip <- factor(pheno$chip)
pheno$gender <- factor(pheno$gender)
rownames(pheno) <- pheno$geo_accession
pheno <- pheno[,-1] #### remove geo_accession column, is now the rownames

### can remove WB.noob now as don't need phenotype info
rm(WB.noob)

#### quick look at the balance of disease_sate across the chips
table(pheno[,c("chip","disease_state")])

#### remove snp and non CpG probes
betas.clean = rmSNPandCH(betas.rcp,  mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY= TRUE)
nCpG = dim(betas.clean)[1]
nCpG

##### first let's look at one CpG
CpG.name = "cg09234453"
CpG.level <- betas.clean[CpG.name,]

knitr::kable(cbind(Min   = round( tapply(CpG.level,pheno$disease_state,min   ),3),
                   Mean  = round( tapply(CpG.level,pheno$disease_state,mean  ),3), 
                   Median= round( tapply(CpG.level,pheno$disease_state,median),3),
                   Max   = round( tapply(CpG.level,pheno$disease_state,max   ),3),
                   SD    = round( tapply(CpG.level,pheno$disease_state,sd    ),3),
                   N     = table( pheno$disease_state )))

#' boxplot by disease state
par(mfrow=c(1,1))
boxplot(CpG.level ~ pheno$disease_state, main=paste0("Beta-values\n", CpG.name), col=c("blue","red"))

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

#' Disease State (chron's disease or not) as predictor  
#' note that CpGassoc is quite fast for running almost a million regressions!

pheno$Chrons = ifelse(pheno$disease_state=="Crohn's Disease",1,0)
system.time(results.basic <- cpg.assoc(betas.clean, pheno$Chrons, covariates=pheno[,c("chip")])) ### very very fast

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

#' FDR significant hits
table(results.adj$results[,5] < 0.05)

print(results.adj) ### now there are no significant results

###################################################################
#### QQ Plots and Volcano plots are two common ways of visualizing the data to check that basic assumptions hold
#' Volcano Plot-results.adj and QQ plot
#' Lambda - this is a summary measure of genomic inflation  
#' ratio of observed vs expected median p-value - is there early departure of the qqline
#' estimated at -log10(median=0.5) ~ 0.3 on the x-axis of a qqplot  
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)

#' with Bonferroni threshold and current FDR
par(mfrow=c(1,2))
plot(results.adj$coefficients[,4],-log10(results.adj$results[,3]), 
     xlab="Estimate", ylab="-log10(Pvalue)", main="Volcano Plot\nadjusted for cell proportions",ylim=c(0,8))
#Bonferroni threshold & FDR threshold
abline(h = -log10(0.05/(nCpG)), lty=1, col="red", lwd=2)

plot(results.adj)

#' Lambda before adjustments
lambda(results.basic$results[,3])
#' Lambda after adjustments
lambda(results.adj$results[,3])

basic.merge <- merge(results.basic$results, results.basic$coefficients, by.x="CPG.Labels", by.y=0)
full.merge <- merge(results.adj$results, results.adj$coefficients, by.x="CPG.Labels", by.y=0)

par(mfrow=c(1,1))
plot(basic.merge$effect.size[basic.merge$P.value < 0.01], full.merge$effect.size[basic.merge$P.value < 0.01], 
     xlab="Full Estimate", ylab="Basic Estimate", main="Comparison of Beta for Full and Basic Model")
#Bonferroni threshold & FDR threshold
abline(a=0, b=1, col="red", lty="dashed")
points(basic.merge$effect.size[basic.merge$P.value < 1E-6],full.merge$effect.size[basic.merge$P.value < 1E-6], 
       col="red", pch=19, cex=(1+(-log10(full.merge$P.value[basic.merge$P.value < 1E-6])/5) ) )



####################################################################
##### map results to genome
#' Map the results to the epigenetic annotation
IlluminaAnnot<-as.data.frame(getAnnotation(Gbeta))

#' annotate results
results.anno <- results.adj$results

#' check that results and effect estimates in same order
identical(results.anno$CPG.Labels, rownames(results.adj$coefficients))

results.anno <- cbind(results.anno, results.adj$coefficients)

#' Restrict to good quality probes and order data frames
IlluminaAnnot <- IlluminaAnnot[IlluminaAnnot$Name %in% results.anno$CPG.Labels,]
#### 791 probes in results but not annotation
results.anno <- results.anno[results.anno$CPG.Labels %in% IlluminaAnnot$Name,]
IlluminaAnnot <- IlluminaAnnot[match(results.anno$CPG.Labels, IlluminaAnnot$Name),]

#' Check that CpGs are align
identical(IlluminaAnnot$Name,results.anno$CPG.Labels)

#### switch this to using a merge if CpGs don't align
datamanhat <- data.frame(CpG=results.anno$CPG.Labels, Chr=IlluminaAnnot$chr,
                         Mapinfo=IlluminaAnnot$pos, UCSC_RefGene_Name=IlluminaAnnot$UCSC_RefGene_Name, 
                         Pval=results.anno$P.value, Eff.Size = results.anno$effect.size, Std.Error = results.anno$std.error)

#' see where the top hits are
head(datamanhat[order(datamanhat$Pval), ],n=7)

#'## Manhattan plot for cell-type adjusted EWAS  
#' Reformat the variable Chr (so we can simplify and use a numeric x-axis)
datamanhat$Chr <- as.numeric(sub("chr","",datamanhat$Chr))

#' the function manhattan needs data.frame including CpG, Chr, MapInfo and Pvalues
#' ##### fix this

par(mfrow=c(1,1))
qqman::manhattan(datamanhat,"Chr","Mapinfo", "Pval", "CpG", 
          genomewideline = -log10(0.05/(nCpG)), suggestiveline = FALSE,
          main = "Manhattan Plot \n Adjusted Model",ylim=c(0,8))


###### Regional analyses
#### Regional analyses can be more powerful than individual CpG analyses as they aggregate signals from a region
#' First we need to set up a model

model = model.matrix( ~Chrons+factor(chip),data=pheno)

#'Let's run the regional analysis using the Beta-values from our preprocessed data
#' First annotate the data so that the regions can be determined based on nearby probes
myannotation <- cpg.annotate("array", betas.clean, analysis.type="differential",arraytype="EPIC",
                             what="Beta",design=model, coef=2)

#'Regions are now agglomerated from groups of significant probes 
#'where the distance to the next consecutive probe is less than lambda nucleotides away
dmr.chrons <- dmrcate(myannotation, lambda=1000, C=2)

#'Let's look at the results
head(dmr.chrons$results)

#'Visualizing the data can help us understand where the region lies 
#'relative to promoters, CpGs islands or enhancers

#' Let's extract the genomic ranges and annotate to the genome
results.ranges <- extractRanges(dmr.chrons, genome = "hg19")

#' if you are interested in plotting genomic data the Gviz is extremely useful
#'Let's look at the first region
results.ranges[1]

# set up the grouping variables and colours
cols = c("magenta","red")[pheno$Chrons]
names(cols) = levels(pheno$Chrons)[pheno$Chrons]

#'Draw the plot for the top DMR\
#+ fig.width=8, fig.height=6, dpi=300
DMR.plot(ranges=results.ranges, dmr=1, CpGs=betas.clean, phen.col=cols, what = "Beta",
         arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
         genome="hg19", samps=1:nrow(pheno))

#' cleanup
rm(tx.hg19,tx.hg38,tx.mm10,snpsall,myBetas,myannotation,crosshyb,XY.probes);gc()


#'Extracting CpGs-names and locations
coord = dmr.chrons$results$coord[1]
coord = stri_match(coord,regex="^(chr.+):(\\d+)-(\\d+)$")

chr = coord[2]
start = as.integer(coord[3])
end = as.integer(coord[4])

#'CpG ID and individual metrics

cpgs = subset(dmr.chrons$input, CHR == chr & pos >= start & pos <= end)
knitr::kable(cpgs)

##### enrichement analysis with MissMethyl
library(missMethyl)

gst <- gometh(sig.cpg=basic.merge$CPG.Labels[basic.merge$FDR<0.1], all.cpg=basic.merge$CPG.Labels, collection="GO", array.type="EPIC", prior.prob=TRUE)
gst <- subset(gst, N > 2) ### require 3 or more CpGs
head(gst[order(gst$P.DE),])

## Mendelian Randomization
#' examine a Mendelian Randomization analysis 

library(TwoSampleMR)
library(MRInstruments)

data("aries_mqtl")

adult_mqtl_basic <- subset(aries_mqtl, cpg%in%basic.merge$CPG.Labels[basic.merge$P.value < 1E-4])
aries_exp_basic <- format_aries_mqtl(subset(adult_mqtl_basic, age=="Middle age"))
basic <- subset(basic, CPG.Labels%in%adult_mqtl_basic$cpg)
### clump SNPs that are in LD
aries_exp_basic <- clump_data(aries_exp_basic)

ao <- available_outcomes()

chrons <- subset(ao, grepl("Crohn's", ao$trait) & mr==1)
chrons <- chrons[order(chrons$sample_size, decreasing=TRUE),]

### extract relevant instruments
#' extracting from [https://www.ncbi.nlm.nih.gov/pubmed/26192919]
chrons_instruments <- extract_outcome_data(outcomes=12, snps=unique(c(aries_exp_basic$SNP, aries_exp_full$SNP)))

### if you have trouble with the above scripts can load the pre-extracted instruments and exposures
# MR Aries mqtl exposure.RData
# MR Chrons Instruments.RData

### clump SNPs that are in LD
dat_basic <- harmonise_data(exposure_dat = aries_exp_basic,
                            outcome_dat = chrons_instruments)
dat_basic <- power.prune(dat_basic,method.size=T)

mr_basic <- mr(dat_basic)
mr_basic$exposure <- gsub(" \\(Middle age\\)","",mr_basic$exposure)

mr_basic <- merge(mr_basic, basic.merge[,c("CPG.Labels","P.value","effect.size")], by.x="exposure", by.y="CPG.Labels")
mr_basic <- merge(mr_basic, IlluminaAnnot[,c("Name","chr","pos","UCSC_RefGene_Name","Relation_to_Island")], by.x="exposure", by.y="Name")
mr_basic <- mr_basic[order(mr_basic$pval),]


mr_basic[mr_basic$pval < 0.05,]

### check memeory usage
pryr::mem_used()