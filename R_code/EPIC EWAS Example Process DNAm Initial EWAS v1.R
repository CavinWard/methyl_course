### code to process DNA methylation data

rm(list=ls())

home=TRUE
test = FALSE
if(home)
{
  if(test)
  {
    path <- "EPIC/test"
  } else
  {
    path <- "EPIC"
  }
  
  setwd("C:/Users/cavin/Desktop/methyl_course/GSE99788_Data/")
} else
{
  path <- "EPIC Data/Fibroblasts_GSE99788/tar_zip/test"
  setwd("M:/Meetings Presentations And Travel/SIU Meeting/Data Analysis Class/")
}


suppressMessages(library(minfi)) # popular package for methylation data
library(shinyMethyl) # for visualizing quality control
library(pryr) # for monitoring memory use
suppressMessages(library(matrixStats)) # for calculating summary statistics
library(ENmix) # probe type adjustment "rcp"
suppressMessages(library(limma)) # for MDS plots
suppressMessages(library(reshape,scales)) # reshape data and graphig 
suppressMessages(require(sva)) # for addressing batch effects
library(IlluminaHumanMethylationEPICmanifest) # manifest for Illumina's EPIC methylation 
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) # annotation for Illumina's EPIC methylation arrays.
library(FlowSorted.Blood.EPIC)
library(GEOquery)


gse=getGEO(filename="GSE99788 Pheno/GSE99788_series_matrix.txt.gz", getGPL=FALSE)

pheno <- gse@phenoData@data
rm(gse); gc()


cols_to_keep <- c("geo_accession","gender:ch1","disease state:ch1","disease stage:ch1","montreal score:ch1","passage:ch1",
                  "sampleID:ch1","thiopurine usage:ch1","description","characteristics_ch1.4")

pheno <- pheno[,cols_to_keep]
names(pheno)[which(names(pheno)=="characteristics_ch1.4")] <- "age"
pheno$age <- as.integer(gsub("age: ","",pheno$age))

colnames(pheno) <- gsub(":ch1","",colnames(pheno))
colnames(pheno) <- gsub(" ","_",colnames(pheno))

#### have to set the base name of the idat files
pheno$Basename <- paste0(pheno$geo_accession, "_", pheno$description)

### trick to force factor to character
pheno$description <- paste0(pheno$description)

### call str to make sure everying is what you expect including the variable types
str(pheno)

#### read in the meth data
WB <- read.metharray.exp(base=path, targets=NULL, verbose=T) # read the idat file one by one, specifying targets=pheno is safer as it forces you to specify with 'Basename' what all the files are 
ncol(WB)

### set phenotype data
pData(WB) <- as(subset(pheno, Basename%in%WB@colData@rownames),"DataFrame") ### sometimes have to change 'Basename' to 'geo_accession'

### visusalize QC with shiny
summaryqc <- shinySummarize(WB)
# runShinyMethyl(summaryqc)
# rm(covariates, colorSet, current.control.type, current.density.type, current.probe.type, 
#   genderCutoff, mouse.click.indices, sampleColors, summaryqc)


### can take a lot of RAM
cellprop <- estimateCellCounts2(WB, compositeCellType = "Blood", referencePlatform = "IlluminaHumanMethylationEPIC")

cell_table <- cellprop[[1]]
cell_table <- cbind(rownames(cell_table), cell_table)
colnames(cell_table)[1] <- "samp"

write.table(cell_table, file="GSE99788 Pheno/GSE99788_cellprops.csv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")

rm(FlowSorted.Blood.EPIC);gc()

#### background correct using noob
WB.noob <- preprocessNoob(WB)

#### comparing betas before and after preprocessing with noob
densityPlot(WB, main = "density plots before and after preprocessing", pal="#440154FF", ylim=c(0,4.5))
densityPlot(WB.noob, add = F, pal = "#FDE725FF")
# Add legend
legend("topleft", c("Noob","Raw"), 
       lty=c(1,1), title="Normalization", 
       bty='n', cex=1.3, col=c("#FDE725FF","#440154FF"))

#### detection p-values give a p-value associated with the deviation of a probe's inensity from background
#### only want to keep those with low p-value (high intensity) to avoid contamination with noise
detect.p <- minfi::detectionP(WB, type = "m+u")
#' Let's look at the median detection P-values
#+ fig.width=8, fig.height=6, dpi=300
barplot(colMeans(detect.p), col=rainbow(dim(detect.p)[2]), las=2, 
        cex.names=0.7, main="Mean detection P by sample",cex.axis=0.8, ylim=c(0,3e-4))

#### we will restrict to probes with detection p-value <= 0.01
sum(detect.p > 0.01) ### < 20,000 probes
#'Restrict data to good probes only:
detect.p[detect.p > 0.01] <- NA
detect.p <- na.omit(detect.p)
intersect <- intersect(rownames(getAnnotation(WB)), rownames(detect.p))
length(intersect)
#' Filter bad probes from our methylset
nrow(WB.noob)
WB.noob <- WB.noob[rownames(getAnnotation(WB.noob)) %in% intersect,]
nrow(WB.noob)
# cleanup
rm(intersect, detect.p);gc()

### probe type adjustment adjusts for the Type I vs Type II probes
#' RCP with EnMix: Regression on Correlated Probes [Niu et al. Bioinformatics 2016](http://www.ncbi.nlm.nih.gov/pubmed/27153672)
betas.rcp <- rcp(WB.noob)
dim(betas.rcp)
#' note that this package takes beta values out of the minfi object - result is a matrix
class(betas.rcp) 


## Annotation of Infinium type for each probe (I vs II)
typeI <-   minfi::getProbeInfo(WB.noob,type="I")$Name
typeII <-  minfi::getProbeInfo(WB.noob,type="II")$Name
onetwo <- rep(1, nrow(betas.rcp))
onetwo[rownames(betas.rcp) %in% typeII] <- 2
# vast majority of probes are type II
#knitr::kable(t(table(onetwo)))
t(table(onetwo))

#' Density plots by Infinium type: before and after RCP calibration  
#' Probe-type bias adjustment before and after RCP
#+ fig.width=15, fig.height=7, dpi=300
par(mfrow=c(1,2)) # Side-by-side density distributions 
densityPlot(WB.noob[rownames(getAnnotation(WB.noob)) %in% typeI,],pal = "#FDE725FF",main='Beta density',ylim=c(0,6.5))
densityPlot(WB.noob[rownames(getAnnotation(WB.noob)) %in% typeII,],add = F, pal = "#440154FF")
densityPlot(betas.rcp[rownames(getAnnotation(WB.noob)) %in% typeI,],pal = "#FDE725FF",main='Beta density probe-type adjusted',ylim=c(0,6.5))
densityPlot(betas.rcp[rownames(getAnnotation(WB.noob)) %in% typeII,],add = F, pal = "#440154FF")
legend("center", c("Infinium I","Infinium II"), 
       lty=c(1,1), title="Infinium type", 
       bty='n',col=c("#FDE725FF","#440154FF"))
#' notice that the type I and II peaks are more closely aligned after rcp adjustment  
#' (particularly in the higher peak)
rm(onetwo, typeI, typeII)

##### add chip and position to pdata
pData(WB.noob)$chip <- sapply(strsplit(pData(WB.noob)$description, "_"), "[", 1)
pData(WB.noob)$pos <- sapply(strsplit(pData(WB.noob)$description, "_"), "[", 2)

#### calculate some principal components to see major sources of variability
PCobject = prcomp(t(betas.rcp), retx = T, center = T, scale. = T)

#' Extract the Principal Components from SVD
PCs <- PCobject$x
#' Proportion of variance explained by each additional PC
cummvar <- summary(PCobject)$importance["Cumulative Proportion", 1:10]
t(as.matrix(cummvar))

#### is PC1 associatd with the chip used and the position on the chip?
par(mfrow = c(1, 2))
boxplot(PCs[, 1] ~ pData(WB.noob)$chip,
        ylab = "PC1",las=2, main="Chip Used",col=rainbow(8))
boxplot(PCs[, 1] ~ pData(WB.noob)$pos,
        ylab = "PC1",las=2, main="Position on chip",col=rainbow(8))

#### combination of both
par(mfrow = c(1, 1))
boxplot(PCs[, 1] ~ pData(WB.noob)$description,
        ylab = "PC1",las=2, main="Chip Used")

summary(lm(PCs[, 1] ~ pData(WB.noob)$chip))
summary(lm(PCs[, 1] ~ pData(WB.noob)$pos))

### use combat to remove position batch effects and adjust for chip in analyses 

Mvals <- log2(betas.rcp)-log2(1-betas.rcp)
#' ComBat eBayes adjustment using a known variable of interest (here we use row)
Mvals.ComBat <- ComBat(Mvals, batch = pData(WB.noob)$pos)
# Convert M-values back to beta-values
betas.rcp <- 2^Mvals.ComBat/(1+2^Mvals.ComBat)

PCobject <- prcomp(t(betas.rcp), retx = T, center = T, scale. = T)
PCs <- PCobject$x
#' The not difference in first PC across the positions
boxplot(PCs[, 1] ~ pData(WB.noob)$pos,
        ylab = "PC1",las=2, main="Position on chip",col=rainbow(8))

#### still a difference by chip
boxplot(PCs[, 1] ~ pData(WB.noob)$chip,
        ylab = "PC1",las=2, main="chip used",col=rainbow(8))

summary(lm(PCs[, 1] ~ pData(WB.noob)$pos))

#cleanup
rm(PCs, Mvals, cummvar, PCobject)

#'## predict sex from methylation
Gbeta <- mapToGenome(WB)  #map to the genome
#' getSex predicts sex based on X and Y chromosome methylation intensity
getSex(Gbeta) 
#' we see that our predictions match the phenodata
table(pData(WB)$gender,getSex(Gbeta)$predictedSex)

#### save out some data so that you don't have to recreate it
save(WB,WB.noob,file="Processed Data/WB.noob.RData")
save(betas.rcp,file="Processed Data/betas.rcp.RData")
save(Gbeta,file="Processed Data/Gbeta.RData")

pryr::mem_used()
