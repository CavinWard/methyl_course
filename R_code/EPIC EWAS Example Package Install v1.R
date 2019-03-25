#### building exmple code for epigenetics boot camp

rm(list=ls())

cran_packages <- c("dplyr","Rcpp","openssl","CpGassoc","ggplot2", "matrixStats","reshape","glmnet","statmod","XML",
                      "pryr", "data.table", "qqman", "RPMM", "MASS", "sandwich", "lmtest","foreach", "stringi","doParallel","devtools","BiocManager")

cran_install = setdiff(cran_packages, installed.packages()[,1])
cat(paste(length(cran_install),"out of",length(cran_packages), "CRAN packages need to be installed\n"))

if(length(cran_install) >= 1)
{
  install.packages(cran_install ,dependencies=TRUE)
}

### if TwoSampleMR not already installed
if(!any(grepl("TwoSampleMR", installed.packages()[,1])))
{
  if(any(grepl("devtools",installed.packages()[,1])) ) ## if devtools installed and TwoSampleMR still not
  {
    library(devtools)
    devtools::install_github("MRCIEU/TwoSampleMR")
    devtools::install_github("MRCIEU/MRInstruments")
  } else
  {
    cat("devtools not installed; cannot install TwoSampleMR")
  }
}

cat("CRAN package install attempt done\n")

cran_install = setdiff(cran_packages, installed.packages()[,1])
if(length(cran_install) > 0)
{
  cat(paste(length(cran_install),"out of",length(cran_packages), "still not installed\n"))
  cat(cran_install)
} else
{
  cat("all CRAN packages installed\n")
}

bioc_packages <- c("bumphunter","minfi", "FlowSorted.Blood.450k", "missMethyl", "ENmix","IlluminaHumanMethylation450kanno.ilmn12.hg19",
                      "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylationEPICmanifest",
                      "sva", "IlluminaHumanMethylationEPICanno.ilm10b2.hg19", 
                      "DMRcate", "shinyMethyl","wateRmelon","FDb.InfiniumMethylation.hg19", "GEOquery")

bioc_install <- setdiff(bioc_packages, installed.packages()[,1])

cat(paste(length(bioc_install),"out of",length(bioc_packages), "Biocondutor packages need to be installed\n"))

if(length(bioc_install >= 1)) 
{
  source("https://bioconductor.org/biocLite.R")
  biocLite(bioc_install, suppressUpdates = TRUE)
  
} 

### if FlowSorted.Blood.EPIC not already installed
if(!any(grepl("FlowSorted.Blood.EPIC", installed.packages()[,1])))
{
  if(any(grepl("BiocManager",installed.packages()[,1])) ) ## if BiocManager installed 
  {
    BiocManager::install("FlowSorted.Blood.EPIC", version = "3.8")
  } else
  {
    cat("BiocManager not installed; cannot install FlowSorted.Blood.EPIC")
  }
}

cat("BioConductor package install attempt done\n")
bioc_install <- setdiff(bioc_packages, installed.packages()[,1])

if(length(bioc_install) > 0)
{
  
  cat(paste(length(bioc_install),"out of",length(bioc_packages), "still not installed\n"))
  cat(bioc_install)
} else
{
  cat("All Bioconductor packages installed\n")
}

