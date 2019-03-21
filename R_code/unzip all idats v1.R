rm(list=ls())
library(plyr)
library(GEOquery)

home=TRUE
if(home)
{
  setwd("C:/Users/cavin/Desktop/methyl_course/GSE99788_Data/EPIC/")
  
} else
{
  setwd("M:/Meetings Presentations and Travel/SIU Meeting/Data Analysis Class/EPIC Data/Fibroblasts_GSE99788/tar_zip")
  
}

zipF <- list.files(pattern = "*idat.gz", full.names = TRUE)

ldply(.data=zipF, .fun=gunzip)


