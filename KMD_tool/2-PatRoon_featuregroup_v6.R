#############################################
## Title: PatRoon - feature group only 
###############################################################################
## version:2.8
## Date: February 2024
## Author: Boris Droz 
## Modified from Tutorial and Handbook on https://github.com/rickhelmus/patRoon
## Depends:
##        R(>=4.3.1)
##        patRoon(>=2.3.0)
##
## Description:
## Performed peak picking
## 
## Should have a input file with a raw, mzxml folder and the sample list (csv)
## need to create an empty output folder
## 
###############################################################################
## Parameters, directories, and parameters. Change if necessary below
############
## library
##########
# load patRoon profile
# source("D:/Patroon_NTS/.RData.R")
library(patRoon) # v2.3.0
library(xcms)
library(dplyr)
library(glue)

## path
workPath <- getwd()
workdir <- glue("{workPath}/")

## Enter your sample name list here
sample.list <- "sample_list_KMD.csv"


## Optimized XCMS parameters for peak picking
opt.ppm = 17
opt.pw = c(15, 140) # peak width min and max

## Parameter for filtering check patroon  help(filter)
min.intensity.thr = 200## absMinIntensity, typical range between 100 - 1000
rp.feature = 1 #relMinReplicateAbundance 
bk.sa.thr = 3 # blankThreshold - never go under 3

## Adduct and formula search parameter
adduct <- c( "[M-H]-")

##   ----  SCRIPT STARTS HERE ---- #############################################
################################################################################ 
# #################################################################################################################################################
# FUNCTION check and produced subDir folder
###########################################
#February 2017 -- mod on the 2022-08-18 

creat.subDir <- function (mainDir,subDir)
{
  if ( dir.exists(paste(mainDir,"/",subDir, sep="") ) ){
    
    i <- 1
    while( file.exists( paste(mainDir,"/",subDir,"_",i, sep="") ) )
    {i <-i+1}
    
    dir.create(file.path(mainDir, paste(subDir,"_",i, sep="") ))
    outpath <- file.path(mainDir, paste(subDir,"_",i, sep=""))
    
  } else {
    dir.create(file.path(mainDir, subDir))
    outpath <- file.path(mainDir, subDir)
  }
  
  return(outpath)
}

############################################################################
############################################################################
# -------------------------
# initialization
# ------------------------
# Set outpath folder
# date <- Sys.Date()
# folder <- paste("/",date,"_NTA_patRoon", sep="")
# outpath <- creat.subDir(paste(workdir,"/output",sep=""), folder)
outpath <- paste(workdir,"output",sep="")
inpath <- paste(workdir,"input",sep="")

## load data info
setwd(inpath)
df <- read.csv(sample.list, sep=",",header=TRUE)

setwd(workdir) # set directory

anaInfo <- data.frame(cbind(path = df$path,
                                analysis =df$filename,
                                group = df$group,
                                blank = df$blank) )
print(anaInfo)



# save parameter of the script
f.info <- paste(outpath,"/AA_XCMS_FEATUREGROUP_PARAMETERS_LOG.txt",sep="")
cat( paste("*** patRoon parameter for the run....", Sys.Date()), 
     file= f.info, append=TRUE, sep="\n")
cat( "#########################################################", 
     file= f.info, append=TRUE,sep="\n")
cat( paste("SampleList: ", workPath,"/input/",sample.list,sep=""), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("XCMS_ppm:", opt.ppm), file= f.info, append=TRUE,sep="\n")
cat( paste("XCMS_peakwidth:", opt.pw), file= f.info, append=TRUE,sep="\n")
cat( paste("absMinIntensity:", min.intensity.thr), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("relMinReplicateAbundance :", rp.feature), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("blankThreshold:", bk.sa.thr), file= f.info, append=TRUE,sep="\n")
cat( paste("adduct:", adduct), file= f.info, append=TRUE,sep="\n")



# -------------------------
# features
# -------------------------
# to disable parallelization when required:
# BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

# Find all features
# NOTE: see the XCMS manual for many more options
## --> first need to optimized parameter ......
param.xcms <- xcms::CentWaveParam(ppm = opt.ppm,
                              peakwidth = opt.pw,
                              snthresh = 3,
                              prefilter = c(3, 100),
                              noise = 0 ) # noise could be sensitive

fListPos <- findFeatures(anaInfo, "xcms3", param = param.xcms)

fList <- makeSet(fListPos, adducts = adduct )  

# performed RT alignement and group feature
fGroups <- groupFeatures(fList, "xcms3")
                        
# Basic rule based filtering
fGroups <- patRoon::filter(fGroups,  
                            absMinIntensity = min.intensity.thr, 
                            relMinReplicateAbundance = rp.feature, 
                            maxReplicateIntRSD = 0.75,
                            blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                            retentionRange = NULL, mzRange = NULL)

# -------------------------
# reporting
# -------------------------
## export groupfeature as table

df.fGroups <- as.data.table(fGroups)

write.table(df.fGroups, file=paste(outpath,"/featureGroupsXCMS_neg_example.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE )

