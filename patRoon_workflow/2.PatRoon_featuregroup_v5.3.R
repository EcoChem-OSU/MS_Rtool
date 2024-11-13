###############################################################################
## Title: PatRoon - feature group only
###############################################################################
## version:5.3
## Date: September 2024
## Author: Boris Droz 
## Modified from Tutorial and Handbook on https://github.com/rickhelmus/patRoon
## Depends:
##        R(>=4.3.2)
##        patRoon(>=2.3.0)
##        Java 21 download at https://www.oracle.com/java/technologies/downloads/#jdk21-windows
###############################################################################
## Description:
###############
## Performed peak picking
## 
## Should have a input file with a raw, mzxml folder and the sample list (csv)
## need to create an empty output folder
##
################################################################################
## Parameter -- MODIFIED IF NEEDED
############
## path
# workPath <- "D:/Patroon_NTS"
workPath <- "C:/Users/drozditb/Documents/OSU_data_analysis/WWT_testcode"

## Input data - 
sample.list <- "sample_list_testOMS.csv"

## Optimized XCMS parameters for peak picking
opt.ppm = 25
opt.pw = c(3, 143) # peak width min and max

## Parameter for filtering check patroon  help(filter)
min.intensity.thr = 200## absMinIntensity, typical range between 100 - 1000
rp.feature = 1 #relMinReplicateAbundance 
bk.sa.thr = 3 # blankThreshold - never go under 3

## Adduct and formula search parameter
adduct <- "[M+H]+"
                      
################################################################################
## load library -- DO NOT MODIFIED
library(patRoon) # v2.3.0
###############################################################################
###############################################################################

## DO NOT MODIFIED BELOW JUST RUN THE CODE ####

##   ----  SCRIPT START HERE ----
################################################################################

# ##############################################################################
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
#######
#######  SCRIPT START HERE
###########################################################################
# -------------------------
# initialization
# ------------------------
# Set outpath folder
# date <- Sys.Date()
# folder <- paste("/",date,"_NTA_patRoon", sep="")
#outpath <- creat.subDir(paste(workPath,"/output",sep=""), folder)
outpath <- paste(workPath,"output",sep="")
inpath <- paste(workPath,"input",sep="")

setwd(workPath) # set directory

# save parameter of the script
f.info <- paste(outpath,"/AA_INFO_RUN_README.txt",sep="")
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

## load data info
df <- read.csv(paste(workPath,"/input/",sample.list,sep=""),
                sep=",",header=TRUE)

anaInfo <- data.frame(cbind(path = df$path, 
                                analysis =df$filename,
                                group = df$group,
                                blank = df$blank) )
# -------------------------
# features
# -------------------------
# Find all features
param.xcms <- xcms::CentWaveParam(ppm = opt.ppm,
                              peakwidth = opt.pw,
                              snthresh = 10,
                              prefilter = c(3, 100),
                              noise = 0 )

fListPos <- findFeatures(anaInfo, "xcms3", param = param.xcms)

fList <- makeSet(fListPos, adducts = adduct)  

# performed RT alignement and group feature
fGroups <- groupFeatures(fList, "xcms3")
                      
# Basic rule based filtering
fGroups <- patRoon::filter(fGroups,  
                            absMinIntensity = min.intensity.thr, 
                            relMinReplicateAbundance = rp.feature, 
                            blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                            retentionRange = NULL, mzRange = NULL)

# -------------------------
# reporting
# -------------------------
## export averaged groupfeature as table
df.fGroups <- patRoon::as.data.table(fGroups, average = TRUE)
                
write.table(df.fGroups, file=paste(outpath,"/featureGroupsXCMS.csv", sep=""),
            append = FALSE, quote = TRUE, sep = ",",
            row.names = FALSE,col.names = TRUE )
