###############################################################################
## Title: PatRoon - feature peak picking comp
###############################################################################
## version:1.0
## Date: July 2026
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

# ## Optimized XCMS parameters for peak picking
# opt.ppm = 25
# opt.pw = c(3, 143) # peak width min and max
## --> current version used defaut parameter for ech peak picking

## Parameter for filtering check patroon  help(filter)
min.intensity.thr = 1000## absMinIntensity, typical range between 100 - 1000
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
outpath <- paste(workPath,"/output",sep="")
inpath <- paste(workPath,"/input",sep="")

setwd(workPath) # set directory

# save parameter of the script
f.info <- paste(outpath,"/AA_INFO_RUN_README.txt",sep="")
cat( paste("*** patRoon parameter for the run....", Sys.Date()), 
     file= f.info, append=TRUE, sep="\n")
cat( "#########################################################", 
     file= f.info, append=TRUE,sep="\n")
cat( paste("SampleList: ", workPath,"/input/",sample.list,sep=""), 
     file= f.info, append=TRUE,sep="\n")
#cat( paste("XCMS_ppm:", opt.ppm), file= f.info, append=TRUE,sep="\n")
#cat( paste("XCMS_peakwidth:", opt.pw), file= f.info, append=TRUE,sep="\n")
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
###################
fListBRU <- findFeatures(anaInfo, "bruker")
fListOMS <- findFeatures(anaInfo, "openms") # OpenMS, with default settings
fListXCMS <- findFeatures(anaInfo, "xcms", ppm = 10) # XCMS
fListXCMS3 <- findFeatures(anaInfo, "xcms3", 
                           CentWaveParam(peakwidth = c(5, 15))) # XCMS3
fListEP <- findFeatures(anaInfo, "envipick", minint = 1E3) # enviPick
fListSIRIUS <- findFeatures(anaInfo, "sirius") # SIRIUS
fListKPIC2 <- findFeatures(anaInfo, "kpic2", kmeans = TRUE, level = 1E4) # KPIC2
fListSAF <- findFeatures(anaInfo, "safd")



# fList <- makeSet(fListPos, adducts = adduct)  

# performed RT alignement and group feature
fGroupsOMS <- groupFeatures(fListOMS, "openms") # OpenMS grouping, default settings
fGroupsXCMS <- groupFeatures(fListXCMS, "openms", maxGroupRT = 6) # group XCMS features with OpenMS, adjusted grouping parameter
# group enviPick features with XCMS3, disable minFraction
fGroupsXCMS3 <- groupFeatures(fListEP, "xcms3",
                             xcms::PeakDensityParam(sampleGroups = analInfo$replicate,
                                                    minFraction = 0))
# group with KPIC2 and set some custom grouping/aligning parameters
fGroupsKPIC2 <- groupFeatures(fListKPIC2, "kpic2", groupArgs = list(tolerance = c(0.002, 18)),
                              alignArgs = list(move = "loess"))
# greedy algorithm with custom tolerances and weights
fGroupsGreedy <- groupFeatures(fListXCMS3, "greedy", rtWindow = 5, mzWindow = 0.003,
                               scoreWeights = c(retention = 0.5, mz = 3, mobility = 1, intensity = 1))
fGroupsSIRIUS <- groupFeatures(anaInfo, "sirius") # find/group features with SIRIUS

# # export raw data for control
# df.fGroups <- as.data.table(fGroups)
# df.fGroups <- na.omit(df.fGroups)
# 
# write.table(df.fGroups, file=paste(outpath,"/raw_aligned_grouped.txt", sep=""),
#             append = FALSE, quote = FALSE, sep = "\t",
#             row.names = FALSE,col.names = TRUE )
#                       
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
df.fGroups <- as.data.table(fGroups)

write.table(df.fGroups, file=paste(outpath,"/featureGroups.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

## export averaged groupfeature as table
df.fGroups <- as.data.table(fGroups, average = TRUE)

write.table(df.fGroups, file=paste(outpath,"/featureGroups_averaged.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )


