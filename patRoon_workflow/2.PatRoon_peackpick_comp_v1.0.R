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
## Performed peak picking with different algorythm
## 
## Should have a input file with a raw, mzxml folder and the sample list (csv)
## need to create an empty output folder

#### --- Warning some option are not available depending on the type of input file ---
# Algorythm	file types
# bruker	Bruker‘.d’ files
# openms	mzML
# xcms	mzML or mzXML
# xcms3	mzML or mzXML
# envipick	mzXML
# sirius	mzML or mzXML
# kpic2	mzML or mzXML
# safd	mzXML
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
ppm.dev = 10 # Allowed mass deviation (ppm) for trace detection
min.intensity.thr = 1000## absMinIntensity, typical range between 100 - 1000
rp.feature = 1 #relMinReplicateAbundance 
bk.sa.thr = 3 # blankThreshold - never go under 3

## Adduct and formula search parameter
# adduct <- "[M+H]+"
                      
################################################################################
## load library -- DO NOT MODIFIED
library(patRoon) # v2.3.0
library(stats)
library(BiocParallel)
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
cat( "PatRoon - feature peak picking comparions between algorythm", 
     file= f.info, append=TRUE, sep="\n")
cat( paste("SampleList: ", workPath,"/input/",sample.list,sep=""), 
     file= f.info, append=TRUE,sep="\n")
#cat( paste("XCMS_ppm:", opt.ppm), file= f.info, append=TRUE,sep="\n")
#cat( paste("XCMS_peakwidth:", opt.pw), file= f.info, append=TRUE,sep="\n")
cat( paste("Allowed mass deviation (ppm) for trace detection:", ppm.dev), file= f.info, append=TRUE,sep="\n")
cat( paste("absMinIntensity:", min.intensity.thr), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("relMinReplicateAbundance :", rp.feature), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("blankThreshold:", bk.sa.thr), file= f.info, append=TRUE,sep="\n")
# cat( paste("adduct:", adduct), file= f.info, append=TRUE,sep="\n")

## load data info
df <- read.csv(paste(workPath,"/input/",sample.list,sep=""),
                sep=",",header=TRUE)

anaInfo <- data.frame(cbind(path = df$path, 
                                analysis =df$filename,
                                group = df$group,
                                blank = df$blank) )

## exemple data set from patroon
# path <- list.files(patRoonData::exampleDataPath(), full.names=TRUE)
# 
# anaInfo <-data.frame(cbind(path = patRoonData::exampleDataPath(), 
#                            analysis = c(paste(rep("solvent-pos", 3),c("-1","-2","-3"), sep=""), 
#                                         paste(rep("standard-pos", 3),c("-1","-2","-3"), sep="")),
#                            group = c(rep("solvent-pos", 3), rep("standard-pos", 3)),
#                            blank = c(rep("solvent-pos", 6)) ))
  
# -------------------------
# features
# -------------------------
# Find all features  --setting for peak homogenized other 
##########################################################
#fListBRU <- findFeatures(anaInfo, "bruker")
fListOMS <- findFeatures(anaInfo, "openms", mzPPM=ppm.dev) # OpenMS, with default settings
fListXCMS <- findFeatures(anaInfo, "xcms",method = "centWave" ) # Old XCMS --> doesn't work
fListXCMS3 <- findFeatures(anaInfo, "xcms3",xcms::CentWaveParam(ppm = ppm.dev,
                                              peakwidth = c(1, 30)) )# NEW XCMS
fListEP <- findFeatures(anaInfo, "envipick", minint = min.intensity.thr, dmzgap=ppm.dev) # enviPick
fListSIRIUS <- findFeatures(anaInfo, "sirius") # SIRIUS
fListKPIC2 <- findFeatures(anaInfo, "kpic2", kmeans = TRUE, level = 1000) # KPIC2
fListSAF <- findFeatures(anaInfo, "safd", minInt = min.intensity.thr)

# fList <- makeSet(fListPos, adducts = adduct)  

# performed RT alignement and group feature
###########################################
# fGroupsBRU<- groupFeatures(fListOMS, "openms")
# fGroupsOMS <- groupFeatures(fListOMS, "openms") # OpenMS grouping, default settings
# fGroupsXCMS <- groupFeatures(fListXCMS, "xcms")
# fGroupsXCMS3 <- groupFeatures(fListXCMS3, "xcms3")
# fGroupsSIRIUS <- groupFeatures(fListSIRIUS, "sirius")
# fGroupsKPIC2 <- groupFeatures(fListKPIC2, "kpic2")

# do grouping with same algo
############################
fGroupsOMS <- groupFeatures(fListOMS, "xcms3")
fGroupsXCMS <- groupFeatures(fListXCMS, "xcms3")
fGroupsXCMS3 <- groupFeatures(fListXCMS3, "xcms3")
fGroupsEP <- groupFeatures(fListEP, "xcms3")
fGroupsSIRIUS <- groupFeatures(fListSIRIUS, "xcms3")
fGroupsKPIC2 <- groupFeatures(fListKPIC2, "xcms3")
fGroupsSAF <- groupFeatures(fListSAF, "xcms3")

## Do similar filtering for all data
####################################
#                       
# Basic rule based filtering
fGroupsOMS <- patRoon::filter(fGroupsOMS,  
                            absMinIntensity = min.intensity.thr, 
                            relMinReplicateAbundance = rp.feature, 
                            blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                            retentionRange = NULL, mzRange = NULL)

fGroupsXCMS <- patRoon::filter(fGroupsXCMS,  
                              absMinIntensity = min.intensity.thr, 
                              relMinReplicateAbundance = rp.feature, 
                              blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                              retentionRange = NULL, mzRange = NULL)

fGroupsXCMS3 <- patRoon::filter(fGroupsXCMS3,  
                              absMinIntensity = min.intensity.thr, 
                              relMinReplicateAbundance = rp.feature, 
                              blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                              retentionRange = NULL, mzRange = NULL)

fGroupsEP <- patRoon::filter(fGroupsEP,  
                              absMinIntensity = min.intensity.thr, 
                              relMinReplicateAbundance = rp.feature, 
                              blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                              retentionRange = NULL, mzRange = NULL)

fGroupsSIRIUS <- patRoon::filter(fGroupsSIRIUS,  
                              absMinIntensity = min.intensity.thr, 
                              relMinReplicateAbundance = rp.feature, 
                              blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                              retentionRange = NULL, mzRange = NULL)

fGroupsKPIC2 <- patRoon::filter(fGroupsKPIC2,  
                              absMinIntensity = min.intensity.thr, 
                              relMinReplicateAbundance = rp.feature, 
                              blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                              retentionRange = NULL, mzRange = NULL)

fGroupsSAF <- patRoon::filter(fGroupsSAF,  
                              absMinIntensity = min.intensity.thr, 
                              relMinReplicateAbundance = rp.feature, 
                              blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                              retentionRange = NULL, mzRange = NULL)

# -------------------------
# reporting
# -------------------------
## export averaged groupfeature as table
df.fGroupsOMS <- as.data.table(fGroupsOMS)
write.table(df.fGroupsOMS, file=paste(outpath,"/featureGroups_OMS.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

df.fGroupsXCMS <- as.data.table(fGroupsXCMS)
write.table(df.fGroupsXCMS, file=paste(outpath,"/featureGroups_XCMS.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

df.fGroupsEP <- as.data.table(fGroupsEP)
write.table(df.fGroupsEP, file=paste(outpath,"/featureGroups_envipick.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

df.fGroupsSIRIUS <- as.data.table(fGroupsSIRIUS)
write.table(df.fGroupsSIRIUS, file=paste(outpath,"/featureGroups_SIRIUS.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

df.fGroupsKPIC2 <- as.data.table(fGroupsKPIC2)
write.table(df.fGroupsKPIC2, file=paste(outpath,"/featureGroups_KPIC2.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

df.fGroupsSAF <- as.data.table(fGroupsSAF)
write.table(df.fGroupsXCMS3, file=paste(outpath,"/featureGroups_SAF.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

## export averaged groupfeature as table
df.fGroupsOMS <- as.data.table(fGroupsOMS, average = TRUE)
write.table(df.fGroupsOMS, file=paste(outpath,"/featureGroups_OMS_averaged.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

df.fGroupsXCMS <- as.data.table(fGroupsXCMS, average = TRUE)
write.table(df.fGroupsXCMS, file=paste(outpath,"/featureGroups_XCMS_averaged.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

df.fGroupsEP <- as.data.table(fGroupsEP, average = TRUE)
write.table(df.fGroupsEP, file=paste(outpath,"/featureGroups_envipick_averaged.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

df.fGroupsSIRIUS <- as.data.table(fGroupsSIRIUS, average = TRUE)
write.table(df.fGroupsSIRIUS, file=paste(outpath,"/featureGroups_SIRIUS_averaged.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

df.fGroupsKPIC2 <- as.data.table(fGroupsKPIC2, average = TRUE)
write.table(df.fGroupsKPIC2, file=paste(outpath,"/featureGroups_KPIC2_averaged.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

df.fGroupsSAF <- as.data.table(fGroupsSAF, average = TRUE)
write.table(df.fGroupsXCMS3, file=paste(outpath,"/featureGroups_SAF_averaged.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )




