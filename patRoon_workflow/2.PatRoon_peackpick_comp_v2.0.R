###############################################################################
## Title: PatRoon - feature peak picking comp
###############################################################################
## version: 2.0
## Date: October 2026
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
workPath <- "E:/UKL_2023"

## Input data - 
sample.list <- "20251024_sample_list_UKL2023_mzml.csv"

in.fn.type <-"pos_mzML" ## just text to make the ouputfolder name

# ## Optimized parameters for peak picking
## --> current version used defaut parameter for each peak picking

## Parameter for filtering check patroon  help(filter)
ppm.dev = 10 # Allowed mass deviation (ppm) for trace detection
min.intensity.thr = 100## absMinIntensity, typical range between 100 - 1000
rp.feature = 1 #relMinReplicateAbundance 
bk.sa.thr = 3 # blankThreshold - never go under 3

## Adduct and formula search parameter
adduct <- "[M+H]+"
                      
################################################################################
## load library -- DO NOT MODIFIED
library(patRoon) # v2.3.0
library(stats)
library(BiocParallel)

patRoon::verifyDependencies()
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
date <- Sys.Date()
inpath <- paste(workPath,"/input",sep="")
proj <- sub(".*/", "", workPath)

folder <- paste("/",date,"_peak_picking_comp_", proj,"_",in.fn.type,sep="")
outpath <- creat.subDir(paste(workPath,"/output",sep=""), folder)
#outpath <- paste(workPath,"/output",sep="")

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
cat("current version used defaut parameter for each peak picking", file= f.info, append=TRUE,sep="\n")
cat( paste("Allowed mass deviation (ppm) for trace detection:", ppm.dev), file= f.info, append=TRUE,sep="\n")
cat( paste("absMinIntensity:", min.intensity.thr), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("relMinReplicateAbundance :", rp.feature), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("blankThreshold:", bk.sa.thr), file= f.info, append=TRUE,sep="\n")
cat( paste("adduct:", adduct), file= f.info, append=TRUE,sep="\n")
cat( paste("file_type:", in.fn.type), file= f.info, append=TRUE,sep="\n")

## load data info
df <- read.csv(paste(workPath,"/input/",sample.list,sep=""),
                sep=",",header=TRUE)

anaInfo <- data.frame(cbind(path = df$path, 
                                analysis =df$filename,
                                replicate = df$group,
                                blank = df$blank) )
##test data
# anaInfo <- anaInfo[1:8,]

# check file extentsion to make decision
fn.ext <- unique(tools::file_ext(list.files(anaInfo$path, full.names = TRUE)))

mzxml <- fn.ext=="mzXML"
mzml <- fn.ext=="mzML"

feat.summ <-data.frame(NULL) 

## exemple data set from patroon
# path <- list.files(patRoonData::exampleDataPath(), full.names=TRUE)
# 
# anaInfo <-data.frame(cbind(path = patRoonData::exampleDataPath(), 
#                          analysis = c(paste(rep("solvent-pos", 3),c("-1","-2","-3"), sep=""), 
#                                      paste(rep("standard-pos", 3),c("-1","-2","-3"), sep="")),
#                           group = c(rep("solvent-pos", 3), rep("standard-pos", 3)),
#                          blank = c(rep("solvent-pos", 6)) ))
  
# -------------------------
# Brucker
# -------------------------
#fListBRU <- findFeatures(anaInfo, "bruker")

if (mzml) {
          # -------------------------
          # OpenMS
          # -------------------------
          fListOMS <- findFeatures(anaInfo, "openms", mzPPM=ppm.dev) # OpenMS, with default settings
          fGroupsOMS <- groupFeatures(fListOMS, "openms")
          
          # Basic rule based filtering
          fGroupsOMS <- patRoon::filter(fGroupsOMS,
                                      absMinIntensity = min.intensity.thr,
                                      relMinReplicateAbundance = rp.feature,
                                      blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                                      retentionRange = NULL, mzRange = NULL)
          
          df.fGroupsOMS <- as.data.table(fGroupsOMS,areas = TRUE)
          write.table(df.fGroupsOMS, file=paste(outpath,"/featureGroups_OMS.txt", sep=""),
                      append = FALSE, quote = FALSE, sep = "\t",
                      row.names = FALSE,col.names = TRUE )
          
          ## export averaged groupfeature as table
          df.fGroupsOMS <- as.data.table(fGroupsOMS, average = TRUE, areas = TRUE)
          write.table(df.fGroupsOMS, file=paste(outpath,"/featureGroups_OMS_averaged.txt", sep=""),
                      append = FALSE, quote = FALSE, sep = "\t",
                      row.names = FALSE,col.names = TRUE )
          
          feat.summ <-data.frame(c(feat.summ, nb_feature_OpenMS= nrow(df.fGroupsOMS))) 
          }

if (mzml|mzxml) {
          # -------------------------
          # XCMS
          # -------------------------
          fListXCMS <- findFeatures(anaInfo, "xcms",method = "centWave" ) # Old XCMS
          fListXCMS <- makeSet(fListXCMS, adducts = adduct)  
          fGroupsXCMS <- groupFeatures(fListXCMS, "xcms3")
          
          fGroupsXCMS <- groupFeatures(fListXCMS, "xcms3")
          
          fGroupsXCMS <- patRoon::filter(fGroupsXCMS,  
                                         absMinIntensity = min.intensity.thr, 
                                         relMinReplicateAbundance = rp.feature, 
                                         blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                                         retentionRange = NULL, mzRange = NULL)
          df.fGroupsXCMS <- as.data.table(fGroupsXCMS, areas = TRUE)
          write.table(df.fGroupsXCMS, file=paste(outpath,"/featureGroups_XCMS.txt", sep=""),
                      append = FALSE, quote = FALSE, sep = "\t",
                      row.names = FALSE,col.names = TRUE )
          
          df.fGroupsXCMS <- as.data.table(fGroupsXCMS, average = TRUE, areas = TRUE)
          write.table(df.fGroupsXCMS, file=paste(outpath,"/featureGroups_XCMS_averaged.txt", sep=""),
                      append = FALSE, quote = FALSE, sep = "\t",
                      row.names = FALSE,col.names = TRUE )
          
          feat.summ <-data.frame(c(feat.summ, nb_feature_XCMS= nrow(df.fGroupsXCMS))) 
        }

if (mzml|mzxml) {
        # -------------------------
        # XCMS3
        # -------------------------
        fListXCMS3 <- findFeatures(anaInfo, "xcms3",xcms::CentWaveParam(ppm = ppm.dev,
                                                    peakwidth = c(1, 30)) )# NEW XCMS
        fListXCMS3 <- makeSet(fListXCMS3, adducts = adduct)  
        fGroupsXCMS3 <- groupFeatures(fListXCMS3, "xcms3")
        
        fGroupsXCMS3 <- patRoon::filter(fGroupsXCMS3,  
                                        absMinIntensity = min.intensity.thr, 
                                        relMinReplicateAbundance = rp.feature, 
                                        blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                                        retentionRange = NULL, mzRange = NULL)
        
        df.fGroupsXCMS3 <- as.data.table(fGroupsXCMS3, areas = TRUE)
        write.table(df.fGroupsXCMS3, file=paste(outpath,"/featureGroups_XCMS3.txt", sep=""),
                    append = FALSE, quote = FALSE, sep = "\t",
                    row.names = FALSE,col.names = TRUE )
        
        df.fGroupsXCMS3 <- as.data.table(fGroupsXCMS3, average = TRUE, areas = TRUE)
        write.table(df.fGroupsXCMS3, file=paste(outpath,"/featureGroups_XCMS3_averaged.txt", sep=""),
                    append = FALSE, quote = FALSE, sep = "\t",
                    row.names = FALSE,col.names = TRUE )
        
        feat.summ <-data.frame(c(feat.summ, nb_feature_XCMS3= nrow(df.fGroupsXCMS3))) 
        }

if (mzxml) {
          # -------------------------
          # enviPick
          # -------------------------
          fListEP <- findFeatures(anaInfo, "envipick", minint = min.intensity.thr, dmzgap=ppm.dev) # enviPick
          fListEP <- makeSet(fListEP, adducts = adduct) 
          fGroupsEP <- groupFeatures(fListEP, "xcms3")
          
          fGroupsEP <- patRoon::filter(fGroupsEP,  
                                       absMinIntensity = min.intensity.thr, 
                                       relMinReplicateAbundance = rp.feature, 
                                       blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                                       retentionRange = NULL, mzRange = NULL)
          
          df.fGroupsEP <- as.data.table(fGroupsEP, areas = TRUE)
          write.table(df.fGroupsEP, file=paste(outpath,"/featureGroups_envipick.txt", sep=""),
                      append = FALSE, quote = FALSE, sep = "\t",
                      row.names = FALSE,col.names = TRUE )
          
          df.fGroupsEP <- as.data.table(fGroupsEP, average = TRUE, areas = TRUE)
          write.table(df.fGroupsEP, file=paste(outpath,"/featureGroups_envipick_averaged.txt", sep=""),
                      append = FALSE, quote = FALSE, sep = "\t",
                      row.names = FALSE,col.names = TRUE )
          
          feat.summ <-data.frame(c(feat.summ, nb_feature_enviPick= nrow(df.fGroupsEP))) 
        }
# 
# if (mzml|mzxml) {
#             # -------------------------
#             # SIRIUS
#             # ------------------------
#             fListSIRIUS <- findFeatures(anaInfo, "sirius") # SIRIUS
#             # fListSIRIUS <- makeSet(fListSIRIUS, adducts = adduct)
#             # fGroupsSIRIUS <- groupFeatures(fListSIRIUS, "xcms3")
# 
#             fGroupsSIRIUS <- groupFeatures(anaInfo, "sirius")
# 
#             fGroupsSIRIUS <- patRoon::filter(fGroupsSIRIUS,
#                                              absMinIntensity = min.intensity.thr,
#                                              relMinReplicateAbundance = rp.feature,
#                                              blankThreshold = bk.sa.thr, removeBlanks = TRUE,
#                                              retentionRange = NULL, mzRange = NULL)
# 
#             df.fGroupsSIRIUS <- as.data.table(fGroupsSIRIUS, areas = TRUE)
#             write.table(df.fGroupsSIRIUS, file=paste(outpath,"/featureGroups_SIRIUS.txt", sep=""),
#                         append = FALSE, quote = FALSE, sep = "\t",
#                         row.names = FALSE,col.names = TRUE )
# 
#             df.fGroupsSIRIUS <- as.data.table(fGroupsSIRIUS, average = TRUE, areas = TRUE)
#             write.table(df.fGroupsSIRIUS, file=paste(outpath,"/featureGroups_SIRIUS_averaged.txt", sep=""),
#                         append = FALSE, quote = FALSE, sep = "\t",
#                         row.names = FALSE,col.names = TRUE )
# 
#             feat.summ <-data.frame(c(feat.summ, nb_feature_SIRIUS= nrow(df.fGroupsSIRIUS)))
#         }

if (mzml|mzxml) {
            # -------------------------
            # KPIC2
            # -------------------------
            fListKPIC2 <- findFeatures(anaInfo, "kpic2", kmeans = TRUE, level = 1000) # KPIC2
            fListKPIC2 <- makeSet(fListKPIC2, adducts = adduct) 
            fGroupsKPIC2 <- groupFeatures(fListKPIC2, "xcms3")
            
            fGroupsKPIC2 <- patRoon::filter(fGroupsKPIC2,  
                                            absMinIntensity = min.intensity.thr, 
                                            relMinReplicateAbundance = rp.feature, 
                                            blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                                            retentionRange = NULL, mzRange = NULL)
            
            df.fGroupsKPIC2 <- as.data.table(fGroupsKPIC2, areas = TRUE)
            write.table(df.fGroupsKPIC2, file=paste(outpath,"/featureGroups_KPIC2.txt", sep=""),
                        append = FALSE, quote = FALSE, sep = "\t",
                        row.names = FALSE,col.names = TRUE )
            
            df.fGroupsKPIC2 <- as.data.table(fGroupsKPIC2, average = TRUE, areas = TRUE)
            write.table(df.fGroupsKPIC2, file=paste(outpath,"/featureGroups_KPIC2_averaged.txt", sep=""),
                        append = FALSE, quote = FALSE, sep = "\t",
                        row.names = FALSE,col.names = TRUE )
            
            feat.summ <-data.frame(c(feat.summ, nb_feature_KPIC2= nrow(df.fGroupsKPIC2)))
            }

# if (mzxml) {
#             # -------------------------
#             # SAF
#             # -------------------------
#             fListSAF <- findFeatures(anaInfo, "safd", minInt = min.intensity.thr)
#             fListSAF <- makeSet(fListSAF, adducts = adduct) 
#             fGroupsSAF <- groupFeatures(fListSAF, "xcms3")
#             
#             fGroupsSAF <- patRoon::filter(fGroupsSAF,  
#                                           absMinIntensity = min.intensity.thr, 
#                                           relMinReplicateAbundance = rp.feature, 
#                                           blankThreshold = bk.sa.thr, removeBlanks = TRUE,
#                                           retentionRange = NULL, mzRange = NULL)
#             
#             df.fGroupsSAF <- as.data.table(fGroupsSAF, areas = TRUE)
#             write.table(df.fGroupsXCMS3, file=paste(outpath,"/featureGroups_SAF.txt", sep=""),
#                         append = FALSE, quote = FALSE, sep = "\t",
#                         row.names = FALSE,col.names = TRUE )
#             
#             df.fGroupsSAF <- as.data.table(fGroupsSAF, average = TRUE, areas = TRUE)
#             write.table(df.fGroupsXCMS3, file=paste(outpath,"/featureGroups_SAF_averaged.txt", sep=""),
#                         append = FALSE, quote = FALSE, sep = "\t",
#                         row.names = FALSE,col.names = TRUE )
#             
#             feat.summ <-data.frame(c(feat.summ, nb_feature_SAF= nrow(df.fGroupsSAF)))
#             }

write.table(feat.summ, file=paste(outpath,"/feature_summary.txt", sep=""),
                        append = FALSE, quote = FALSE, sep = "\t",
                        row.names = FALSE,col.names = TRUE )



