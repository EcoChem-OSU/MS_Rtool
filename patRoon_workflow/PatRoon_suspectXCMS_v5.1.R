###############################################################################
## Title: PatRoon - suspect screening 
###############################################################################
## version:5.0
## Date: August 2024
## Author: Boris Droz 
## Modified from Tutorial and Handbook on https://github.com/rickhelmus/patRoon
## Depends:
##        R(>=4.3.2)
##        patRoon(>=2.3.0)
##        Java 21 download at https://www.oracle.com/java/technologies/downloads/#jdk21-windows
###############################################################################
## Description:
###############
## Performed suspect analysis with MS2 library search and insilico Metfrag
## 
## Should have a input file with a raw, mzxml folder and the sample list (csv)
##
## REQUIERE:
## YOU WILL NEED to update the suspect list, MetfragPuchem DB and MS2 library 
## before starting
## just make sure to used the right R version and the R assoiated R-library
## 
###############################################################################
## Parameter -- MODIFIED IF NEEDED
############
## path
# workPath <- "D:/Patroon_NTS"
workPath <- "C:/Users/drozditb/Documents/OSU_data_analysis/NIST_20240501/"

## Input data - 
sample.list <- "sample_list_nist.csv"

## Optimized XCMS parameters for peak picking
opt.ppm = 25
opt.pw = c(3, 143) # peak width min and max

## Parameter for filtering check patroon  help(filter)
min.intensity.thr = ## absMinIntensity, typical range between 100 - 1000
rp.feature = 1 #relMinReplicateAbundance 
bk.sa.thr = 3 # blankThreshold 

## Adduct and formula search parameter
adduct <- "[M+H]+"
form.ele <- "CHNOPSFCl" # which element is considered in formula search
                        # "CHNOPSCl" +BrF are common considered elements for pollutant

################################################################################
## load library
library(patRoon) # v2.3.0
# patRoon::verifyDependencies() ## Verifying the installation of dependency 
# library(xcms) 
library(BiocParallel) # to parralel calculation
library(dplyr)
library(webchem)

## Parameter folder link
## generally set once and do not need to be change
# PatRoon.directory -> where suspect, MS2, ... are
PatRoon.dir <- "C:/Users/drozditb/Documents/general_library/patRoon-install"

## load suspect list 
fns <- paste(PatRoon.dir,"/suspect_list/neg_nist_mds2-2387.csv",sep="")
## set path to 
options(patRoon.path.obabel = "C:/Program Files/OpenBabel-3.1.1/") # open babel exe
options(patRoon.path.MetFragCL = 
          paste(PatRoon.dir,"/MetFrag/MetFragCommandLine-2.5.0.jar", sep="")) # metfrag exe
options(patRoon.path.MetFragPubChemLite = 
          paste(PatRoon.dir,"/MetFrag/PubChem_OECDPFAS_largerPFASparts_20220324.csv",sep="") )
          
# # select and load MS2 library
mslibraryNIST <- loadMSLibrary(paste(PatRoon.dir,"/MS2_library/MassBank_NIST.msp", sep=""), "msp")
mslibraryRIKEN <- loadMSLibrary(paste(PatRoon.dir,"/MS2_library/MassBank_RIKEN.msp", sep=""), "msp") 
# # 
mslibraryM <- merge(mslibraryNIST, mslibraryRIKEN) # merge library

mslibraryM <- loadMSLibrary(paste(PatRoon.dir,"/MS2_library/DIMSpecForPFAS_2023-10-03.msp", sep=""), "msp")

###############################################################################
###############################################################################

## DO NOT MODIFIED BELOW JUST RUN THE CODE ####

##   ----  SCRIPT START HERE ----
###############################################################################

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
#######
#######  SCRIPT START HERE
###########################################################################
# -------------------------
# initialization
# ------------------------
# Set outpath folder
date <- Sys.Date()
folder <- paste("/",date,"_NTA_patRoon", sep="")
outpath <- creat.subDir(paste(workPath,"/output",sep=""), folder)
inpath <- paste(workPath,"input",sep="")

setwd(workPath) # set directory

# save parameter of the script




## load data info
df <- read.csv(paste(workPath,"/input/",sample.list,sep=""),
                sep=",",header=TRUE)

anaInfo <- data.frame(cbind(path = df$path, #paste(df$path,"/",df$folder,sep=""),
                                analysis =df$filename,
                                group = df$group,
                                blank = df$blank) )

# -------------------------
# features
# -------------------------
# Find all features
# NOTE: see the XCMS manual for many more options
param.xcms <- xcms::CentWaveParam(ppm = opt.ppm,
                              peakwidth = opt.pw,
                              snthresh = 10,
                              prefilter = c(3, 100),
                              noise = 0 ) # noise could be sensitive do not change  

fListPos <- findFeatures(anaInfo, "xcms3", param = param.xcms)

fList <- makeSet(fListPos, adducts = adduct)  

df.fList <- as.data.table(fList)
df.fList <- na.omit(df.fList)

write.table(df.fList, file=paste(outpath,"/raw_unaligned_ungrouped.csv", sep=""),
            append = FALSE, quote = TRUE, sep = ",",
            row.names = FALSE,col.names = TRUE )

# performed RT alignement and group feature
fGroups <- groupFeatures(fList, "xcms3")
                         
# save raw  data for control
## export data
df.fGroups <- as.data.table(fGroups)
df.fGroups <- na.omit(df.fGroups)

write.table(df.fGroups, file=paste(outpath,"/raw_aligned_grouped.csv", sep=""),
            append = FALSE, quote = TRUE, sep = ",",
            row.names = FALSE,col.names = TRUE )

# Basic rule based filtering
fGroups <- patRoon::filter(fGroups,  
                            absMinIntensity = min.intensity.thr, 
                            relMinReplicateAbundance = rp.feature, 
                            #maxReplicateIntRSD = 0.5,
                            blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                            retentionRange = NULL, mzRange = NULL)

# -------------------------
# reporting
# -------------------------
## export groupfeature as table
df.fGroups <- as.data.table(fGroups)

write.table(df.fGroups, file=paste(outpath,"/featureGroups.csv", sep=""),
            append = FALSE, quote = TRUE, sep = ",",
            row.names = FALSE,col.names = TRUE )

## export averaged groupfeature as table
df.fGroups <- as.data.table(fGroups, average = TRUE)
                
write.table(df.fGroups, file=paste(outpath,"/featureGroups_averaged.csv", sep=""),
            append = FALSE, quote = TRUE, sep = ",",
            row.names = FALSE,col.names = TRUE )
# -------------------------
# Annotation 
# -------------------------
# Retrieve MS peak lists 
avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.002)
mslists <- generateMSPeakLists(fGroups, "mzr", 
                               maxMSRtWindow = 5, 
                               precursorMzWindow = 4,
                               avgFeatParams = avgPListParams, 
                               avgFGroupParams = avgPListParams)

# Rule based filtering of MS peak lists.
mslists <- patRoon::filter(mslists, 
                           withMSMS = TRUE, absMSIntThr = NULL, 
                            absMSMSIntThr = NULL, relMSIntThr = NULL, 
                            relMSMSIntThr = 0.05,
                            topMSPeaks = NULL, topMSMSPeaks = 25)

formulas <- generateFormulasGenForm(fGroups, mslists,
                                    elements = form.ele, 
                                    relMzDev = 5,
                                    absAlignMzDev = 0.002,
                                    topMost = 10,
                                    calculateFeatures = TRUE,
                                    featThresholdAnn = 1,
                                    MSMode ="both")

## MetFrag -----
################
# Calculate compound structure candidates
# scoreTypes to be defined, e.g.: "fragScore","metFusionScore","...
compsMF <- generateCompounds(
                    fGroups, mslists,
              "metfrag",
              method = "CL",
              dbRelMzDev = 5 ,
              fragRelMzDev = 5,
              fragAbsMzDev = 0.002,
              database = "pubchemlite",
              setThresholdAnn = 1,
              scoreTypes = c("individualMoNAScore", "fragScore", "score"),
              scoreWeights = 1,
              maxCandidatesToStop = 100,
              timeoutRetries = 20 )

# Summary of MetFrag Results in a a Single Table
MFsummary <- as.data.table(compounds)
outputSummary <- paste(outpath, "MFsummary.csv", sep = "/")
write.csv(MFsummary, outputSummary)

# Annotation with the Library MS2 algorithm
#########################################
compsLib <- generateCompounds(fGroups, mslists, "library", 
                              MSLibrary = mslibraryM, minSim = 0.4)

## Suspect screening 
#####################
dat <- read.csv(fns, header=TRUE) #open suspect list
    
suspects <- data.frame(name = dat$ID,
                         SMILES =dat$SMILES,
                         stringsAsFactors = FALSE) 
  
fGroupsSusp <- screenSuspects(fGroups, suspects, 
                              # rtWindow = 60, #in sec
                              mzWindow = 0.005,
                              onlyHits = TRUE) # remove any non-hits immediately
# -------------------------
# Final Annotation 
# -------------------------
fGroupsSusp_MF <- annotateSuspects(
                    fGroupsSusp,
                    MSPeakLists = mslists,
                    formulas = formulas,
                    compounds = compsMF,
                    absMzDev = 0.005,
                    specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
                    checkFragments = c("mz", "formula", "compound"),
                    formulasNormalizeScores = "max",
                    compoundsNormalizeScores = "max") #

fGroupsSusp_Lib <- annotateSuspects(
                    fGroupsSusp,
                    MSPeakLists = mslists,
                    formulas = formulas,
                    compounds = compsLib,
                    absMzDev = 0.005,
                    specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
                    checkFragments = c("mz", "formula", "compound"),
                    formulasNormalizeScores = "max",
                    compoundsNormalizeScores = "max") #

# final filtration
fGroupsSusp_MF <- patRoon::filter(fGroupsSusp_MF,
                        selectHitsBy="level", # take the highest Level conf. if same suspect match to multiple feature groups
                        selectBestFGroups=TRUE) # keep only best identification score

fGroupsSusp_Lib <- patRoon::filter(fGroupsSusp_Lib,
                                  selectHitsBy="level", 
                                  selectBestFGroups=TRUE)
# -------------------------
# reporting
# -------------------------
## export groupfeature as table
df.fGroupsSusp_MF <-as.data.table(fGroupsSusp_MF, average = TRUE, 
                               collapseSuspects = NULL,
                               onlyHits = TRUE) 

df.fGroupsSusp_Lib <-as.data.table(fGroupsSusp_Lib, average = TRUE, 
                                  collapseSuspects = NULL,
                                  onlyHits = TRUE)

# combine the two suspect tables, delete duplicate with keep higest LC rank
df.fGroupsSusp <- rbind( cbind(df.fGroupsSusp_MF, Source="metfrag"), 
			          	cbind(df.fGroupsSusp_Lib, Source="Lib"), fill=TRUE )

# Sort merged table based on Value column in descending order
df.fGroupsSusp <- df.fGroupsSusp[order(df.fGroupsSusp$susp_estIDLevel, decreasing = FALSE), ]
df.fGroupsSusp <- df.fGroupsSusp[order(df.fGroupsSusp$group, decreasing = TRUE), ]

# Identify fully identical rows 
fully_identical <- duplicated(df.fGroupsSusp$group) & duplicated(df.fGroupsSusp$susp_estIDLevel)

# Select the first occurrence of each group
max_conf <- !duplicated(df.fGroupsSusp$group)

# Update the "Source" column for the row with the highest level of conf among fully identical rows
df.fGroupsSusp$Source <- ifelse(fully_identical & !max_conf, "", paste(df.fGroupsSusp$Source, c(df.fGroupsSusp$Source[-1], NA) ))

# Keep only the unique rows based on ID, keeping the highest value
df.fGroupsSusp <- df.fGroupsSusp[max_conf, ]

# get chemical names based on InChIKey
# check if present in pubchem
dmol <-NULL
for (i in 1:nrow(df.fGroupsSusp) )
	{
	  get.id <- df.fGroupsSusp$susp_InChIKey[i]
	  
      fcid <- try ( get_cid(get.id, from =  'inchikey', match='first', verbose = TRUE, arg = NULL) , silent = TRUE )
      
      # if present in pubchem get
     if ( length(fcid$cid)==2 | fcid$cid==0 | is.na(fcid$cid) )  { 
       dmol <- rbind(dmol, data.frame(CID=NA, IUPACName=NA,XLogP=NA))
       }else{
         # data properties from pubchem
         prop <- pc_prop(fcid$cid, properties =  c('IUPACName','XLogP')
                         , verbose = TRUE)
         # check length of prop
         if (length(prop)==3){dmol <- rbind(dmol,prop)
           }else{  if (all(names(prop) == c("CID","XLogP") ))
                     { dmol <- rbind(dmol,c(prop$CID,IUPACName=NA ,prop$XLogP) )
                    }else{ dmol <- rbind(dmol,c(prop,XLogP=NA) ) }
              }
       }
    }

df.fGroupsSusp <- cbind(df.fGroupsSusp,dmol)

write.table(df.fGroupsSusp, file=paste(outpath, "/SuspectScreening_all.csv", sep=""),
            append = FALSE, quote = TRUE, sep = ",",
            row.names = FALSE,col.names = TRUE )


# create a working table for sample only 
df.fGroupsSusp <- data.frame(df.fGroupsSusp)
names(df.fGroupsSusp) <- sub("^X", "", names(df.fGroupsSusp)) # rename if name start by X -- 

data <- df.fGroupsSusp[ ,names(df.fGroupsSusp) %in% unique(df$group[df$sampletype=="SA"]) ]

# select minimal info
df.data <- cbind(data.frame(cbind(group= df.fGroupsSusp$group,
                  ret=df.fGroupsSusp$ret,
                 mz=df.fGroupsSusp$mz,
                 InChIKey=df.fGroupsSusp$susp_InChIKey,
                 estIDLevel=df.fGroupsSusp$susp_estIDLevel,
                 IUPACName=df.fGroupsSusp$IUPACName,
                 XLogP=df.fGroupsSusp$LogP)), data)

# final check if all value = zero 
index <- rowSums( df.data[,8:ncol(df.data)]) >0 #######CHECK the 8 here####

df.data <- df.data[index,]

write.table(df.data, file=paste(outpath, "/SuspectScreening_sample.csv", sep=""),
            append = FALSE, quote = TRUE, sep = ",",
            row.names = FALSE,col.names = TRUE )
