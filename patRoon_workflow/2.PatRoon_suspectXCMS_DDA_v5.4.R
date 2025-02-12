###############################################################################
## Title: PatRoon - suspect screening 
###############################################################################
## version:5.4
## Date: January 2025
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
## just make sure to used the right R version and the R associated R-library
## 
################################################################################
## Parameter -- MODIFIED IF NEEDED
############
## path
# workPath <- "D:/Patroon_NTS"
workPath <- "C:/Users/drozditb/Documents/OSU_data_analysis/deconvoltest_midcal"

## Input data - 
sample.list <- "sample_list_midcal_mzML.csv"

# check for ISTD - option are "YES" or "NO"
check.istd <- "NO"
istd.list <- "istd_list.csv"

## Optimized XCMS parameters for peak picking
opt.ppm = 25
opt.pw = c(3, 143) # peak width min and max

## Parameter for filtering check patroon  help(filter)
min.intensity.thr = 200## absMinIntensity, typical range between 100 - 1000
rp.feature = 1 #relMinReplicateAbundance 
bk.sa.thr = 3 # blankThreshold - never go under 3

# Mass defect filtering - option are "YES" or "NO"
MD.filter <- "YES" # used the suspect list to mass filtering.
#MD.minmax <-c(-0.49,0.50) # based on the merge suspect list (OCDE,NIST,NORMAN,2EPA)
#MD.minmax <-c(-0.25,0.1) # from Zwiener paper based on OECD suspect list

## Adduct and formula search parameter
adduct <- "[M-H]-"
form.ele <- "CHNOPSFClBrF" #  CHNOPSCl" +BrF are common considered elements for pollutant
                      
################################################################################
## load library -- DO NOT MODIFIED
library(patRoon) # v2.3.0
library(xcms) 
library(BiocParallel) 
library(dplyr)
library(webchem)
###############################################################################
## Parameter Path
## generally set once and do not need to be change
# PatRoon.directory -> where suspect, MS2, ... are
PatRoon.dir <- "C:/Users/drozditb/Documents/general_library/patRoon-install"

## set suspect list, MS2 and Metfrag all located in PatRoon.dir
fns <- paste(PatRoon.dir,"/suspect_list/neg_Targets_std_List_20240719_Peter_mod.csv",sep="")
MS2.lib <- c("Fluoros_2.5_editedV4.msp") 
fn.metfrag <- paste(PatRoon.dir,"/MetFrag/PubChem_OECDPFAS_largerPFASparts_20220324.csv",sep="")

## set path to -- GENERALLY DO NOT NEED TO MODIFY
options(patRoon.path.obabel = "C:/Program Files/OpenBabel-3.1.1/") # open babel exe
options(patRoon.path.MetFragCL = 
          paste(PatRoon.dir,"/MetFrag/MetFragCommandLine-2.5.0.jar", sep="")) # metfrag exe
options(patRoon.path.MetFragPubChemLite = fn.metfrag )
          
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
folder <- paste("/",date,"_NTA_patRoon", sep="")
outpath <- creat.subDir(paste(workPath,"/output",sep=""), folder)
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
if (check.istd =="YES"){
  cat( paste("ISTDList: ", workPath,"/input/",istd.list,sep=""), 
       file= f.info, append=TRUE,sep="\n")
}else{}
cat( paste("XCMS_ppm:", opt.ppm), file= f.info, append=TRUE,sep="\n")
cat( paste("XCMS_peakwidth:", opt.pw), file= f.info, append=TRUE,sep="\n")
cat( paste("absMinIntensity:", min.intensity.thr), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("relMinReplicateAbundance :", rp.feature), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("blankThreshold:", bk.sa.thr), file= f.info, append=TRUE,sep="\n")
if (MD.filter =="YES"){
  cat("Used suspect list Mass defect filtering ",
#cat(paste("Mass defect filtering from", MD.minmax[1],"to", MD.minmax[2]),
    file= f.info,append=TRUE, sep="\n")
  }else{}
cat( paste("adduct:", adduct), file= f.info, append=TRUE,sep="\n")
cat( paste("formula:", form.ele), file= f.info, append=TRUE,sep="\n")
cat( paste("suspect list:", fns), file= f.info, append=TRUE,sep="\n")
cat( paste("MS2 library:", MS2.lib), file= f.info, append=TRUE,sep="\n")
cat( paste("MetFrag list:", fn.metfrag), file= f.info, append=TRUE,sep="\n")

## load MS2 library
for (k in 1:length(MS2.lib))
{
  if (k==1) {
    mslibraryM <- loadMSLibrary(paste(PatRoon.dir,"/MS2_library/",MS2.lib[k] , sep=""), "msp") 
  }else{
    mslibraryM <- append(mslibraryM, 
                         list(loadMSLibrary(paste(PatRoon.dir,"/MS2_library/",MS2.lib[k] , sep=""), "msp") ))
  }
}

if (length(MS2.lib)>1) {mslibraryM <- Reduce(function(x, y) merge(x, y, all = FALSE), mslibraryM)} else {}

## read suspect list
dat <- read.csv(fns, header=TRUE) #open suspect list

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

df.fList <- as.data.table(fList)
df.fList <- na.omit(df.fList)

write.table(df.fList, file=paste(outpath,"/raw_unaligned_ungrouped.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

## check ISTD on unaligned
if (check.istd=="YES") {
  
df.istd <- read.csv(paste(workPath,"/input/",istd.list,sep=""),
                    sep=",",header=TRUE) #open istd list

istd <- data.frame(name = df.istd$name,
                   formula = df.istd$formula,
                   rt = df.isstd$rt,
                   stringsAsFactors = FALSE) 

fGroupsISTD <- screenSuspects(fList, istd, 
                              rtWindow = 60,
                              mzWindow = 0.005,
                              onlyHits = TRUE)

# ## export ISTD intensity data
df.fGroupsISTD <- as.data.table(fGroupsISTD)
df.fGroupsISTD <- na.omit(df.fGroupsISTD)

write.table(df.fGroups, file=paste(outpath,"/ISTD_check.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )
}else{}

# performed RT alignement and group feature
fGroups <- groupFeatures(fList, "xcms3")
                         
# export raw data for control
df.fGroups <- as.data.table(fGroups)
df.fGroups <- na.omit(df.fGroups)

write.table(df.fGroups, file=paste(outpath,"/raw_aligned_grouped.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )


# Basic rule based filtering
fGroups <- patRoon::filter(fGroups,  
                            absMinIntensity = min.intensity.thr, 
                            relMinReplicateAbundance = rp.feature, 
                            blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                            retentionRange = NULL, mzRange = NULL)

if (MD.filter=="YES") # Mass defect filtration
  {
  MD <- dat$MONOISOTOPIC_MASS-floor(dat$MONOISOTOPIC_MASS) # to follow patRoon def of mass defect.
  MD.minmax <- c(min(MD),max(MD)) 
  
  fGroups <- patRoon::filter(fGroups ,mzDefectRange = MD.minmax ) 
  
  }else{}

# -------------------------
# reporting
# -------------------------
## export groupfeature as table
df.fGroups <- as.data.table(fGroups)

write.table(df.fGroups, file=paste(outpath,"/featureGroups.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

## export averaged groupfeature as table
df.fGroups <- as.data.table(fGroups, average = TRUE)
                
write.table(df.fGroups, file=paste(outpath,"/featureGroups_averaged.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
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
                                    calculateFeatures = FALSE,
                                    featThresholdAnn = 1,
                                    MSMode ="both")
## MetFrag -----
################
# Calculate compound structure candidates
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
    
suspects <- data.frame(name = dat$ID,
                         SMILES =dat$SMILES,
                         stringsAsFactors = FALSE) 
  
fGroupsSusp <- screenSuspects(fGroups, suspects, 
                              mzWindow = 0.005,
                              onlyHits = TRUE) 
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

# Take the highest Level conf. between MF and Lib
fGroupsSusp_MF <- patRoon::filter(fGroupsSusp_MF,
                        selectHitsBy="level", 
                        selectBestFGroups=TRUE) 

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

# Select the first "BEST" occurrence of each group
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

write.table(df.fGroupsSusp, file=paste(outpath, "/SuspectScreening_all.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

# create a working table for sample only 
df.fGroupsSusp <- data.frame(df.fGroupsSusp)
names(df.fGroupsSusp) <- sub("^X", "", names(df.fGroupsSusp)) # rename if name start by X -- 
data <- df.fGroupsSusp[ ,names(df.fGroupsSusp) %in% unique(df$group[df$sampletype=="SA"]) ]

# select minimal info
df.data <- cbind(data.frame(cbind(group= df.fGroupsSusp$group,
                  ret= as.numeric(df.fGroupsSusp$ret),
                  mz= as.numeric(df.fGroupsSusp$mz),
                  InChIKey=df.fGroupsSusp$susp_InChIKey,
                  estIDLevel=df.fGroupsSusp$susp_estIDLevel,
                  IUPACName=df.fGroupsSusp$IUPACName,
                  LogP= as.numeric(df.fGroupsSusp$LogP) )), data)

# final check if all value = zero 
if (ncol(df.data)==8){
      df.data <- df.data[df.data[,8]!=0,]
    }else{
      index <- rowSums( df.data[,8:ncol(df.data)]) >0 
      df.data <- df.data[index,]
    }

write.table(df.data, file=paste(outpath, "/SuspectScreening_sample.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

