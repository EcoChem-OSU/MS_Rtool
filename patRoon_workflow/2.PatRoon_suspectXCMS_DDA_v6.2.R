###############################################################################
## Title: PatRoon - suspect screening DON'T use it for the moment please us
# v5.4 of the workflow
###############################################################################
## version:6.2
## Date: August 2026
## Author: Boris Droz 
## Modified from Tutorial and Handbook on https://github.com/rickhelmus/patRoon
## Depends:
##        R(>=4.5.3)
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
script.name <- "2.PatRoon_suspectXCMS_DDA_v6.1.R"
## path
workPath <- "E:/Patroon_NTS/BioCrab/"

## Input data - 
# sample.list <- "sample_list_midcal_mzML.csv"
sample.list <- "sample_biocrab_2025_B1.csv"

# check for ISTD - option are "YES" or "NO"
check.istd <- "NO"
istd.list <- "istd_list.csv"

## Optimized XCMS parameters for peak picking
mode <- c("pos","neg") # should be in order pos then neg
opt.ppm.pos = 14
opt.pw.pos = c(7, 63) # peak width min and max
opt.ppm.neg = 14
opt.pw.neg = c(7, 63) # peak width min and max

## Parameter for filtering check patroon  help(filter)
min.intensity.thr = 200 # remove features <intensity, typical range between 100 - 1000
rp.feature = 1 # features should be in all analysis of replicate groups
max.Rep.RSD = 0.75 # remove features with intensity RSD in replicates >ratio - 0.75 good cutoff
bk.sa.thr = 3 # remove features <x intensity of (average) blank intensity - never go under 3
RT.range = c(0.2, 10*60) #in sec

# Mass defect filtering - option are "YES" or "NO"
MD.filter <- "NO" # used the suspect list to mass filtering.

## Abduct and formula search parameter
adduct <- c("[M+H]+","[M-H]-")
form.search <- "NO" #search formula is very long sometimes not necessary 
form.ele <- "CHNOPSClFBr" #  NO Necessary NEED to modify it -- been re-filter based on you suspect list
ppm.error <- 10 ## qtof typically 5-10 - check you QC

pubchem.search <-"NO" # Kow and add para.search in pubchem
                      
################################################################################
## load library -- DO NOT MODIFIED
library(patRoon) # v2.3.0
library(xcms) 
library(BiocParallel) 
library(dplyr)
library(webchem)
library(stats)
###############################################################################
## Parameter
############
## path -- generally set once and do not need to be change
# PatRoon.directory -> where suspect, MS2, ... are
PatRoon.dir <- "E:/Patroon_NTS/library/"

## set suspect list, MS2 and Metfrag all located in PatRoon.dir
fns <- paste(PatRoon.dir,"/suspect_list/pos_SLPFAS_nist_susdat_ocde_master_struct_v2.1.csv",sep="")

MS2.lib <- c("2025-06-04_JF_OSU.msp","DIMSpecForPFAS_2023-10-03.msp","Fluoros_2.5_ed.msp",
             "MassBank_NIST_2024.06.msp", "MassBank_RIKEN_2024.06.msp") 

## set path to -- GENERALLY DO NOT NEED TO MODIFY
options(patRoon.path.obabel = "C:/Program Files/OpenBabel-3.1.1/") # open babel exe
options(patRoon.path.MetFragCL = 
          paste(PatRoon.dir,"/MetFrag/MetFragCommandLine-2.5.0.jar", sep="")) # metfrag exe

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

#load suspect list
dat <- read.csv(fns, header=TRUE) #open suspect list

# creat a formula element list if form.search
if (form.search=="YES")
  {
    all_formulas <- paste(dat$FORMULA, collapse = "") # Extract and concatenate all formulas
    elements <- sort( unique(unlist(regmatches(all_formulas, gregexpr("[A-Z][a-z]?", all_formulas)))) ) # Extract unique element
    
    form.ele <- unlist(regmatches(form.ele, gregexpr("[A-Z][a-z]?", form.ele))) # create character list
    
    filtered_list <- na.omit(form.ele[form.ele %in% elements]) #re-filter based on selected element
    form.ele <- paste(filtered_list, collapse = "") 
  }else{}

# generate Metfrag database
# make customized compound database for MetFrag bases on suspect list
db <- data.frame(
      identifier = dat$ID,
      formula = dat$FORMULA,
      SMILES = dat$SMILES,
      InChIKey = dat$INCHIKEY,
      stringsAsFactors = FALSE)

# Remove empty or all-NA columns
db <- db[, colSums(is.na(db) | db == "") < nrow(db), drop = FALSE]

# Write to CSV
write.csv(db, file = paste(outpath,"/custom_metfrag_db.csv",sep=""), row.names = FALSE, quote = TRUE)

# save parameter of the script
f.info <- paste(outpath,"/AA_INFO_RUN_README.txt",sep="")
cat( paste("*** patRoon parameter for the run....", Sys.Date()), 
     file= f.info, append=TRUE, sep="\n")
cat( "#########################################################", 
     file= f.info, append=TRUE,sep="\n")
cat( paste("Scriptname: ", script.name,sep=""), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("SampleList: ", workPath,"/input/",sample.list,sep=""), 
     file= f.info, append=TRUE,sep="\n")
if (check.istd =="YES"){
  cat( paste("ISTDList: ", workPath,"/input/",istd.list,sep=""), 
       file= f.info, append=TRUE,sep="\n")
}else{}
cat( "Peak-picking parameter", 
     file= f.info, append=TRUE, sep="\n")
cat( "#########################################################", 
     file= f.info, append=TRUE,sep="\n")

if (any(grepl("pos", mode)))
    {
    cat( paste("XCMS_pos_ppm:", opt.ppm.pos), file= f.info, append=TRUE,sep="\n")
    cat( paste("XCMS_pos_peakwidth:", opt.pw.pos), file= f.info, append=TRUE,sep="\n")
    }

if (any(grepl("neg", mode)))
    {
      cat( paste("XCMS_neg_ppm:", opt.ppm.neg), file= f.info, append=TRUE,sep="\n")
      cat( paste("XCMS_neg_peakwidth:", opt.pw.neg), file= f.info, append=TRUE,sep="\n")
    }

cat( "Filtering parameter", 
     file= f.info, append=TRUE, sep="\n")
cat( "#########################################################", 
     file= f.info, append=TRUE,sep="\n")
cat( paste("absMinIntensity:", min.intensity.thr), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("relMinReplicateAbundance:", rp.feature), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("maxReplicateIntRSD:", max.Rep.RSD), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("blankThreshold:", bk.sa.thr), file= f.info, append=TRUE,sep="\n")
cat( paste("Retention time range (s):", RT.range[1], "to",RT.range[2] ),
     file= f.info, append=TRUE,sep="\n")
if (MD.filter =="YES"){
  cat("Used suspect list Mass defect filtering ",
    file= f.info,append=TRUE, sep="\n")
    }
cat( "Annotation setting", 
     file= f.info, append=TRUE, sep="\n")
cat( "#########################################################", 
     file= f.info, append=TRUE,sep="\n")
cat( paste("adduct:", adduct), file= f.info, append=TRUE,sep="\n")
cat( paste("formula search:", form.search), file= f.info, append=TRUE,sep="\n")
cat( paste("formula:", form.ele), file= f.info, append=TRUE,sep="\n")
cat( paste("ppm error:", ppm.error), file= f.info, append=TRUE,sep="\n")
cat( paste("suspect list:", fns), file= f.info, append=TRUE,sep="\n")
cat( paste("MS2 library:", MS2.lib), file= f.info, append=TRUE,sep="\n")
cat( "MetFrag list: custom_metfrag_db.csv", file= f.info, append=TRUE,sep="\n")
#cat( "MetFrag list: Pubchem", file= f.info, append=TRUE,sep="\n")
cat( paste("Pubchem para search:",pubchem.search), file= f.info, append=TRUE,sep="\n")

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


## load data info
df <- read.csv(paste(workPath,"/input/",sample.list,sep=""),
                sep=",",header=TRUE)

df <-df[df$folder=="C18_pH7_pos"|df$folder=="C18_pH7_neg",]

anaInfo <- data.frame(cbind(path = #df$path, 
                              paste(df$path_workdir,"/",df$folder,sep=""),
                                analysis =df$filename,
                                group = df$group,
                                blank = df$blank,
                            mode=df$mode))
##########################
##test case
##
anaInfo <- anaInfo[anaInfo$group=="NIST_B"|anaInfo$group=="BK",]

# -------------------------
# features
# -------------------------
# Find all features pos, neg or both

if (any(grepl("pos", mode)))
      {
      param.xcms.pos <- xcms::CentWaveParam(ppm = opt.ppm.pos,
                                    peakwidth = opt.pw.pos,
                                    snthresh = 10,
                                    prefilter = c(3, 100),
                                    noise = 50 )
      
      fListPos <- findFeatures(anaInfo[anaInfo$mode=="pos",], "xcms3", param = param.xcms.pos)
      
      # performed RT alignement and group feature
      fGroupsPos <- groupFeatures(fListPos, "xcms3")
      
      # Basic rule based filtering
      fGroupsPos <- patRoon::filter(fGroupsPos,  
                                 absMinIntensity = min.intensity.thr, 
                                 relMinReplicateAbundance = rp.feature,
                                 maxReplicateIntRSD =max.Rep.RSD,
                                 blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                                 retentionRange = RT.range, mzRange = NULL)
      
      if (MD.filter=="YES") # Mass defect filtration
      {
        MD <- dat$MONOISOTOPIC_MASS-floor(dat$MONOISOTOPIC_MASS) # to follow patRoon def of mass defect.
        fGroupsPos <- patRoon::filter(fGroupsPos, mzDefectRange = c(min(MD),max(MD)),
                                   negate = TRUE)
      }else{ }
      
    }

if (any(grepl("neg", mode)))
      {
        param.xcms.neg <- xcms::CentWaveParam(ppm = opt.ppm.neg,
                                              peakwidth = opt.pw.neg,
                                              snthresh = 10,
                                              prefilter = c(3, 100),
                                              noise = 50 )
        
        fListNeg <- findFeatures(anaInfo[anaInfo$mode=="neg",], "xcms3", param = param.xcms.neg)
        # performed RT alignement and group feature
        fGroupsNeg <- groupFeatures(fListNeg, "xcms3")
        
        # Basic rule based filtering
        fGroupsNeg <- patRoon::filter(fGroupsNeg,  
                                   absMinIntensity = min.intensity.thr, 
                                   relMinReplicateAbundance = rp.feature,
                                   maxReplicateIntRSD =max.Rep.RSD,
                                   blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                                   retentionRange = RT.range, mzRange = NULL)
        
        if (MD.filter=="YES") # Mass defect filtration
        {
          MD <- dat$MONOISOTOPIC_MASS-floor(dat$MONOISOTOPIC_MASS) # to follow patRoon def of mass defect.
          fGroupsNeg <- patRoon::filter(fGroupsNeg, mzDefectRange = c(min(MD),max(MD)),
                                     negate = TRUE)
        }else{ }
      }

# has_pos <- any(grepl("pos", mode))
# has_neg <- any(grepl("neg", mode))
# 
# if (has_pos && has_neg) {
#       fList <- makeSet(fListPos, fListNeg, adducts = adduct)
#     } else if (has_pos) {
#       fList <- makeSet(fListPos, adducts = adduct)
#     } else if (has_neg) {
#       fList <- makeSet(fListNeg, adducts = adduct)
#     }

# df.fList <- as.data.table(fList)
# df.fList <- na.omit(df.fList)
# 
# feat.summ <-data.frame(nb.r.unali_feat =nrow(df.fList))
# 
# write.table(df.fList, file=paste(outpath,"/raw_unaligned_ungrouped.txt", sep=""),
#             append = FALSE, quote = FALSE, sep = "\t",
#             row.names = FALSE,col.names = TRUE )

## check ISTD on unaligned
# if (check.istd=="YES") {
#   
# df.istd <- read.csv(paste(workPath,"/input/",istd.list,sep=""),
#                     sep=",",header=TRUE) #open istd list
# 
# istd <- data.frame(name = df.istd$name,
#                    formula = df.istd$formula,
#                    rt = df.isstd$rt,
#                    stringsAsFactors = FALSE) 
# 
# fGroupsISTD <- screenSuspects(fList, istd, 
#                               rtWindow = 60,
#                               mzWindow = 0.005,
#                               onlyHits = TRUE)
# 
# # ## export ISTD intensity data
# df.fGroupsISTD <- as.data.table(fGroupsISTD, areas = TRUE)
# df.fGroupsISTD <- na.omit(df.fGroupsISTD)
# 
# write.table(df.fGroups, file=paste(outpath,"/ISTD_check.txt", sep=""),
#             append = FALSE, quote = FALSE, sep = "\t",
#             row.names = FALSE,col.names = TRUE )
# }else{}


                         
# # export raw data for control
# df.fGroups <- as.data.table(fGroups, areas = TRUE)
# df.fGroups <- na.omit(df.fGroups)
# 
# feat.summ <-data.frame(c(feat.summ, nb.al_gfeat =nrow(df.fGroups)) )
# 
# write.table(df.fGroups, file=paste(outpath,"/raw_aligned_grouped.txt", sep=""),
#             append = FALSE, quote = FALSE, sep = "\t",
#             row.names = FALSE,col.names = TRUE )


#components <- generateComponents(fGroups, "camera")
#fGroups <- selectIons(fGroups, components, prefAdduct=c("[M+H]+", "[M-H]-"))

## merged pos and neg group together
has_pos <- any(grepl("pos", mode))
has_neg <- any(grepl("neg", mode))

if (has_pos && has_neg) {
  fGroups <- makeSet(fGroupsPos, fGroupsNeg, adducts = adduct)
    } else if (has_pos) {
      fGroups <- makeSet(fGroupsPos, adducts = adduct)
    } else if (has_neg) {
      fGroups <- makeSet(fGroupsNeg, adducts = adduct)
    }

# -------------------------
# reporting
# -------------------------
## export groupfeature as table
df.fGroups <- as.data.table(fGroups, areas = TRUE)

feat.summ <-data.frame(c(feat.summ, nb.fil_feat =nrow(df.fGroups)) )

write.table(df.fGroups, file=paste(outpath,"/featureGroups.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

## export averaged groupfeature as table
df.fGroups <- as.data.table(fGroups, average = TRUE,areas = TRUE)

# write.table(df.fGroups, file=paste(outpath,"/featureGroups_averaged.txt", sep=""),
#             append = FALSE, quote = FALSE, sep = "\t",
#             row.names = FALSE,col.names = TRUE )
# -------------------------
# Annotation 
# -------------------------
## to speed up keep on sample here....
SA.group <- unique(df$group[df$sampletype=="SA"])
fGroups <- patRoon::filter(fGroups, rGroups = SA.group)

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

# get info for groupfeature which one as MS2
df.mslists <- as.data.table(mslists)

df.mslists <- df.mslists %>%
              group_by(group) %>%
              summarise(
                type = if (any(type == "MSMS", na.rm = TRUE)) "MSMS" else "MS",
                .groups = "drop" )

df.fGroups <- df.fGroups %>%
              left_join(df.mslists, by = "group") %>%
              mutate(type = coalesce(type, "MS"))

feat.summ <-data.frame(c(feat.summ, nb.feat.MS2 =nrow(df.fGroups[df.fGroups$type=="MSMS",])) )

write.table(df.fGroups, file=paste(outpath,"/featureGroups_averaged.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

# -------------------------
## Suspect screening 
#-------------------------
suspects <- data.frame(name = dat$ID,
                       SMILES =dat$SMILES,
                       stringsAsFactors = FALSE) 

fGroupsSusp <- screenSuspects(fGroups, suspects, 
                              mzWindow = 0.005,
                              onlyHits = TRUE) 
if (form.search=="YES") { 

  formulas <- generateFormulasGenForm(fGroupsSusp, mslists,
                                    elements = form.ele,
                                   relMzDev = ppm.error,
                                  absAlignMzDev = 0.002,
                                 topMost = 10,
                                calculateFeatures = FALSE,
                               featThresholdAnn = 1,
                              MSMode ="both")
  }else {}

## MetFrag -----
################
options(patRoon.path.MetFragPubChemLite = paste(outpath,"/custom_metfrag_db.csv",sep="") )

# Calculate compound structure candidates
compsMF <- generateCompounds(
              fGroupsSusp, mslists,
              "metfrag",
              method = "CL",
              timeout = 300,
              maxCandidatesToStop = 2000,
              timeoutRetries = 20,
              errorRetries = 2,
              topMost = 100,
              dbRelMzDev = ppm.error ,
              fragRelMzDev = 5,
              fragAbsMzDev = 0.002,
              database = "pubchemlite",
              setThresholdAnn = 1,
              scoreTypes = c("individualMoNAScore", "fragScore", "score"),
              scoreWeights = 1)

# Summary of MetFrag Results in a a Single Table
# MFsummary <- as.data.table(compounds)
# outputSummary <- paste(outpath, "MFsummary.csv", sep = "/")
# write.csv(MFsummary, outputSummary)

# Annotation with the Library MS2 algorithm
#########################################
compsLib <- generateCompounds(fGroups, mslists, "library", 
                              MSLibrary = mslibraryM, minSim = 0.4)
# -------------------------
# Final Annotation 
# -------------------------

if (form.search=="YES") { 

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
                    compoundsNormalizeScores = "max")
}else{
  fGroupsSusp_MF <- annotateSuspects(
    fGroupsSusp,
    MSPeakLists = mslists,
    compounds = compsMF,
    absMzDev = 0.005,
    specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
    checkFragments = c("mz", "formula", "compound"),
    formulasNormalizeScores = "max",
    compoundsNormalizeScores = "max") #
  
  fGroupsSusp_Lib <- annotateSuspects(
    fGroupsSusp,
    MSPeakLists = mslists,
    compounds = compsLib,
    absMzDev = 0.005,
    specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
    checkFragments = c("mz", "formula", "compound"),
    formulasNormalizeScores = "max",
    compoundsNormalizeScores = "max")
}

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
df.fGroupsSusp_MF <-as.data.table(fGroupsSusp_MF, average = TRUE, areas = TRUE,
                               collapseSuspects = NULL,
                               onlyHits = TRUE) 

df.fGroupsSusp_Lib <-as.data.table(fGroupsSusp_Lib, average = TRUE, areas = TRUE,
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

if (pubchem.search=="YES"){
  
  # get chemical names based on InChIKey
  # check if present in pubchem
  dmol <-NULL
  for (i in 1:nrow(df.fGroupsSusp) )
  {
    get.id <- df.fGroupsSusp$susp_InChIKey[i]
    
    fcid <- try ( get_cid(get.id, from =  'inchikey', match='first', verbose = TRUE, arg = NULL) , silent = TRUE )
    
    # if present in pubchem get
    if ( length(fcid$cid)==2 | fcid$cid==0 | is.na(fcid$cid) )  { 
      dmol <- rbind(dmol, c(CID=NA, IUPACName=NA,XLogP=NA, CAS=NA))
    }else{
      # data properties from pubchem
      prop <- pc_prop(fcid$cid, properties =  c('IUPACName','XLogP')
                      , verbose = TRUE)
      cas <- cir_query(fcid$query, representation = "cas", 
                       match = "first")
      
      # check length of prop 
      if (length(prop)==3){dmol <- rbind(dmol,c(CID=prop$CID, IUPACName=prop$IUPACName, XLogP=prop$XLogP, CAS= cas$cas) )
      }else{  if (all(names(prop) == c("CID","XLogP") ))
      { dmol <- rbind(dmol,c(CID=prop$CID, IUPACName=NA , XLogP=prop$XLogP, CAS=cas$cas) )
      }else{ dmol <- rbind(dmol,c(CID=prop$CID, IUPACName=prop$IUPACName,XLogP=NA,CAS=cas$cas) ) }
      }
    }
  } 
  df.fGroupsSusp <- cbind(df.fGroupsSusp,dmol )
  
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
                                    LogP= as.numeric(df.fGroupsSusp$LogP),
                                    CAS=df.fGroupsSusp$CAS )), data)
  
  # final check if all value = zero 
  if (ncol(df.data)==9){
    df.data <- df.data[df.data[,9]!=0,]
  }else{
    index <- rowSums( df.data[,9:ncol(df.data)]) >0 
    df.data <- df.data[index,]
  }
  
  write.table(df.data, file=paste(outpath, "/SuspectScreening_sample.txt", sep=""),
              append = FALSE, quote = FALSE, sep = "\t",
              row.names = FALSE,col.names = TRUE )
}else{
  
  df.fGroupsSusp <- df.fGroupsSusp %>%
    left_join(df.mslists, by = "group") %>%
    mutate(type = coalesce(type, "MS"))
  
  write.table(df.fGroupsSusp, file=paste(outpath, "/SuspectScreening_all.txt", sep=""),
              append = FALSE, quote = FALSE, sep = "\t",
              row.names = FALSE,col.names = TRUE )
  
feat.summ <-data.frame(c(feat.summ, nb.Annot.MS2 =nrow(df.fGroupsSusp[df.fGroupsSusp$type=="MSMS",]) )) 
  
}

save.image(file=paste(outpath,'/',date,'_suspect_NTA_session_DDA.RData',sep="") )



