
########### DON'T use it for the moment please us
# v5.4 of the workflow

###############################################################################
## Title: PatRoon - suspect screening 
###############################################################################
## version:6.3
## Date: August 2026
## Author: Boris Droz 
## Modified from Tutorial and Handbook on https://github.com/rickhelmus/patRoon
## Depends:
##        R(>=4.6.1)
##        patRoon(>=3.1.0)
##        Java 21 download at https://www.oracle.com/java/technologies/downloads/#jdk21-windows
###############################################################################
## Description:
###############
## Performed suspect analysis with MS2 library search and insilico Metfrag
## 
## Should have a input file with a raw, mzxml folder and the sample list (csv)
##
## REQUIERE:
## YOU WILL NEED to update the suspect list, Metfrag Puchem DB and MS2 library 
## before starting
## just make sure to used the right R version and the R associated R-library
## 
################################################################################
## Parameter -- MODIFIED IF NEEDED
############
script.name <- "2.PatRoon_suspectXCMS_DDA_v6.3.R"
## path
workPath <- "D:/patRoon_NTS/BioCrab/"

## Input data - 
# sample.list <- "sample_list_midcal_mzML.csv"
sample.list <- "sample_biocrab_2025_B1.csv"

# check for ISTD - option are TRUE or FALSE 
## option doesn't work on this version yet.
check.istd <- FALSE
# istd.list <- "istd_list.csv"

## Optimized XCMS parameters for peak picking
mode <- c("pos") #,"neg") # should be in order pos then neg
opt.ppm.pos = 10
opt.pw.pos = c(5, 15) # peak width min and max for pos
opt.ppm.neg = 10
opt.pw.neg = c(5, 15) # peak width min and max for neg

## Parameter for filtering check patroon  help(filter) - NULL to remove filter
min.intensity.thr = 200 # remove features <intensity, typical range between 100 - 1000
rp.feature = 1 # 1= features should be in all analysis of replicate groups
max.Rep.RSD = 0.75 # remove features with intensity RSD in replicates >ratio - 0.75 good cutoff
bk.sa.thr = 3 # remove features <x intensity of (average) blank intensity - never go under 3
RT.range = c(0, Inf) #in sec 0 Inf to not select

# Mass defect filtering - option are TRUE or FALSE
MD.filter <- FALSE # used the suspect list to mass filtering.

## Adduct and formula search parameter
custom.metfrag <-TRUE #generate custom metfrag 
adduct <- c("[M+H]+") #,"[M-H]-")
form.search <- TRUE #search formula is very long sometimes not necessary 
form.ele <- "CHNOPSClFBr" #  NO Necessary NEED to modify it -- been re-filter based on you suspect list
ppm.error <- 10 ## qtof typically 5-10 - check you QA/QC to define if smaller

#post analysis
pubchem.search <-FALSE # Kow and add para.search in pubchem - take time
                      
###############################################################################
## Path of the script
####################
## path -- generally set once and do not need to be change
# PatRoon.directory -> where suspect, MS2, ... are
PatRoon.dir <- "D:/patRoon_NTS/library/"

## set suspect list, MS2 and Metfrag all located in PatRoon.dir
f.susp.list <- paste(PatRoon.dir,"/suspect_list/pos_susdat_2024-01-12-082351_ed.csv",sep="")
f.target.list <- paste(PatRoon.dir,"/suspect_list/", sep="")

## set path to -- GENERALLY DO NOT NEED TO MODIFY
options(patRoon.path.obabel = "C:/Program Files/OpenBabel-3.1.1/") # open babel exe
options(patRoon.path.MetFragCL = 
          paste(PatRoon.dir,"/MetFrag/MetFragCommandLine-2.5.0.jar", sep="")) # metfrag exe

################################################################################
## load library -- DO NOT MODIFIED
library(patRoon) # v2.3.0
library(xcms) 
library(BiocParallel) 
library(dplyr)
library(webchem)
library(stats)
library(Rcpp)
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
folder <- paste("/",date,"_SuspAnal_patRoon", sep="")
outpath <- creat.subDir(paste(workPath,"/output",sep=""), folder)
inpath <- paste(workPath,"input",sep="")

setwd(workPath) # set directory

#load suspect list
dat <- read.csv(f.susp.list, header=TRUE) #open suspect list
# d.target

# creat a formula element list if form.search
if (form.search)
  {
    all_formulas <- paste(dat$FORMULA, collapse = "") # Extract and concatenate all formulas
    elements <- sort( unique(unlist(regmatches(all_formulas, gregexpr("[A-Z][a-z]?", all_formulas)))) ) # Extract unique element
    
    form.ele <- unlist(regmatches(form.ele, gregexpr("[A-Z][a-z]?", form.ele))) # create character list
    
    filtered_list <- na.omit(form.ele[form.ele %in% elements]) #re-filter based on selected element
    form.ele <- paste(filtered_list, collapse = "") 
  }

if (custom.metfrag) # generate kustome Metfrag database
    { # make customized compound database for MetFrag bases on suspect list
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
    }
################################################################################
## load MS2 library
###################
MS2.lib <- list.files(paste(PatRoon.dir, "/MS2_library",sep=""),pattern= "msp" )#load all library contained here... 
# check if compilated version exist
check.dir <- length(list.dirs(paste(PatRoon.dir, "/MS2_library",sep=""),recursive = FALSE ))==0

if (check.dir){
 cat("COMPIL MS2 library is running...")
  
  lib.summary <- data.frame(Library = character(),N_Spectra = integer(),  stringsAsFactors = FALSE)
    
  for (k in 1:length(MS2.lib))
    {
    lib <- loadMSLibrary(paste(PatRoon.dir, "/MS2_library/", MS2.lib[k], sep = ""), "msp")
    
    lib.summary <- rbind(lib.summary,
                    data.frame( Library = basename(MS2.lib[k]),
                                N_Spectra = length(lib), # count spectra
                                stringsAsFactors = FALSE) )
  #   if (k==1) {
  #     mslibraryM <- loadMSLibrary(paste(PatRoon.dir,"/MS2_library/",MS2.lib[k] , sep=""), "msp")
  # 
  #   }else{
  #     mslibraryM <- append(mslibraryM,
  #                          list(loadMSLibrary(paste(PatRoon.dir,"/MS2_library/",MS2.lib[k] , sep=""), "msp") ))
  #   }
  # }
  
  if (k == 1) { mslibraryM <- lib
      } else {
        mslibraryM <- append(mslibraryM, list(lib))
      }
    }

  if (length(MS2.lib)>1) {mslibraryM <- Reduce(function(x, y) merge(x, y, all = FALSE), mslibraryM)} else {}
  
  # total
  lib.summary <- rbind(lib.summary,
              data.frame(Library = "TOTAL compiled spectra",
                        N_Spectra = length(mslibraryM),
                      stringsAsFactors = FALSE) )

  folder <- paste("/","compil_MS2lib", sep="")
  compil.lib.path <- creat.subDir(paste(PatRoon.dir,"/MS2_library",sep=""), folder)

  patRoon::export(mslibraryM, type="msp",
                  out = paste(compil.lib.path,"/",date,"_compil_MS2lib.msp",sep="") )
  
  write.csv( lib.summary,
            file = paste(compil.lib.path, "/MS2_library_summary.csv", sep = ""),
          row.names = FALSE )

  saveRDS(MS2.lib, file=paste(compil.lib.path,"/","MS2.lib.rds",sep="") )

} else {

  # check if the compil we have all MS2.lib

  old.MS2.lib <- readRDS(file=paste(compil.lib.path,"/","MS2.lib.rds",sep=""))

  lib.changed <- !identical(old.MS2.lib, MS2.lib)

  if(lib.changed){ #compil again
    cat("Library as been changed COMPIL MS2 is running again...")
    
    lib.summary <- data.frame(Library = character(),N_Spectra = integer(),  stringsAsFactors = FALSE)

    for (k in 1:length(MS2.lib))
    {
      lib <- loadMSLibrary(paste(PatRoon.dir, "/MS2_library/", MS2.lib[k], sep = ""), "msp")
      
      lib.summary <- rbind(lib.summary,
                           data.frame( Library = basename(MS2.lib[k]),
                                       N_Spectra = length(lib), # count spectra
                                       stringsAsFactors = FALSE)  )
      # if (k==1) {
      #   mslibraryM <- loadMSLibrary(paste(PatRoon.dir,"/MS2_library/",MS2.lib[k] , sep=""), "msp")
      # 
      # }else{
      #   mslibraryM <- append(mslibraryM,
      #                        list(loadMSLibrary(paste(PatRoon.dir,"/MS2_library/",MS2.lib[k] , sep=""), "msp") ))
      # }
      
      if (k == 1) { mslibraryM <- lib
                } else {
                  mslibraryM <- append(mslibraryM, list(lib))
                }
    }

    if (length(MS2.lib)>1) {mslibraryM <- Reduce(function(x, y) merge(x, y, all = FALSE), mslibraryM)} else {}
    
    # total
    lib.summary <- rbind(lib.summary,
                         data.frame(Library = "TOTAL compiled spectra",
                                    N_Spectra = length(mslibraryM),
                                    stringsAsFactors = FALSE)  )

    folder <- paste("/","compil_MS2lib", sep="")
    #compil.lib.path <- creat.subDir(paste(PatRoon.dir,"/MS2_library",sep=""), folder)

    patRoon::export(mslibraryM, type="msp",
                    out = paste(compil.lib.path,"/",date,"_compil_MS2lib.msp",sep="") )
    
    write.csv( lib.summary,
               file = paste(compil.lib.path, "/MS2_library_summary.csv", sep = ""),
               row.names = FALSE )

    saveRDS(MS2.lib, file=paste(compil.lib.path,"/","MS2.lib.rds",sep="") )
  }else{MS2.lib <-old.MS2.lib}
}

###############################
# save parameter of the script
###############################
f.info <- paste(outpath,"/AA_INFO_RUN_README.txt",sep="")
cat( paste("*** patRoon parameter for the run....", Sys.Date()), 
     file= f.info, append=TRUE, sep="\n")
cat( "#########################################################", 
     file= f.info, append=TRUE,sep="\n")
cat( paste("Scriptname: ", script.name,sep=""), 
     file= f.info, append=TRUE,sep="\n")
cat( paste("SampleList: ", workPath,"/input/",sample.list,sep=""), 
     file= f.info, append=TRUE,sep="\n")
if (check.istd){
  cat( paste("ISTDList: ", workPath,"/input/",istd.list,sep=""), 
       file= f.info, append=TRUE,sep="\n")
  }
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
if (MD.filter){
  cat("Used suspect list Mass Defect filtering ",
    file= f.info,append=TRUE, sep="\n")
    }
cat( "Annotation setting", 
     file= f.info, append=TRUE, sep="\n")
cat( "#########################################################", 
     file= f.info, append=TRUE,sep="\n")
if (custom.metfrag){
  cat( "MetFrag list: custom_metfrag_db.csv", file= f.info, append=TRUE,sep="\n")
}else{
  cat( "MetFrag list: Pubchem", file= f.info, append=TRUE,sep="\n")
}
cat( paste("adduct:", adduct), file= f.info, append=TRUE,sep="\n")
cat( paste("formula search:", form.search), file= f.info, append=TRUE,sep="\n")
if (form.search)
  {cat( paste("formula:", form.ele), file= f.info, append=TRUE,sep="\n")}
cat( paste("ppm error:", ppm.error), file= f.info, append=TRUE,sep="\n")
cat( paste("suspect list:", f.susp.list), file= f.info, append=TRUE,sep="\n")
cat( paste("MS2 library:", MS2.lib), file= f.info, append=TRUE,sep="\n")

if (pubchem.search){
  cat( "Pubchem post annotation is done to search Kow", file= f.info, append=TRUE,sep="\n")
}



## load data info
# df <- read.csv(paste(workPath,"/input/",sample.list,sep=""),
#                 sep=",",header=TRUE)
# 
# df <-df[df$folder=="C18_pH7_pos"|df$folder=="C18_pH7_neg",]
# 
# anaInfo <- data.frame(cbind(path = #df$path, 
#                               paste(df$path_workdir,"/",df$folder,sep=""),
#                                 analysis =df$filename,
#                                 group = df$group,
#                                 blank = df$blank,
#                             mode=df$mode))
##########################
##test case
#############

#anaInfo <- anaInfo[anaInfo$group=="NIST_B"|anaInfo$group=="BK",]

### from patRoon
# get analysis information for example data, which can be used directly for feature detection
df <-cbind(patRoonData::exampleAnalysisInfo("positive")[, c("analysis", "replicate", "blank", "path_centroid")],
                mode="pos",
                sampletype=c("BK","BK","BK","SA","SA","SA") )

anaInfo <- df

# anaInfo <- rbind(anaInfo,
#                  cbind(patRoonData::exampleAnalysisInfo("negative")[, c("analysis", "replicate", "blank", "path_centroid")],
#                        mode="neg") )

################
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
      }

if (any(grepl("neg", mode)))
      {
        param.xcms.neg <- xcms::CentWaveParam(ppm = opt.ppm.neg,
                                              peakwidth = opt.pw.neg,
                                              snthresh = 10,
                                              prefilter = c(3, 100),
                                              noise = 50 )
        
        fListNeg <- findFeatures(anaInfo[anaInfo$mode=="neg",], "xcms3", param = param.xcms.neg)
      }

## merged pos and neg group together
has_pos <- any(grepl("pos", mode))
has_neg <- any(grepl("neg", mode))

if (has_pos && has_neg) {
  fList <- makeSet(fListPos, fListNeg, adducts =  adduct)
  } else if (has_pos) {
    fList <- makeSet(fListPos, adducts =  adduct)
  } else if (has_neg) {
    fList <- makeSet(fListNeg, adducts =  adduct)
  }

# performed RT alignement and group feature
fGroups <- groupFeatures(fList, "xcms3")

df.fGroups <- as.data.table(fGroups, areas = TRUE) # count all feature
feat.summ <-data.frame(nb.all_grp.feat =nrow(df.fGroups)) 
if (has_pos) {
  feat.summ <-data.frame(c(feat.summ, nb.pos_grp.feat =sum(!is.na(df.fGroups$`adduct-positive`)) ) ) }
if (has_neg) {
  feat.summ <-data.frame(c(feat.summ, nb.neg_grp.feat =sum(!is.na(df.fGroups$`adduct-negative`)) ) ) }

# Basic rule based filtering
fGroups <- patRoon::filter(fGroups,  
                          absMinIntensity = min.intensity.thr, 
                          relMinReplicateAbundance = rp.feature,
                          maxReplicateIntRSD =max.Rep.RSD,
                          blankThreshold = bk.sa.thr, removeBlanks = TRUE,
                          retentionRange = RT.range, mzRange = NULL)

if (MD.filter==TRUE) # Mass defect filtration
    {
      MD <- dat$MONOISOTOPIC_MASS-floor(dat$MONOISOTOPIC_MASS) # to follow patRoon def of mass defect.
      fGroups <- patRoon::filter(fGroups, mzDefectRange = c(min(MD),max(MD)),
                                    negate = TRUE)
    }
# -------------------------
# reporting
# -------------------------
df.fGroups <- as.data.table(fGroups, areas = TRUE) # count all feature
feat.summ <-data.frame(c(feat.summ, nb.filtered_feat =nrow(df.fGroups)) )
if (has_pos) {
  feat.summ <-data.frame(c(feat.summ, nb.pos.filtered_feat =sum(!is.na(df.fGroups$`adduct-positive`)) ) ) }
if (has_neg) {
  feat.summ <-data.frame(c(feat.summ, nb.neg.filtered_feat =sum(!is.na(df.fGroups$`adduct-negative`)) ) ) }

write.table(df.fGroups, file=paste(outpath,"/featureGroups.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

## export averaged groupfeature as table
df.fGroups <- as.data.table(fGroups, average = TRUE,areas = TRUE)

write.table(df.fGroups, file=paste(outpath,"/featureGroups_averaged.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )
# -------------------------
# Annotation 
# -------------------------
## to speed up keep on sample here....
SA.group <- unique(df$replicate[df$sampletype=="SA"])
fGroups <- patRoon::filter(fGroups, replicate = SA.group)

# Retrieve MS peak lists 
mslists <- generateMSPeakLists(fGroups) 

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

MS2.feature <- df.fGroups[df.fGroups$type=="MSMS",]
feat.summ <-data.frame(c(feat.summ, nb.feat.MS2 =nrow(MS2.feature)) )
if (has_pos) {
  feat.summ <-data.frame(c(feat.summ, nb.pos.feat.MS2 =sum(!is.na(MS2.feature$`adduct-positive`)) ) ) }
if (has_neg) {
  feat.summ <-data.frame(c(feat.summ, nb.neg.feat.MS2 =sum(!is.na(MS2.feature$`adduct-negative`)) ) ) }

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
                              mzWindow = 0.010,
                              onlyHits = TRUE) 
if (form.search==TRUE) { 

  formulas <- generateFormulasGenForm(fGroupsSusp, 
                                      mslists,
                                    elements = form.ele,
                                    relMzDev = ppm.error,
                                    absAlignMzDev = 0.005,
                                    topMost = 10,
                                    calculateFeatures = FALSE,
                                    featThresholdAnn = 0,
                                    MSMode ="both")
  }else {}

## MetFrag -----
################
if (custom.metfrag){
      options(patRoon.path.MetFragPubChemLite = paste(outpath,"/custom_metfrag_db.csv",sep="") ) 
      MF.database= "pubchemlite"}else{ MF.database= "pubchem"}

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
              fragAbsMzDev = 0.001,
              database = MF.database,
              setThresholdAnn =0)

# df.fGroupsSusp_MF <-as.data.table(fGroupsSusp_MF, average = TRUE, areas = TRUE,
#                                   collapseSuspects = NULL,
#                                   onlyHits = TRUE)

# Annotation with the Library MS2 algorithm
#########################################
compsLib <- generateCompounds(fGroups, mslists, "library", 
                              MSLibrary = mslibraryM, minSim = 0.2)
# -------------------------
# Final Annotation 
# -------------------------

if (form.search) { 

fGroupsSusp_MF <- estimateIDConfidence(fGroupsSusp,
                                         MSPeakLists = mslists,
                                         formulas = formulas,
                                         compounds = compsMF)#

df.fGroupsSusp_MF <-as.data.table(fGroupsSusp_MF, average = TRUE, areas = TRUE,
                                  collapseSuspects = NULL,
                                  onlyHits = TRUE)

fGroupsSusp_Lib <- estimateIDConfidence(fGroupsSusp,
                    MSPeakLists = mslists,
                    formulas = formulas,
                    compounds = compsLib)

df.fGroupsSusp_Lib <-as.data.table(fGroupsSusp_Lib, average = TRUE, areas = TRUE,
                                  collapseSuspects = NULL,
                                  onlyHits = TRUE)

}else{
  fGroupsSusp_MF <- estimateIDConfidence(fGroupsSusp,
                    MSPeakLists = mslists,
                    compounds = compsMF)
  
  fGroupsSusp_Lib <- estimateIDConfidence(fGroupsSusp,
                      MSPeakLists = mslists,
                      compounds = compsLib)
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


if (pubchem.search==TRUE){
  
    dmol <- NULL
    for (i in seq_len(nrow(df.fGroupsSusp))) {
      get.id <- df.fGroupsSusp$susp_InChIKey[i]
      
      message(sprintf("Processing %d/%d : %s",
                      i, nrow(df.fGroupsSusp), get.id))
      #------------------------------
      # Get CID
      #------------------------------
      fcid <- tryCatch(
        get_cid(get.id, from = "inchikey", match = "first"),
        error = function(e) {
          Sys.sleep(5)
          tryCatch(
            get_cid(get.id, from = "inchikey", match = "first"),
            error = function(e) NULL
              ) } )
      
      # Skip if CID retrieval failed
      if (is.null(fcid) ||
          is.null(fcid$cid) ||
          length(fcid$cid) != 1 ||
          fcid$cid == 0 ||
          is.na(fcid$cid)) {
        
        dmol <- rbind(dmol,
                      data.frame(CID = NA,
                                 IUPACName = NA,
                                 XLogP = NA,
                                 CAS = NA))
        next
      }
      
      #------------------------------
      # Get properties
      #------------------------------
      prop <- tryCatch(
        pc_prop(fcid$cid,
                properties = c("IUPACName", "XLogP"),
                verbose = TRUE),
        error = function(e) {
          message(sprintf(
            "ERROR: pc_prop failed for CID %s (%s)",
            fcid$cid, e$message))
          NULL
        }
      )
      
      #------------------------------
      # Get CAS
      #------------------------------
      cas <- tryCatch(
        cir_query(fcid$query,
                  representation = "cas",
                  match = "first"),
        error = function(e) {
          message(sprintf(
            "ERROR: CAS lookup failed for %s (%s)",
            get.id, e$message))
          NULL
        }
      )
      
      # If property retrieval failed
      if (is.null(prop)) {
        dmol <- rbind(dmol,
                      data.frame(CID = fcid$cid,
                                 IUPACName = NA,
                                 XLogP = NA,
                                 CAS = if (!is.null(cas)) cas$cas else NA))
        next
      }
      
      dmol <- rbind(
        dmol,
        data.frame(
          CID = if ("CID" %in% names(prop)) prop$CID else NA,
          IUPACName = if ("IUPACName" %in% names(prop)) prop$IUPACName else NA,
          XLogP = if ("XLogP" %in% names(prop)) prop$XLogP else NA,
          CAS = if (!is.null(cas) && "cas" %in% names(cas)) cas$cas else NA
        )
      )
    }
    df.fGroupsSusp <- cbind(df.fGroupsSusp,dmol )
  }

df.fGroupsSusp <- df.fGroupsSusp %>%
  left_join(df.mslists, by = "group") %>%
  mutate(type = coalesce(type, "MS"))




## Annotation summary
MS2.annot <- df.fGroupsSusp[df.fGroupsSusp$type=="MSMS",]
feat.summ <-data.frame(c(feat.summ, nb.Annot.MS2 =nrow(MS2.annot) )) 
if (has_pos) {
  feat.summ <-data.frame(c(feat.summ, nb.pos.Annot.MS2 =sum(!is.na(MS2.annot$`adduct-positive`)) ) ) }
if (has_neg) {
  feat.summ <-data.frame(c(feat.summ, nb.neg.Annot.MS2 =sum(!is.na(MS2.annot$`adduct-negative`)) ) ) }

feat.summ <-data.frame(c(feat.summ, nb.Annot =nrow(df.fGroupsSusp) )) 
if (has_pos) {
  feat.summ <-data.frame(c(feat.summ, nb.pos.Annot =sum(!is.na(df.fGroupsSusp$`adduct-positive`)) ) ) }
if (has_neg) {
  feat.summ <-data.frame(c(feat.summ, nb.neg.Annot =sum(!is.na(df.fGroupsSusp$`adduct-negative`)) ) ) }


write.table(df.fGroupsSusp, file=paste(outpath, "/SuspectScreening_all.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

## save feat.summ
feat.summ <- stack(feat.summ)

write.table(feat.summ, file=paste(outpath, "/feature_annotation_summary.txt", sep=""),
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = TRUE )

## save.image(file=paste(outpath,'/',date,'_suspect_anal_session.RData',sep="") )

