###############################################################################
## Title: PatRoon - TPs from a suspect list
###############################################################################
## version:1.1 - modified from Tutorial version of the PatRoon and Handbook
## Date: 2026-05-19
## Author: Boris Droz 
## Tutorial and Handbook on https://github.com/rickhelmus/patRoon
## Depends:
##        R(>=4.3.1)
##        patRoon(>=2.3.0)
##
################################################################################
## Description:
##  Creat a suspect list of transformation product for given target compounds
## minima need: a list of SMILES

################################################################################
## DEPENDENCY AND FIRST INSTALL:
################################################################################
####  --- NEED TO INSTALL ANACONDA or OTHER PYTHON VERSION FIRST ---------------
################################################################################
#### NEED TO INSTALLL THE FOLLOWING STUFF:
# install.packages("remotes")
# remotes::install_github("LarsRichter82/larspack")
## install_miniconda(path = "C:/Users/drozditb/AppData/Local/miniconda", update = TRUE)

library(larspack) #RDKit for R https://zenodo.org/records/10275225
library(reticulate)
library(stringr)
library (patRoon)
library(webchem)
library(rcdk)
library(rinchi)
library(dplyr)
# 
### test communication with python
# repl_python()
# reticulate::py_config()

##### if don't communicate try to gave the path of the version
# reticulate::use_condaenv(condaenv = "my-rdkit-env", 
#                           conda = "/Users/drozditb/AppData/Local/miniconda/condabin")

# reticulate::use_condaenv(condaenv = "my-rdkit-env", 
#                          conda = "auto")

py_install("rdkit") # install python library
rdkit <- import("rdkit.Chem")
################################################################################
###############################################################################
## Parameter
############

workdir <- "C:/Users/drozditb/Documents/GitHub/MS_Rtool/Suspect_list"
  
fns <- file.choose()  #list with SMILES and more
# fns <- "/Antibiotics_mz_v2.csv"

generation = 2  # how much generation of TPs to consider, default = 2, 
                # 3 is max for CTS and Biotransformer

#### TPS library see help if need to change
CTS.lib <- "combined_photolysis_abiotic_hydrolysis"
BioTRansformer.typ <- "env"
  
################################################################################
################################################################################
#                                   Start Script   
##     --- Don't modify after here just run the script --------
################################################################################
# NEW FUNCTION
# is_empty_or_na <- function(x) { is.na(x) || x == "" || trimws(x) == "" }
is_empty_or_na <- function(x) { 
  if (missing(x) || is.null(x)) return(TRUE)  # Missing or NULL -> TRUE
  
  if (is.data.frame(x)) return(ncol(x) == 0)  # Data frame with 0 cols -> TRUE
  
  if (is.atomic(x) && length(x) == 1 && is.na(x)) return(TRUE)  # Single NA -> TRUE
  
  if (is.character(x) && all(trimws(x) == "")) return(TRUE)  # Empty string or only whitespace -> TRUE
  
  return(FALSE)  # Everything else -> FALSE
}
################################################################################
setwd(workdir)
  
dat <- read.csv(paste(workdir,fns,sep=""), header =TRUE) 
dat <- dat[,names(dat)=="name"|names(dat)=="SMILES"] #select only what you need!
name.list <- tools::file_path_sans_ext(basename(fns))

# 1) fill up your table
######################
i=1
# 1a) unified and creat field if needed
##########################################
if(any( c("ID", "DTXSID", "SUSPECTID", "Norman_SusDat_ID", "DB") %in% names(dat) )) {
  names(dat[names(dat)=="ID" |names(dat)=="DTXSID" |names(dat)== "SUSPECTID"|names(dat)=="Norman_SusDat_ID"|names(dat)=="DB"]) <- "ID"
}else{
  dat <- cbind(dat, ID= paste(name.list[i],"_",seq(from=1,to=nrow(dat),by=1),sep="") )
}

if(any( "name" %in% names(dat) )) {
}else{
  dat <- cbind(dat, name= rep(NA,nrow(dat)) )
}

if(any( c("CAS","CAS_NUMBER") %in% names(dat) )) {
  names(dat[names(dat)=="CAS" | names(dat)=="CAS_NUMBER"]) <- "CAS"
  }else{
  dat <- cbind(dat, CAS= rep(NA,nrow(dat)) )
}

if(any( c("INCHIKEY", "StdInChIKey", "StdInChIKey") %in% names(dat) )) {
  names(dat[names(dat)=="INCHIKEY"| names(dat)=="StdInChIKey"| names(dat)=="StdInChIKey"]) <- "INCHIKEY"
}else{
  dat <- cbind(dat, INCHIKEY= rep(NA,nrow(dat)) )
}

if(any( c("MOLECULAR.FORMULA", "MOLECULAR_FORMULA", "Formula", "FORMULA") %in% names(dat) )) {
  names(dat[names(dat)=="MOLECULAR.FORMULA" | names(dat)=="MOLECULAR_FORMULA"| names(dat)=="Formula"| names(dat)=="FORMULA"]) <- "FORMULA"
}else{
  dat <- cbind(dat, FORMULA= rep(NA,nrow(dat)) )
}

if(any( c("MONOISOTOPIC.MASS", "MONOISOTOPIC_MASS", "Monoisotopic_mass", "FIXEDMASS") %in% names(dat) )) {
  names(dat[names(dat)=="MONOISOTOPIC.MASS" | names(dat)=="MONOISOTOPIC_MASS"| names(dat)=="Monoisotopic_mass"| names(dat)== "FIXEDMASS" ]) <- "MONOISOTOPIC_MASS"
}else{
  dat <- cbind(dat, MONOISOTOPIC_MASS= rep(NA,nrow(dat)) )
}

if(any( c("SMILES","SMILES_MS_ready") %in% names(dat) )) {
  names(dat[names(dat)=="SMILES"| names(dat)=="SMILES_MS_ready"]) <-"SMILES" 
}else{
  dat <- cbind(dat, SMILES= rep(NA,nrow(dat)) )
}

if(any( "parent" %in% names(dat) )) {
  dat$parent <- gsub(":", "_", dat$name)
}else{
  dat <- cbind(dat, parent= gsub(":", "_", dat$name) ) 
}

## creat column for futur data
if(any( "transformation" %in% names(dat) )) {
}else{
  dat <- cbind(dat, transformation= rep("no_transf",nrow(dat))) 
}

if(any( "generation" %in% names(dat) )) {
}else{
  dat <- cbind(dat, generation= rep(0,nrow(dat))) 
}

if(any( "Type" %in% names(dat) )) {
}else{
  dat <- cbind(dat, Type = rep("PC",nrow(dat))) 
}

if(any( "biosystem" %in% names(dat) )) {
}else{
  dat <- cbind(dat, biosystem = rep(NA,nrow(dat))) 
}

if(any( "REF" %in% names(dat) )) {
}else{
  dat <- cbind(dat, REF= rep(NA,nrow(dat))) 
}
################################################################################
# 1b) fill empty field and check data
##########################################
for (j in 1:nrow(dat))
{
  # check if empty line and search for info or compute to fill up...
  ## first fill up SMILES
  if (is_empty_or_na(dat$SMILES[j])) { 
    if (is_empty_or_na(dat$INCHIKEY[j])){
      if (is_empty_or_na(dat$CAS[j])) {extracted_smiles <-"NA" }else{
        cids <- get_cid(dat$CAS[j],from = "cas", match = "first" )
        if(is.na(cids$cid)){extracted_smiles <-"NA"  }else{
          extracted_smiles <- pc_prop(cids$cid, properties = "CanonicalSMILES")}
      }
    }else{
      cids <- get_cid(dat$INCHIKEY[j],from = "inchikey", match = "first" )
      extracted_smiles <- pc_prop(cids$cid, properties = "CanonicalSMILES")
    }
    if(length(extracted_smiles)==2 ) {
      dat$SMILES[j] <- extracted_smiles$CanonicalSMILES}
    
  }else{ # check if smiles is canonical
    mol <- rdkit$MolFromSmiles(dat$SMILES[j])
    dat$SMILES[j] <- rdkit$MolToSmiles(mol)
  } 
  
  ## second get inchikey if you had smile
  if(is_empty_or_na(dat$SMILES[j])){ }else{
    if (is_empty_or_na(dat$INCHIKEY[j])) { 
      dat$INCHIKEY[j] <- sapply(dat$SMILES[j], get.inchi.key)
    }
  }
  ## finally get CAS
  if (is_empty_or_na(dat$CAS[j])) { 
    extracted_cas <- tryCatch(
      {
        result <- cir_query(dat$SMILES[j],
                            "cas", match = "first")
        if (is.null(result)) NA else result  # Assign NA if result is NULL
      },
      error = function(e) {
        message("Error occurred: ", e$message)
        return(NA)  # Assign NA if an error occurs
      }
    )
    if(length(extracted_cas)==2 ) {
      dat$CAS[j] <- extracted_cas$cas}
  }
  # compute formula and monoisotopic mass based on smiles if found before
  if (is_empty_or_na(dat$SMILES[j])) { }else{
    dat$FORMULA[j] <- RChemMass::MolFormFromSmiles.rcdk(dat$SMILES[j])
    
    dat$MONOISOTOPIC_MASS[j] <- MetaboCoreUtils::calculateMass(dat$FORMULA[j])
  }
  
  if (is_empty_or_na(dat$name[j])) { 
    cids <- get_cid(dat$SMILES[j],from = "smiles", match = "first" )
    dat$name[j] <- pc_prop(cids$cid, properties = "Title")$Title
    dat$name[j] <- gsub(":", "_", dat$name[j])
    dat$parent <- dat$name
    }else{ }
}
################################################################################
## Predict transformation Products
####################################
# Obtain transformation products (TPs) with Chemical Transformation Simulator (CTS)
## need to have one per pathway
TPs.CTS <- generateTPsCTS(
  dat,
  CTS.lib,
  generations = generation,
  errorRetries = 5,
  skipInvalid = TRUE,
  prefCalcChemProps = TRUE,
  neutralChemProps = TRUE,
  neutralizeTPs = TRUE,
  calcLogP = "none",
  calcSims = FALSE,
  parallel = TRUE)

TPs.CTS <- as.data.table(TPs.CTS)

dat <-rbind(dat, cbind(name=TPs.CTS$name, SMILES=TPs.CTS$SMILES, ID= TPs.CTS$ID, 
                       CAS=rep(NA,nrow(TPs.CTS)),
                       INCHIKEY=TPs.CTS$InChIKey, FORMULA=TPs.CTS$formula, 
                         MONOISOTOPIC_MASS=TPs.CTS$neutralMass, 
                       parent=TPs.CTS$parent, transformation=TPs.CTS$transformation,
                       generation=TPs.CTS$generation, Type=rep("TP",nrow(TPs.CTS) ), 
                       biosystem=rep(NA,nrow(TPs.CTS) ), 
                       REF=rep(paste("CTS",CTS.lib, sep="_"),nrow(TPs.CTS) ) ) )

### Obtain transformation products (TPs) with BioTransformer
TPs.biotransf <- generateTPsBioTransformer(
        dat,
        type = BioTRansformer.typ,
        generations = generation,
        maxExpGenerations = 2,
        extraOpts = NULL,
        skipInvalid = TRUE,
        prefCalcChemProps = TRUE,
        neutralChemProps = TRUE,
        neutralizeTPs = TRUE,
        calcSims = FALSE,
        MP = TRUE)

TPs.biotransf <- as.data.table(TPs.biotransf)

dat <-rbind(dat, cbind(name=TPs.biotransf$name, SMILES=TPs.biotransf$SMILES,ID= TPs.biotransf$ID, 
                       CAS=rep(NA,nrow(TPs.biotransf)), 
                       INCHIKEY=TPs.biotransf$InChIKey, FORMULA=TPs.biotransf$formula, 
                      MONOISOTOPIC_MASS=TPs.biotransf$neutralMass, 
                       parent=TPs.biotransf$parent, transformation=TPs.biotransf$transformation,
                       generation=TPs.biotransf$generation, Type=rep("TP",nrow(TPs.biotransf) ), 
                       biosystem=TPs.biotransf$biosystem, 
                       REF=rep(paste("BioTransformer",BioTRansformer.typ, sep="_"),nrow(TPs.biotransf) ) ) )

## Obtain transformation products (TPs) from a library
TPs.lib <- generateTPsLibrary(
        parents = dat,
        TPLibrary = NULL, #Pubchem as default
        generations = generation,
        skipInvalid = TRUE,
        prefCalcChemProps = TRUE,
        neutralChemProps = TRUE,
        neutralizeTPs = TRUE,
        matchParentsBy = "InChIKey",
        matchGenerationsBy = "InChIKey",
        calcSims = FALSE,
        fpType = "extended",
        fpSimMethod = "tanimoto"
        )

TPs.lib <- as.data.table(TPs.lib)

dat <-rbind(dat, cbind(name=TPs.lib$name_lib,  SMILES=TPs.lib$SMILES, ID= TPs.lib$ID, 
                       CAS=rep(NA,nrow(TPs.lib)),
                       INCHIKEY=TPs.lib$InChIKey, FORMULA=TPs.lib$formula, 
                       MONOISOTOPIC_MASS=TPs.lib$neutralMass, 
                       parent=TPs.lib$parent, transformation=TPs.lib$transformation,
                       generation=TPs.lib$generation, Type=rep("TP",nrow(TPs.lib) ), 
                       biosystem=TPs.lib$biosystem, 
                       REF= TPs.lib$evidencedoi ) )
################################################################################
################################################################################
# check out data result
#################

combine_duplicates_by_inchikey <- function(dat, key_col = "INCHIKEY", 
                                           combine_cols = c("name", "parent", 
                                              "transformation", "generation", "biosystem", "REF")) {
  # Ensure required packages are available
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install the 'dplyr' package.")
  
  # Combine specified columns for duplicates
  combined <- dat %>%
    group_by(.data[[key_col]]) %>%
    summarise(across(all_of(combine_cols), ~ paste(unique(.), collapse = "; "), .names = "combined_{col}"), .groups = "drop")
  
  # Get first occurrence of other columns
  other_cols <- setdiff(names(dat), c(combine_cols, key_col))
  first_occurrence <- dat %>%
    distinct(.data[[key_col]], .keep_all = TRUE) %>%
    select(all_of(c(key_col, other_cols)))
  
  # Merge combined columns with first occurrence
  final <- left_join(first_occurrence, combined, by = key_col)
  
  # Optionally rename combined columns back to original names
  names(final)[names(final) %in% paste0("combined_", combine_cols)] <- combine_cols
  
  return(final)
}

dat <- combine_duplicates_by_inchikey(dat)

# clean column "parent"
# Get the list of valid parent names from rows where Type == "PC"
valid_parents <- unique(trimws(unlist(strsplit(dat$parent[dat$Type == "PC"], ";"))))

# Function to clean each parent entry
clean_parent <- function(entry, valid_list) {
  parts <- trimws(unlist(strsplit(entry, ";")))
  filtered <- parts[parts %in% valid_list]
  paste(filtered, collapse = "; ")
}

# Apply the cleaning function to the whole parent column
dat$parent <- sapply(dat$parent, clean_parent, valid_list = valid_parents)

#' #You need to set anaconda environment prior to running RDKit code
# reticulate::use_condaenv(condaenv = "my-rdkit-env", conda = "/Users/lrichter/miniconda3/bin/conda")

## check for negatif
dat <- cbind(dat, NumHDonors = rep(NA,nrow(dat)) )

for (i in 1:length(dat$SMILES)) { # you might need to run this a couple of time..
  # for (i in 4741:length(PFAS.list.nist$SMILES)) { 
  tryCatch( {dat$NumHDonors[i] <-rdkit_descriptor(dat$SMILES[i], 
                                                          descriptor ="NumHDonors")  
  } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

## check for positif
dat <- cbind(dat, NumHAcceptors = rep(NA,nrow(dat)) )

for (i in 1:length(dat$SMILES)) { # you might need to run this a couple of time..
  # for (i in 4741:length(PFAS.list.nist$SMILES)) { 
  tryCatch( {dat$NumHAcceptors[i] <-rdkit_descriptor(dat$SMILES[i], 
                                                             descriptor ="NumHAcceptors")  
  } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#######
date <- Sys.Date()

# save the full list
write.csv(dat, paste(workdir,"/",date,"_TPs_suspectList.csv",sep=""),
          quote = TRUE,
          row.names = FALSE)

##############################################################################
# creat neg suspect list

## double check
unique(dat$NumHDonors) #what values
length(dat$NumHDonors)-length(na.omit(dat$NumHDonors) ) #how much NA

dat.neg <- dat[dat$NumHDonors>=1,]
dat.neg <- dat.neg[!is.na(dat.neg$NumHDonors),]

nrow(dat.neg)

write.csv(dat.neg, paste(workdir,"/",date,"_neg_TPs_suspectList.csv",sep=""),
          quote = TRUE,
          row.names = FALSE)

################
# creat pos suspect list
## double check
unique(dat$NumHAcceptors) #what values
length(dat$NumHAcceptors)-length(na.omit(dat$NumHAcceptors) ) #how much NA

dat.pos <- dat[dat$NumHAcceptors>=1,]
dat.pos <- dat.pos[!is.na(dat.pos$NumHAcceptors),]

nrow(dat.pos)

write.csv(dat.pos, paste(workdir,"/",date,"_pos_TPs_suspectList.csv",sep=""),
          quote = TRUE,
          row.names = FALSE)