##############################################################
# MERGE SUSPECT LIST AND CHECK ..
#####################################
## version:4.1
## Date: 2026-05-19
## Author: Boris Droz @ Oregon State University
#############################################################
## DEPENDENCY AND FIRST INSTALL:
################################################################################
####  --- NEED TO INSTALL ANACONDA or OTHER PYTHON VERSION FIRST ---------------
################################################################################
#### NEED TO INSTALLL THE FOLLOWING STUFF:
# install.packages("remotes")
# remotes::install_github("LarsRichter82/larspack")
## install_miniconda(path = "C:/Users/drozditb/AppData/Local/miniconda", update = TRUE)

library(larspack) #RDKit for R https://zenodo.org/records/10275225
library(reticulate) # commande to activate and install python
library(stringr)
library(ChemmineR)
library(webchem)
library(rcdk)
library(rinchi)
library(RChemMass)
library(MetaboCoreUtils)
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
## PARAMETER
############
# folder contining all suspect list you want to merge
workdir <- "C:/Users/drozditb/Documents/GitHub/MS_Rtool/Suspect_list/"

## name for output list
out.list <- "20260519_testSL"

# remove non-PFAS (meanning suspect without F)
rem.noPFAS <- "YES"

## mass range of the mass spectrometer
mz_min= 100
mz_max= 1250

################################################################################
################################################################################
#                                   Start Script                                  
################################################################################
# NEW FUNCTION
is_empty_or_na <- function(x) { 
  if (missing(x) || is.null(x)) return(TRUE)  # Missing or NULL -> TRUE
  
  if (is.data.frame(x)) return(ncol(x) == 0)  # Data frame with 0 cols -> TRUE
  
  if (is.atomic(x) && length(x) == 1 && is.na(x)) return(TRUE)  # Single NA -> TRUE
  
  if (is.character(x) && all(trimws(x) == "")) return(TRUE)  # Empty string or only whitespace -> TRUE
  
  return(FALSE)  # Everything else -> FALSE
}

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
################################################################################
setwd(paste(workdir,"input",sep="") ) # set env

# Set outpath folder
date <- Sys.Date()
folder <- paste("/",date,"_Suspect_list", sep="")
outpath <- creat.subDir(paste(workdir,"/output",sep=""), folder)

fns <- list.files(paste(workdir,"/input",sep=""),pattern=".csv", full.names=TRUE)
name.list <- sub(".csv", "", list.files(paste(workdir,"input",sep=""),pattern=".csv", full.names=FALSE))

merge.list <- NULL
# check if empty data and merge list
for (i in 1:length(fns))
  {
    cat("LIST", name.list[i],"is under process...", "\n")  
  
     dat <- read.csv(fns[i], header=TRUE) #;print(names(dat))
  
      # compiled all data in std form
      compil_SL <- as.data.frame(matrix(NA, nrow = nrow(dat), ncol = 6))
      colnames(compil_SL) <- c("ID", "CAS","INCHIKEY","FORMULA","MONOISOTOPIC_MASS","SMILES")
      
     if(is_empty_or_na(dat[,names(dat)=="ID" |names(dat)=="DTXSID" |names(dat)== "SUSPECTID"|names(dat)=="Norman_SusDat_ID"|names(dat)=="DB"]))
        {compil_SL$ID <- paste(name.list[i],"_",seq(from=1,to=nrow(dat),by=1),sep="") }else
        {
          compil_SL$ID <- coalesce(
            dat$ID,
            dat$DTXSID,
            dat$SUSPECTID,
            dat$Norman_SusDat_ID,
            dat$DB
          )
          #compil_SL$ID <- dat[,names(dat)=="ID" |names(dat)=="DTXSID" |names(dat)== "SUSPECTID"|names(dat)=="Norman_SusDat_ID"|names(dat)=="DB"]}
        }
          if(is_empty_or_na(dat[,names(dat)=="CASRN" |names(dat)== "CAS"]) ) {}else
        {compil_SL$CAS <- dat[,names(dat)=="CASRN" |names(dat)== "CAS"] }
      if(is_empty_or_na(dat[,names(dat)=="INCHIKEY"| names(dat)=="StdInChIKey"| names(dat)=="InChIKey"])) {}else
        {compil_SL$INCHIKEY <- dat[,names(dat)=="INCHIKEY"| names(dat)=="StdInChIKey"| names(dat)=="InChIKey"] } 
      if(is_empty_or_na(dat[,names(dat)=="MOLECULAR.FORMULA" | names(dat)=="MOLECULAR_FORMULA"| names(dat)=="Formula"| names(dat)=="FORMULA"]) ) {}else
        {compil_SL$FORMULA <- dat[,names(dat)=="MOLECULAR.FORMULA" | names(dat)=="MOLECULAR_FORMULA"| names(dat)=="Formula"| names(dat)=="FORMULA"]}
      if(is_empty_or_na(dat[,names(dat)=="MONOISOTOPIC.MASS" | names(dat)=="MONOISOTOPIC_MASS"| 
                     names(dat)=="Monoisotopic_mass"| names(dat)== "FIXEDMASS" | names(dat)== "ExactMass"]) ) {}else
        {compil_SL$MONOISOTOPIC_MASS <- dat[,names(dat)=="MONOISOTOPIC.MASS" | names(dat)=="MONOISOTOPIC_MASS"| 
                                           names(dat)=="Monoisotopic_mass"| names(dat)== "FIXEDMASS"| names(dat)== "ExactMass" ]} 
     if(is_empty_or_na(dat[,names(dat)=="SMILES"| names(dat)=="SMILES_MS_ready"] )) {}else
     {compil_SL$SMILES <- dat[,names(dat)=="SMILES"| names(dat)=="SMILES_MS_ready"]}
  
      for (j in 1:nrow(compil_SL))
          {
            # check if empty line and search for info or compute to fill up...
            ## first fill up SMILES
            if (is_empty_or_na(compil_SL$SMILES[j])) { 
              if (is_empty_or_na(compil_SL$INCHIKEY[j])){
                if (is_empty_or_na(compil_SL$CAS[j])) {extracted_smiles <-"NA" }else{
                  cids <- get_cid(compil_SL$CAS[j],from = "cas", match = "first" )
                if(is.na(cids$cid)){extracted_smiles <-"NA"  }else{
                  extracted_smiles <- pc_prop(cids$cid, properties = "CanonicalSMILES")}
                }
              }else{
                    cids <- get_cid(compil_SL$INCHIKEY[j],from = "inchikey", match = "first" )
                    extracted_smiles <- pc_prop(cids$cid, properties = "CanonicalSMILES")
                  }
              if(length(extracted_smiles)>=2 ) {
                  if ( any(names(extracted_smiles)=="CanonicalSMILES") )
                      {compil_SL$SMILES[j] <- extracted_smiles$CanonicalSMILES}
                    else{compil_SL$SMILES[j] <-extracted_smiles[2] } 
                }
              
            }else{ # check if smiles is canonical
                
                compil_SL$SMILES[j] <- tryCatch({
                  
                  smiles <- as.character(compil_SL$SMILES[[j]])
                  
                  # Skip empty / NA
                  if (is.na(smiles) || smiles == "") return(NA)
                  
                  mol <- rdkit$MolFromSmiles(smiles)
                  
                  # Skip invalid molecules
                  if (is.null(mol)) return(NA)
                  
                  rdkit$MolToSmiles(mol)
                  
                }, error = function(e) {NA})

              
              if (!is.null(mol)) {
                compil_SL$SMILES[j] <- rdkit$MolToSmiles(mol)
              } else {
                compil_SL$SMILES[j] <- NA   # or keep original
              }
            } 
            ## second get inchikey if you had smile
            if(is_empty_or_na(compil_SL$SMILES[j])){ }else{
                if (is_empty_or_na(compil_SL$INCHIKEY[j])) { 
                  compil_SL$INCHIKEY[j] <- sapply(compil_SL$SMILES[j], get.inchi.key)
                  }
             }
            ## finally get CAS
                if (is_empty_or_na(compil_SL$CAS[j])) { 
                  extracted_cas <- tryCatch(
                      {
                          result <- cir_query(compil_SL$SMILES[j],
                                              "cas", match = "first")
                        if (is.null(result)) NA else result  # Assign NA if result is NULL
                      },
                    error = function(e) {
                        message("Error occurred: ", e$message)
                        return(NA)  # Assign NA if an error occurs
                    }
                  )
                     if(length(extracted_cas)==2 ) {
                        compil_SL$CAS[j] <- extracted_cas$cas}
                  }
            
            # compute formula and monoisotopic mass based on smiles if found before
            if (is_empty_or_na(compil_SL$SMILES[j])) { }else{
              
              if (!is.na(compil_SL$SMILES[[j]]) && nzchar(compil_SL$SMILES[[j]])) {
                compil_SL$FORMULA[j] <- tryCatch({
                  RChemMass::MolFormFromSmiles.rcdk(compil_SL$SMILES[[j]])
                }, error = function(e) {NA })
                
                compil_SL$MONOISOTOPIC_MASS[j] <- tryCatch({
                  MetaboCoreUtils::calculateMass(compil_SL$FORMULA[j])
                }, error = function(e) {NA })
                
              } else { compil_SL$FORMULA[j] <- NA 
                        compil_SL$MONOISOTOPIC_MASS[j] <- NA }
            }
      }  # end loop j
      
      compil_SL[] <- lapply(compil_SL, function(col) {
        if (is.list(col)) {
          sapply(col, function(x) paste(x, collapse = ";"))
        } else {
          col
        }
      })
      
      #save temp list 
      write.csv(compil_SL, paste(outpath,"/mod_",name.list[i],".csv",sep=""),
                row.names = FALSE)
      
      merge.list <- rbind(merge.list,compil_SL,
                  stringsAsFactors = FALSE) 
  } # end loop i

################################################################################
# save merge and edit manually then reload the data
################################################################################
merge.list -> PFAS.merged

PFAS.merged <- PFAS.merged %>% mutate(across(where(is.list), ~sapply(., toString)))

PFAS.merged <- PFAS.merged[!is.na(PFAS.merged$SMILES),]

# Match PFAS mass with HRMS analysis
PFAS.merged$MONOISOTOPIC_MASS <-as.numeric(PFAS.merged$MONOISOTOPIC_MASS)
PFAS.merged <- PFAS.merged[PFAS.merged$MONOISOTOPIC_MASS<=mz_max
                           & PFAS.merged$MONOISOTOPIC_MASS>=mz_min,]

# check out for duplicate based on
PFAS.merged <- PFAS.merged[ !duplicated(PFAS.merged$INCHIKEY),]

# remove PFAs associated with 13C or other label std
PFAS.merged <- PFAS.merged[!str_detect(PFAS.merged$FORMULA,fixed("[") ), ] 
PFAS.merged <- PFAS.merged[!str_detect(PFAS.merged$FORMULA,fixed(" ") ), ]
PFAS.merged <- PFAS.merged[!str_detect(PFAS.merged$FORMULA,fixed("D") ), ]

# remove nonPFAS (mean suspect without F)
if (rem.noPFAS=="YES") {PFAS.merged <- PFAS.merged[str_detect(PFAS.merged$FORMULA,fixed("F") ), ]}

#' #You need to set anaconda environment prior to running RDKit code
# reticulate::use_condaenv(condaenv = "my-rdkit-env", conda = "/Users/lrichter/miniconda3/bin/conda")

## check for negatif
PFAS.merged <- cbind(PFAS.merged, NumHDonors = rep(NA,nrow(PFAS.merged)) )

for (i in 1:length(PFAS.merged$SMILES)) { # you might need to run this a couple of time..
  # for (i in 4741:length(PFAS.list.nist$SMILES)) { 
  tryCatch( {PFAS.merged$NumHDonors[i] <-rdkit_descriptor(PFAS.merged$SMILES[i], 
                               descriptor ="NumHDonors")  
          } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

## check for positif
PFAS.merged <- cbind(PFAS.merged, NumHAcceptors = rep(NA,nrow(PFAS.merged)) )

for (i in 1:length(PFAS.merged$SMILES)) { # you might need to run this a couple of time..
  # for (i in 4741:length(PFAS.list.nist$SMILES)) { 
  tryCatch( {PFAS.merged$NumHAcceptors[i] <-rdkit_descriptor(PFAS.merged$SMILES[i], 
                                                          descriptor ="NumHAcceptors")  
  } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#################################################################################

# check with list
for (i in 1:length(fns))
{
  dat <- read.csv(fns[i], header=TRUE) #;print(names(dat))
  
  nam.SL <- name.list[i]
  
  #check if present in the list
  inchikey <- dat[ , names(dat) =="INCHIKEY" | names(dat) =="InChIKey" | names(dat) == "StdInChIKey"]
  
  te <- rep(0,nrow(PFAS.merged))
  
  te[PFAS.merged$INCHIKEY %in% inchikey]<- 1
  
  name.out <- c(names(PFAS.merged),nam.SL)
  
  PFAS.merged <- cbind(PFAS.merged,te)
  
  names(PFAS.merged) <- name.out                 
}

write.csv(PFAS.merged, paste(outpath,"/",date,"_SL_allcomp.csv",sep=""),
          quote = TRUE,
          row.names = FALSE, col.names = TRUE)

###############################################################################
# creat neg suspect list

## double check
unique(PFAS.merged$NumHDonors) #what values
length(PFAS.merged$NumHDonors)-length(na.omit(PFAS.merged$NumHDonors) ) #how much NA

PFAS.merged.neg <- PFAS.merged[PFAS.merged$NumHDonors>=1,]
PFAS.merged.neg <- PFAS.merged.neg[!is.na(PFAS.merged.neg$NumHDonors),]

nrow(PFAS.merged.neg)

write.csv(PFAS.merged.neg, paste(outpath,"/",date,"_neg_SL.csv",sep=""),
        quote = TRUE,
          row.names = FALSE, col.names = TRUE)

################
# creat pos suspect list

## double check
unique(PFAS.merged$NumHAcceptors) #what values
length(PFAS.merged$NumHAcceptors)-length(na.omit(PFAS.merged$NumHAcceptors) ) #how much NA

PFAS.merged.pos <- PFAS.merged[PFAS.merged$NumHAcceptors>=1,]
PFAS.merged.pos <- PFAS.merged.pos[!is.na(PFAS.merged.pos$NumHAcceptors),]

nrow(PFAS.merged.pos)

write.csv(PFAS.merged.pos, paste(outpath,"/",date,"_pos_SL.csv",sep=""),
          quote = TRUE,
          row.names = FALSE, col.names = TRUE)

