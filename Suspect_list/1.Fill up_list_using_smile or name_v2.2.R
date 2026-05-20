##############################################################
# Modified version of MERGE SUSPECT LIST AND CHECK
#####################################
## version:2.2
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
library(reticulate)
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
# workdir <- "C:/Users/drozditb/Documents/GitHub/MS_Rtool/Suspect_list"

# setwd(workdir)

fns <- file.choose()

################################################################################
################################################################################
#                                   Start Script                                  
################################################################################
# NEW FUNBCTION B.Droz chatGbt 2025-030-06
# is_empty_or_na <- function(x) { is.na(x) || x == "" || trimws(x) == "" }
is_empty_or_na <- function(x) { 
  if (missing(x) || is.null(x)) return(TRUE)  # Missing or NULL -> TRUE
  
  if (is.data.frame(x)) return(ncol(x) == 0)  # Data frame with 0 cols -> TRUE
  
  if (#is.atomic(x) && length(x) == 1 && 
    is.na(x)) return(TRUE)  # Single NA -> TRUE
  
  if (## is.character(x) && 
      all(trimws(x) == "")) return(TRUE)  # Empty string or only whitespace -> TRUE
  
  return(FALSE)  # Everything else -> FALSE
}

################################################################################
  
     
     dat <- read.csv(fns, header=TRUE) ;print(names(dat))
  
      # compiled all data in std form
      compil_SL <- dat
     
      ### FIRST make the table and creat column if missing....
      missing_cols <- setdiff(
                  c("SMILES", "INCHIKEY", "CAS", "FORMULA", "MONOISOTOPIC_MASS"),
                  names(compil_SL))
      
      compil_SL[missing_cols] <- NA
    
      for (j in 1:nrow(compil_SL))
          {
            # check if empty line and search for info or compute to fill up...
            ## first fill up SMILES
            if (is_empty_or_na(compil_SL$SMILES[j])) {
              if (is_empty_or_na(compil_SL$INCHIKEY[j])){
                if (is_empty_or_na(compil_SL$CAS[j])) {
                  if (is_empty_or_na(compil_SL$Name[j])){extracted_smiles <-"NA"
                    
                  }else{
                    cids <- get_cid(compil_SL$Name[j],from = "name", match = "first" )
                    if(is.na(cids$cid)){extracted_smiles <-"NA"  }else{
                      extracted_smiles <- pc_prop(cids$cid, properties = "CanonicalSMILES")}
                  }
                  
                   }else{
                  cids <- get_cid(compil_SL$CAS[j],from = "cas", match = "first" )
                if(is.na(cids$cid)){extracted_smiles <-"NA"  }else{
                  extracted_smiles <- pc_prop(cids$cid, properties = "CanonicalSMILES")}
                }
              }else{
                    cids <- get_cid(compil_SL$INCHIKEY[j],from = "inchikey", match = "first" )
                    extracted_smiles <- pc_prop(cids$cid, properties = "CanonicalSMILES")
                  }
              if(length(extracted_smiles)==2 ) {
                      compil_SL$SMILES[j] <- extracted_smiles[2]}
              
            }else{ 
              compil_SL$SMILES[j] <-gsub(" ", "", compil_SL$SMILES[j]) #check if any space on the smiles code if so remove it

            } 
            ## second get inchikey if you had smile
            if(is_empty_or_na(compil_SL$SMILES[j])){ }else{
                if (is_empty_or_na(compil_SL$INCHIKEY[j])) { 
                  compil_SL$INCHIKEY[j] <- as.character(sapply(compil_SL$SMILES[j], get.inchi.key))
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
                compil_SL$FORMULA[j] <- RChemMass::MolFormFromSmiles.rcdk(as.character(compil_SL$SMILES[j]))
            
                compil_SL$MONOISOTOPIC_MASS[j] <- MetaboCoreUtils::calculateMass(compil_SL$FORMULA[j])
            }
      }  # end loop j
      
     
 # check if table is proper     
      
      as_scalar_chr <- function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        if (is.list(x)) x <- x[[1]]
        if (length(x) > 1) x <- paste(as.character(x), collapse = "; ")
        as.character(x)
      }
      
      for (nm in names(compil_SL)) {
        if (is.list(compil_SL[[nm]])) {
          compil_SL[[nm]] <- vapply(compil_SL[[nm]], as_scalar_chr, character(1))
        }
      }
      
      write.csv(compil_SL, fns, quote = TRUE, row.names = FALSE)
      

