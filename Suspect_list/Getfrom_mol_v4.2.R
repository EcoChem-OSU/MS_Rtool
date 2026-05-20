################################################################################
###    			
###  Get molecular info from *.mol file to creat a suspect list  ###   
###
################################################################################
################################################################################
## v3.0 April 2020 - B.Droz unistra
## v3.1 
# v4.1 more generic version
##########
## Historic
#---------
#  tested with R 3.5.3
# v2.0 Debugged by Nan Xiao April 2020
###############################################
####################
################
## Descriptive
##############
# Batch transform mol file fromaninput folder (or any molecule identification) 
## into smile
# calculate INCHIKEY, FORMULA,	MONOISOTOPIC_MASS and	SMILES

################################################################################
## Dong-Sheng Cao, Nan Xiao, Qing-Song Xu, and Alex F. Chen. (2015). Rcpi: R/Bioconductor 
## package to generate various descriptors of proteins, compounds and their interactions. Bioinformatics 31 (2), 279-281.
## https://github.com/nanxstats/Rcpi
##
#################################################################
## ---- NECESSARY TO GET THIS FOLLOWING LIBRARY FROM THIS SOURCE
#######################################################################
# install.packages("BiocManager")
# BiocManager::install("Rcpi")
# BiocManager::install("Rcpi", dependencies = c("Imports", "Enhances"))
# BiocManager::install("ChemmineOB")
#
# find local dowload and unpack it ## if needed
# utils:::menuInstallLocal()
##
##
##############################################################################################
#####################
#Loading libraries
#####################
packages <- c("callr","ChemmineR","Rcpi", "RCurl","yaml","rlang","tibble", 
              "stringi", "ChemmineOB", "rinchi","rcdk", "RChemMass","MetaboCoreUtils")

# load library
for(pkg in packages){library(pkg,character.only=TRUE) }
 
####################################################################################################################################
## SCRIPT PARAMETER
#####################
rm(list=ls()) # remove all object in workspace to start fresh

# Set folder 
workdir <- "C:/Users/drozditb/Documents/OSU_data_analysis/biomimetic_mol/"

out.name <- "EPA_mol_bimimetic.csv"

setwd(workdir) # set dir

###############################################################################
###############################################################################

## DO NOT MODIFIED BELOW JUST RUN THE CODE ####

##   ----  SCRIPT START HERE ----
###############################################################################

  ## get file path and id
  fns <- list.files(paste(workdir,"/input",sep="") ,pattern=".mol$",full.names = TRUE) 
  fid <- gsub(".mol", "", list.files(paste(workdir,"/input",sep=""),pattern=".mol$",full.names = FALSE))# ; print(fid)
  
  ### FROM Nan Xiao chunks mode
  #############################
  convert <- function (fns, idx) {
    callr::r(function (fns, idx) {
      smiles <- c()
      for (k in idx) {
        Rcpi::convMolFormat(infile = fns[k], outfile = "temp.smi", from = "mol", to = "smiles")
        smiles <- c(smiles, Rcpi::readMolFromSmi(smifile = "temp.smi", type = "text")[1])
      }
      smiles
    }, args = list(fns, idx))
  }
  
  p <- length(fns)
  chunks <- split(1:p, ceiling(seq_along(1:p)/400))
  smi <- rep(NA, p)
  for (k in 1:length(chunks)) smi[chunks[[k]]] <- convert(fns, chunks[[k]])
  ###################################################################################
  inchikey <- sapply(smi, get.inchi.key) # get inchikey from smiles
  
  formula <- NULL
  for (i in 1:length(smi)) { formula <-  c(formula, MolFormFromSmiles.rcdk(smi[i])) }
  
  mono_mass <- calculateMass(formula)
  
  d.out <- data.frame(ID=fid , SMILES=smi, INCHIKEY=inchikey, FORMULA=formula,
                        MONOISOTOPIC_MASS=mono_mass) 

# write the output table
	write.csv(d.out,file=paste(workdir,"/",out.name ,sep=""),
            append=FALSE,row.names=FALSE,col.names=TRUE,quote=TRUE)
	#############################################################################
	