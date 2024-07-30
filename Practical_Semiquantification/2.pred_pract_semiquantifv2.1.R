####################################################
## Title: Practical Semi-quantification prediction
####################################################
## version:2.0
## Date: 2024-03-18
## Author: Dunping Cao & Chenyang Duan @ Oregon state University modified by B.Droz 
##
## Description:
# Predict semi-quantification using calibration model 

## Reference to cite
## Cao, D.; et al. J. Am. Soc. Mass Spectrom. 2023, 34 (5), 939–947. 
## https://doi.org/10.1021/jasms.3c00019.
###############################################################################
## PARAMETER -- MODIFIED IF NEEDED
############
## set your working directory --> refer to the model folder
workdir <- "R:/Boris Droz/script_HRMS/Practical Semiquantification/2024-07-30NEG_semiquant_model"

setwd(workdir) # do not change it

################################################################################
###############################################################################

## DO NOT MODIFIED BELOW JUST SELECT AND RUN THE CODE ####

################################################################################
###############################################################################
library("dplyr")
#####################Step 1: Import Dataset ##########################################
#Read data
df.cal <- read.csv(file.choose(), header= TRUE)
df.compounds <- read.csv("target_surrogatet_list_copy.csv", header= TRUE)

names(df.cal) <- sub("^X", "", names(df.cal)) # remove X if start by number

# homogenous sample names that the two df name fit
df.compounds$Name <- gsub(":", ".", df.compounds$Name)
df.compounds$Name <- gsub("-", ".", df.compounds$Name)
df.compounds$Name <- gsub(" ", ".", df.compounds$Name)

# x for concentration and response for A/A_average
data_set <- NULL
name_set <- NULL
# n.susp <- df.compounds$Name[!df.compounds$Type=="surrogate"]
n.sur <-  df.compounds$Name[df.compounds$Type=="surrogate"]

######################################################################################
## check if data set contain the same number of surrogate than used in the calibration
if (sum(names(df.cal)%in%n.sur) == length(n.sur)) {

    av.sur <- apply(df.cal[,names(df.cal)%in%n.sur],1,mean)
    
    susp <- select_if(df.cal[,!names(df.cal)%in%n.sur],is.numeric) # select suspect only
    
    susp <- susp/av.sur
    
    load(paste(workdir,"/cali_lm_weight",sep="")) # load the model
    
    beta <- coef(cal)["x"]
    
    # Predict concentration values
    predicted <- susp / beta
    
    text <- select_if(df.cal,is.character)
    
    predicted <- cbind(text,predicted)
    
   
    write.table(predicted, file=paste(workdir,"/suspect_pred_conc.csv",sep=""),
                sep=",", append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)

}else{
  
cat("WARNING: MODEL DO NOT HAVE SAME NUMBER OF SURROGATE WITH THE DATASET!!!!!!")  
}
