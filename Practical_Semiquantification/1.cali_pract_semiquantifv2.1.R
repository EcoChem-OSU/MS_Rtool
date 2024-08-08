####################################################
## Title: Practical Semi-quantification calibration
####################################################
## version:2.0
## Date: 2024-03-18
## Author: Dunping Cao & Chenyang Duan @ Oregon state University modified by B.Droz 
##
## Description:
################
## Compute an semiquantitative average calibration using the area of target
## divided by the average area of surrogate (deutered intern stad) areas in 
## function of target concentration in units of nmoles/L.
## Use a weighted linear regression models.

################################################################################
## Reference to cite
## Cao, D.; et al. J. Am. Soc. Mass Spectrom. 2023, 34 (5), 939–947. 
## https://doi.org/10.1021/jasms.3c00019.
################################################################################
################################################################################
##################################
## PARAMETER -- MODIFIED IF NEEDED
############
## set your working directory
workdir <- "C:/Users/drozditb/Documents/GitHub/MS_Rtool/Practical_Semiquantification"

setwd(workdir) # do not change it

ccal <- "060121_ccal.csv" # file with data of your calibration curve - 
                          ## concentration in function of Area count 
                          ## concentration should be the first column
unit.label <- "nmol/L" # unit of the concentration label for graph do not 
                      # influence the calculation but should be specified

compounds <-"target_surrogatet_list.csv" # file with compounds information -
                      ## Name of compounds --> should be similar in all files ccal and compounds
                      ## Chemical formula --> be aware of the format  
                      ##                  type "?calculateMass" 
                      ## Type: Target (std) or surrogate (internal std)
  
In.mass <- "YES" # option for mass in gram or substance in mole based unit 
                # Yes if input unit from the calicurve are in ng/L, g/L,.... default
                # No if input unit are in mmol/L, mol/L,.... 

alpha <- 0.05 # confidence interval for the incertitude calculation 
              ## default 0.05

output.folder.name <- "NEG_semiquant_model" #name the output folder

################################################################################
###############################################################################

## DO NOT MODIFIED BELOW JUST SELECT AND RUN THE CODE ####

################################################################################
###############################################################################

source("FUN_Pract_semiquantifv2.0.R") # hand made function from Cao 2023 and Droz 2024

# creat output
date <- Sys.Date()
output <- creat.subDir(workdir,paste(date,output.folder.name,sep="") )
# check package 
check.lib(c("tidyverse","dplyr") )

# check lib metabocore
libTest <- try(library("MetaboCoreUtils",character.only=TRUE),silent=TRUE)

if(class(libTest)=='try-error'){

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MetaboCoreUtils")}

#####################Step 1: Import Dataset ##########################################
#Read data
df.cal <- read.csv(ccal, header= TRUE)
df.compounds <- read.csv(compounds, header= TRUE)

# save copy of both files
write.table(df.cal, file=paste(output,"/cali_curve_copy.csv",sep=""),
            sep=",", append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)

write.table(df.compounds, file=paste(output,"/target_surrogatet_list_copy.csv",sep=""),
            sep=",", append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)

names(df.cal) <- sub("^X", "", names(df.cal)) # remove X if start by number

# homogenous sample names that the two df name fit
df.compounds$Name <- gsub(":", ".", df.compounds$Name)
df.compounds$Name <- gsub("-", ".", df.compounds$Name)
df.compounds$Name <- gsub(" ", ".", df.compounds$Name)

# x for concentration and response for A/A_average
data_set <- NULL
name_set <- NULL
n.targ <- df.compounds$Name[df.compounds$Type=="target"]
n.sur <-  df.compounds$Name[df.compounds$Type=="surrogate"]

############################################

av.sur <- apply(df.cal[,names(df.cal)%in%n.sur],1,mean)

for (i in 1:ncol(df.cal[,names(df.cal)%in%n.targ]))
{
  response <- df.cal[,names(df.cal)==n.targ[i]]/av.sur
    if (In.mass=="YES"){
      formula <- df.compounds$Chemical.Formula[df.compounds$Name == n.targ[i] ]
      MM <- calculateMass(formula)
      x <- df.cal[,1]/MM
    } else{ x <- df.cal[,1] }
  
  data_set <-rbind(data_set,data.frame(x=x,response=response))
  name_set <- c(name_set, rep(names(df.cal)[names(df.cal)==n.targ[i]],length(x)) )
}

##Step 2: Perform weighted linear regression and calculate uncertainty###################
#########################################################################################
####2.1 Define weight, return k and b
ywt <- 1/(data_set$x)
cal <- lm(formula = response ~ 0 + x, data=data_set, weights = ywt)

# Save Model - R file
save(list=c('cal'),
     file=paste(output,"/cali_lm_weight",sep=''))

#save cali dataset
write.table(data_set, file=paste(output,"/cali_dataset.csv",sep=""),
            sep=",", append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)

png(filename = paste(output,"/Semiquant_cali_plot.png",sep=""),width = 480, height = 480 )

par(mar=c(5, 5, 2, 2) )
####plot this weighted linear regression
plot(data_set$x,data_set$response,
     xlab=paste("Target concentration (", unit.label,")",sep=""),
       ylab= "Response (A_target/Avg A_surrogate)",
     las=1,
     pch=19)
       
abline(cal, lwd=2)

dev.off()
####plot the residual plot of this weighted linear regression
png(filename = paste(output,"/Semiquant_residual_plot.png",sep=""),width = 480, height = 480 )

par(mar=c(5, 5, 2, 2) )
####plot this weighted linear regression
plot(data_set$x,cal$residuals,
     xlab=paste("Target concentration (", unit.label,")",sep=""),
     ylab= "Residual",
     las=1,
     pch=19)

dev.off()

## model values and performances
slope <-  cal$coefficients
r.squared <-summary(cal)$r.squared
RMSE <- sqrt(sum(cal$residuals^2) / cal$df)

f.info <- paste(output,"/Model_Perfomance.txt",sep="")

cat("######################################################",file= f.info,append=TRUE, sep="\n")
cat( paste("Practical Semi-quantification calibration", Sys.Date()), file= f.info, sep="\n")
cat("######################################################",file= f.info,append=TRUE, sep="\n")
cat("R-script 1.cali_pract_semiquantifv2.0",file= f.info,append=TRUE, sep="\n")
cat("######################################################",file= f.info,append=TRUE, sep="\n")
cat(paste("Model name:", output.folder.name), file= f.info,append=TRUE, sep="\n")
cat(paste("Unit:", unit.label), file= f.info,append=TRUE, sep="\n")
cat(paste("nb target:",ncol(df.cal[,names(df.cal)%in%n.targ]) ),
    file= f.info,append=TRUE, sep="\n")
cat(paste("nb surrogate:",ncol(df.cal[,names(df.cal)%in%n.sur]) ),
    file= f.info,append=TRUE, sep="\n") 
if (In.mass=="NO"){ }else{ 
    cat("INPUT WERE IN KG AND CONVERTED IN MOL USING FORMULA FOR EACH MOLECULES",
           file= f.info,append=TRUE, sep="\n") }
  
cat(paste("confidence interval:", alpha),file= f.info,append=TRUE, sep="\n")
cat(paste("Slope:", round(slope,3) ),file= f.info,append=TRUE, sep="\n")
cat(paste("R2:",round(r.squared,4) ),file= f.info,append=TRUE, sep="\n")
cat(paste("RMSE:",round(RMSE,3) ),file= f.info,append=TRUE, sep="\n")

####3.2 Define alpha and x_series (x), return uncertainty (predicted intervall)
###########################################################################
y_series <- cal$fitted.values

pred.int <- table_uncertainty(y_series,alpha,data_set,index=1) # index=2 for 1/x2,  
                                                              # index=1 for 1/x
                                                              # index=0 for no weighting

#save cali dataset
write.table(pred.int, file=paste(output,"/predicted_interval.csv",sep=""),
            sep=",", append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)


