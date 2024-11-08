##############################################################
# KMD suspect list
#####################################
## version:1.1
## Date: 2024-02-29
## Author: Boris Droz @ Oregon State University
##############################################################
# DESCRIPTION:
# create a KMD suspect list for a specific base unit

## Dependency:
## residual function modified version of...
## https://github.com/usnistgov/NISTPFAS/tree/main/suspectlist/fn

#############################################################
## PARAMETER
# library
library(MetaboCoreUtils)
library(xlsx)
library(stringr)
workdir <- dirname(rstudioapi::getSourceEditorContext()$path) # get directory work only in rstudio
# workdir <- "R:/Boris Droz/script_HRMS/calc_KMD/" # directory where you put the code
source(paste(workdir,"/func_calc_residual.R",sep=""))
## choose a file
fns <- file.choose()

form_unit <- "CF2" #base unit
##############################################################

workdir <- dirname(fns) # file directory

df.SL <- read.csv(fns)

# caculate MM from formula if not already
# df.SL$MONOISOTOPIC_MASS <- calculateMass(df.SL$FORMULA)

# calculate mass deffect
df.SL$MD <- df.SL$MONOISOTOPIC_MASS - round(df.SL$MONOISOTOPIC_MASS, 0)

min(df.SL$MD)
max(df.SL$MD)

# Calculate KMD
dmz <- calculateMass(form_unit)
KM <- df.SL$MONOISOTOPIC_MASS *round(dmz,0)/dmz
df.SL$KMD <- round(KM,0)- KM

# modulo to classified the HS
df.SL$modulo <- df.SL$MONOISOTOPIC_MASS %% dmz

# calculate residual using nist formula and append to the suspect list
resi <- NULL
for (i in 1:nrow(df.SL))
    {
      resi <- rbind(resi,calculate_residual(df.SL$FORMULA[i],
                  rep_unit = form_unit, df.SL$MONOISOTOPIC_MASS[i]) )
    }

df.SL <-cbind(df.SL, resi)

# remove data if not rep unit present
crit <- "You do not have all of the elements in the repeating units in the elemental formula"
df.SL <-df.SL[!df.SL$rep_unit==crit,]

write.csv(df.SL, paste(workdir,"/KMD_",form_unit,
                                 "_",basename(fns),sep=""),
          quote = TRUE,
          row.names = FALSE, col.names = TRUE)


