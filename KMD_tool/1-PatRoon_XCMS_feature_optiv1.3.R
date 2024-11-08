#############################################
## Title: XCMS Feature optimisation
###############################################################################
## version:1.3 - Modified from NT LuxWater patRoon paper, Jan 2023, DA
## 			 https://doi.org/10.1186/s12302-023-00805-5
## Date: 06 December 2023

## Tutorial and Handbook on https://github.com/rickhelmus/patRoon
## Depends:
##        R(>=4.3.1)
##        patRoon(>=2.3.0)
##
## Description:
# feature parameter optimization steps (based on IPO)
# Optimization data and instrument specific
###############################################################################
###############################
## PARAMETER -- MODIFIED IF NEEDED
#############
## set your working directory
workPath <- "C:/Users/Peter/Documents/OSU/Fields Rotation/GitHub Projects/Example Input/KMD_7600_20240221_Spike_cal/"

## Input data -
sample.list <- "sample_list_KMD.csv"
#sample.list <- "UKL_2023/UKL_2023_favs.csv"

# random select subset of the sample to save time
ratio <-1

## Range of xcms parameter to optimize
ppm.range = c(8,17)
min_peakwidth.range = c(6, 15)
max_peakwidth.range = c(80, 140)
###############################################################################
###############################################################################

## DO NOT MODIFIED BELOW JUST RUN THE CODE ####

##   ----  SCRIPT START HERE ----
###############################################################################
# -------------------------
# initialization
# -------------------------
library(patRoon)

setwd(workPath)
datum <- Sys.Date()

## load data info
df <- read.table(paste("input/",sample.list,sep=""), sep=",",header=TRUE)

anaInfo <- data.frame(cbind(path =df$path, # paste(df$path,"/",df$folder,sep=""),
                                analysis = df$filename,
                                group = df$group,
                                blank = df$blank) )
# select true sample only
anaInfo <- anaInfo[df$sampletype=="SA",]

anaInfo <- anaInfo[sample(nrow(anaInfo),size= 
                    round(nrow(anaInfo)*ratio,0) ), ]

# define parameter sets
pSet <-
  list(
      method = "centWave",
      ppm = ppm.range,
      min_peakwidth = min_peakwidth.range,
      max_peakwidth = max_peakwidth.range)
		
# otimize feature finding parameters
ftOpt <- optimizeFeatureFinding(anaInfo, "xcms3", pSet)

optimizedParameters(ftOpt) 
capt.out <- capture.output(show(ftOpt)) # capture optimized parameter

## save optimized parameter
outpath <- paste(workPath,"/output",sep="")
f.info <- paste(outpath,"/",datum,"_XCMS_optPara.txt",sep="")
cat( paste("*** XCMS optimization parameter....", datum), file= f.info, append=TRUE, sep="\n")
cat( "#########################################################", file= f.info, append=TRUE,sep="\n")
cat( paste("SampleList: ", workPath,"/input/",sample.list,sep=""), file= f.info, append=TRUE,sep="\n")
cat( paste("Ratio of samples: ", ratio, sep=""), file= f.info, append=TRUE,sep="\n")
cat( paste("Number of samples: ", nrow(anaInfo), sep="" ), file= f.info, append=TRUE,sep="\n")
cat( capt.out, file= f.info, append=TRUE,sep="\n")
