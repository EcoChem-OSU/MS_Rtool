#############################################
## Title: XCMS Feature optimisation
###############################################################################
## version:1.2 - Modified from NT LuxWater patRoon paper, Jan 2023, DA
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
workPath <- "C:/Users/drozditb/Documents/OSU_data_analysis/NIST_20240501/"

## Input data -
sample.list <- "input/sample_list_nist.csv"
#sample.list <- "UKL_2023/UKL_2023_favs.csv"

# random select subset of the sample to save time
ratio <-1

## xcms parameter to optimize
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

## load data info
df <- read.table(sample.list, sep=",",header=TRUE)

anaInfo <- data.frame(cbind(path =df$path, # paste(df$path,"/",df$folder,sep=""),
                                analysis = df$filename,
                                group = df$group,
                                blank = df$blank) )
# select true sample only
anaInfo <- anaInfo[df$sampletype=="SA",]

anaInfo <- anaInfo[sample(nrow(anaInfo),size= 
                    round(nrow(anaInfo)*ratio,0) ), ]

print(anaInfo)
nrow(anaInfo)

# define parameter sets
pSet <-
  list(
      method = "centWave",
      ppm = ppm.range,
      min_peakwidth = min_peakwidth.range,
      max_peakwidth = max_peakwidth.range)
		
# otimize feature finding parameters
ftOpt <- optimizeFeatureFinding(anaInfo, "xcms3", pSet)

optimizedParameters(ftOpt) # get optimal values
