################################################################################
#### MALBAC_R
################################################################################
################################################################################
# version:3.3
## Date: September 2024
## Author: Boris Droz @ Oregon State University
################################################################################
## Description:
###############
## Provide metric to evaluate if normalization of the feature is needed
## Test several normalization method for HRMS data
## 
###############################################################################
# Reference:
##  Leach, D. T.; et al. Anal. Chem. 2023, 95 (33), 12195-12199. 
##             https://doi.org/10.1021/acs.analchem.3c01289.
## Ding, X.; et al. Anal. Chem. 2022, 94 (21), 7500-7509. 
##            https://doi.org/10.1021/acs.analchem.1c05502
##
###############################################################################
## Parameter
############
## set your working directory
## workPath <- "R:/Boris Droz/Proj_UKL_2023/data_analysis/"
workPath <- "C:/Users/drozditb/Documents/OSU_data_analysis/"

setwd(workPath) # do not change it

inpath <- paste(workPath ,"/output/",sep="")

fn.df <- "/raw_aligned_grouped.txt" # feature list grouped
fn.sample.list <- "input/sample_list_UKL2023_outliersRemoved.csv" # sample list with all info see example

thr_range <- c(100,500) # range to replace NA and zero should be same than 
                        # threshold used in peak picking filtering
################################################################################
###############################################################################

## DO NOT MODIFIED BELOW JUST SELECT AND RUN THE CODE ####

################################################################################
###############################################################################

###############################################################################
# library
library(malbacR)
library(pmartR)
library(ggplot2)
library(stringr) 
library(factoextra)
library(vegan) # NMDS
library(RColorBrewer)

###############################################################################
# #############################################################################
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
###############################################################################
## START SCRIPT
###############
# Set outpath folder
date <- Sys.Date()
folder <- paste("/",date,"_MALBAC_R", sep="")
outpath <- creat.subDir(inpath, folder)

f.info <- paste(outpath,"/AA_INFO_RUN_README.txt",sep="")
cat( paste("*** malbacR result ****", Sys.Date()), file= f.info, append=TRUE, sep="\n")
cat( "##################################", file= f.info, append=TRUE,sep="\n")
cat( paste("Featuregroup:", fn.df), file= f.info, append=TRUE,sep="\n")
cat( paste("SampleList:", fn.sample.list), file= f.info, append=TRUE,sep="\n")
cat( paste("Zero convertion to random between:", thr_range[1],"-",thr_range[2]), file= f.info, append=TRUE,sep="\n")
cat( "      ", file= f.info, append=TRUE, sep="\n")

## load data
df <- read.table(paste(inpath,fn.df,sep=""), sep="\t", header=TRUE)
sample.list <- read.csv(fn.sample.list, sep=",")

# get SA and QC only
sample.list <- sample.list[sample.list$sampletype=="SA"|sample.list$sampletype=="QC" ,]

names(df) <- sub("^X", "", names(df)) # rename if name start by X -- 
sample.list$filename <- str_replace_all(sample.list$filename,"-",".") #unified name

#select QC                                                      
df.QC <- df[,names(df) %in% sample.list$filename[sample.list$sampletype=="QC"] ]
# Select SA
df.SA <- df[,names(df) %in% sample.list$filename[sample.list$sampletype=="SA"] ]

df.all <- cbind( df.SA,df.QC)

################################################################
#### check if need QC correction
## peak RSD% should be less than 20% if yes no correction needed.
## DOI: 10.1021/acs.analchem.1c05502.

##For QC
sd.QC <- apply(df.QC, 1, sd)
mean.QC <- apply(df.QC, 1, mean)

RSD.QC <- sd.QC/mean.QC*100

png(filename =paste(outpath, "/QC_boplot.png", sep=""),
    width = 480, height = 480, units = "px")
  b.QC <-boxplot(RSD.QC)
dev.off()

### check similar for samples replicate if ...
## the RSD in sample smaller than QC use IS normalization

samp.group <-sample.list$group[!sample.list$group=="QC"]
RSD.SA <- NULL

for (i in 1:length(unique(samp.group)) )
  {
  s.samp <- sample.list$filename[sample.list$group== samp.group[i]]
  df.SA.in <- df.SA[,names(df.SA)==s.samp]
  sd.SA <- apply(df.SA.in, 1, sd)
  mean.SA <- apply(df.SA.in, 1, mean)
  RSD.SA <- cbind(RSD.SA,sd.SA/mean.SA*100)
}

colnames(RSD.SA) <- unique(samp.group)

png(filename =paste(outpath, "/SA_boplot.png", sep=""),
    width = 960, height = 480, units = "px")
  b.SA <-boxplot(RSD.SA,
                 las=2)
dev.off()

## check raw data -- to evaluate the need of normalization
res.pca <- prcomp(t(df.SA), 
                  scale = FALSE)
# predict PCA coordinate
QC.coord <- predict(res.pca, newdata = t(df.QC))             

### a) plot PCA
###############
# plot PCA with QC

png(filename =paste(outpath, "/PCA_raw.png", sep=""),
    width = 480, height = 480, units = "px")

  p <- fviz_pca_ind(res.pca, 
                  geom="point", 
                  pointsize = 4,
                  title = "raw_sample")
  
  fviz_add(p, QC.coord, geom="point", 
          color ="red",
          addlabel = FALSE )
  
dev.off()

#reorganise the data
####################
e_data <- data.frame(cbind(Molecule=df$group,
                           df.all))

names(e_data) <- sub("^X", "", names(e_data)) # rename if name start by X -- 

# remove feature with all zero value
e_data <- e_data[!apply(e_data[,2:ncol(e_data)],1,sum)==0,]

# replace existing zero by 100 to threshold data (random)
seq.num <-seq(from=thr_range[1],to=thr_range[2],by=1)
samp.numb <- sample(seq.num, length(e_data[e_data==0]), replace = TRUE)

e_data[e_data==0] <- samp.numb
#########################################################
f_data <- data.frame(cbind(SampleID=sample.list$filename,
                             Injection_order=sample.list$time, 
                             group= sample.list$sampletype,
                              batch= sample.list$batch )) 

f_data$group[f_data$group=="SA"] <- "1"
f_data$group[f_data$group=="QC"] <- "QC" 

### just reorganized and ....
f_data <-f_data[order(as.numeric(f_data$Injection_order)),]
# rename order and make sur to have an integrer
f_data$Injection_order <-as.integer(seq(from=1, to= nrow(f_data), by=1))

# creat fake batches if necessary
if (length (unique(f_data$batch) ) >1){
        myLetters <- LETTERS[1:26]
  
        f_data$batch <- match(f_data$batch, myLetters)
  
  }else{
    cat( "*** Only one batch --> 3 equal batches are created  ---", file= f.info, append=TRUE, sep="\n") 
   
        incr<- ceiling(nrow(f_data)/3)
        f_data$batch[1:incr] <- 1
        f_data$batch[(incr+1):(2*incr)] <- 2
        f_data$batch[(2*incr+1):nrow(f_data)] <- 3
        f_data$batch <- as.integer(f_data$batch)
}

#######################################
ncol(e_data)-1 # check should be similar
nrow(f_data)###########################
########################################
#########################################################################
# creat an object for analysis
#############################
mymetabo <- as.metabData(
          e_data= e_data,
          f_data= f_data,
          e_meta = NULL,
          edata_cname="Molecule",
          fdata_cname="SampleID",
          emeta_cname = NULL,
          techrep_cname = NULL)

## creat similar object with duplicate QC for TIGER and SERRF...
## duplicate QC values
pos <-f_data$group=="QC"
dupli_rows <- f_data[pos, ]

dupli_rows$SampleID <- paste(dupli_rows$SampleID,"_SA23",sep="")
dupli_rows$group <- "1"

mod_f_data <- rbind(f_data,dupli_rows) # merge row\
mod_f_data <- mod_f_data[order(mod_f_data$Injection_order),] # reorder by Injection

mod_f_data$Injection_order <- seq(from=1,to=nrow(mod_f_data), by=1)# reigni Injection order

# creat mode e data
add_e_data <- e_data[names(e_data) %in% f_data$SampleID[f_data$group=="QC"] ]
names(add_e_data) <- paste(names(add_e_data),"_SA23",sep="")
mod_e_data <- cbind(e_data,add_e_data )

mod_mymetabo <- as.metabData(
                e_data= mod_e_data,
                f_data= mod_f_data,
                e_meta = NULL,
                edata_cname="Molecule",
                fdata_cname="SampleID",
                emeta_cname = NULL,
                techrep_cname = NULL)
        
#need transformation
mymetabo <- group_designation(mymetabo,main_effects = "group",batch_id = "batch")
mod_mymetabo <- group_designation(mod_mymetabo,main_effects = "group",batch_id = "batch")

mymetabo_2log <- edata_transform(mymetabo,
                                 data_scale = "log2")

mymetabo_norm <- normalize_global(mymetabo_2log,subset_fn = "all",
                                  norm_fn = "median",
                                     apply_norm = TRUE,
                                  backtransform = TRUE)

############################# #########  ###########  #########################
# SCALING METHODS
####################
# range scaling
range <- bc_range(omicsData = mymetabo_2log)

# power scaling
power <- bc_power(omicsData = mymetabo_2log)

# pareto scaling
pareto <- bc_pareto(omicsData = mymetabo_2log)

# QUALITY CONTROL METHODS
#########################
# TIGER 
tigerFilt <- tiger_filter(mod_mymetabo,  #pmart_amide,
                          sampletype_cname = "group",test_val = "QC")
pmart_Filt <- apply_tigerFilt(tigerFilt,mod_mymetabo)
tiger_abundance <- bc_tiger(omicsData = pmart_Filt,
                                  sampletype_cname = "group",
                                  test_val = "QC",
                                  injection_cname = "Injection_order",
                                  group_cname = "group")
tiger <- edata_transform(omicsData = tiger_abundance, 
                               data_scale = "log2")

# QC-RLSC
qcrlsc <- tryCatch( { bc_qcrlsc(omicsData = mymetabo_2log,
                          block_cname = "batch",
                          qc_cname = "group", 
                          qc_val = "QC", 
                          order_cname = "Injection_order",
                          missing_thresh = 0.5, 
                          rsd_thresh = 0.3, 
                          backtransform  = FALSE,
                          keep_qc = TRUE)}, 
                  error = function(e) {
                "The first and last sample aquisition run for each batch must be a QC sample" })
if (is.character(qcrlsc) ) {
  cat( "qcrlsc not run because the first and last sample aquisition run for each batch must be a QC sample",
       file= f.info, append=TRUE, sep="\n")
}
  
# SERRF 
serrf_abundance <- bc_serrf(omicsData = mod_mymetabo,
                            sampletype_cname = "group",
                            test_val = "QC",
                            group_cname = "group")
serrf <- edata_transform(omicsData = serrf_abundance, 
                         data_scale = "log2")
# QC-RFSC
qcrfsc_abundance <- bc_qcrfsc(omicsData = mymetabo,
                              qc_cname = "group",
                              qc_val = "QC",
                              order_cname = "Injection_order",
                              group_cname = "group",
                              ntree = 500, keep_qc = TRUE)
qcrfsc <- edata_transform(omicsData = qcrfsc_abundance,
                          data_scale = "log2")

# OTHER METHODS
###############
# ComBat
combat <- bc_combat(omicsData = mymetabo_norm, 
                    use_groups = FALSE)

# EigenMS
eigen <- bc_eigenMS(omicsData = mymetabo_2log)

# WaveICA2.0
wave_abundance <- bc_waveica(omicsData = mymetabo, 
                             batch_cname = "batch",
                                   injection_cname = "Injection_order",
                                   version = "WaveICA2.0", alpha = 0, 
                                   cutoff_injection = 0.1, K = 10)
                                   #negative_to_na = TRUE)
wave <- edata_transform(omicsData = wave_abundance, data_scale = "log2")

####################################################################################################
####################################################################################################
## MERGE ALL NORM MODEL>>>>> Change the next to object if necessary
###############################################################################
all.results <-list(mymetabo,range,power,pareto, tiger,qcrlsc,
                   serrf, qcrfsc,
                    combat, eigen,wave)

# creat similar name list
tile.list <- c("unadjusted","range","power","pareto","tiger","qcrlsc",
                      "serrf", "qcrfsc","combat", "eigen","wave")

# Function to check if a value is NULL or text-only
is_valid <- function(x) {
  !is.null(x) && !is.character(x)}

indices <- which(sapply(all.results, is_valid))

# Filter the list to omit NULL or text-only variables
all.results <- Filter(function(x) is_valid(x), all.results)
tile.list <- tile.list[indices]

# INSPECT RESULTS 
#################

all.RSD <-NULL
all.dist <- NULL
all.Rsquare <-NULL

for (i in 1: length(all.results))
  {
     cat("INSPECT RESULT", tile.list[i],"...", "\n")
  
    all.data <-all.results[[i]]$e_data
    
    if (tile.list[i]=="tiger"|tile.list[i]=="serrf")
        { # removw SA23 from the QC name
        names(all.data) <- sub("_SA23", "", names(all.data))
      }else{}
    
  names(all.data)[1] <-"group"
      
	write.table(all.data, file=paste(outpath, "/featureGroups_normby_",tile.list[i],".csv", sep=""),
            append = FALSE, quote = TRUE, sep = ",",
            row.names = FALSE,col.names = TRUE )
	
	# create average data by group without blank and QC
	sa.group <- unique(sample.list$group[sample.list$sampletype=="SA"])
	 
	av.data <- all.data[,1]
	for (j in 1:length(sa.group))
	    {
	      pos <- match(sample.list$filename[sample.list$group==sa.group[j]], names(all.data))
	      av.data <-cbind(av.data, apply(all.data[,pos],1,mean,na.rm=TRUE))
	    }
	
	colnames(av.data) <- c("group",sa.group)
	
	write.table(av.data, file=paste(outpath, "/featureGroups_normby_",tile.list[i],"_averaged.csv", sep=""),
	              append = FALSE, quote = TRUE, sep = ",",
	              row.names = FALSE,col.names = TRUE )

    SA.data <- all.data[,2:(ncol(df.SA)+1)] 
    QC.data <- all.data[,(ncol(df.SA)+2):ncol(all.data)] 
  
    if ( any(is.na(SA.data)>0) ) {
      #cat(">> NA are present in SA --> convert to between threshold", "\n")
      cat( ">> NA are present in SA --> convert to between threshold", 
            file= f.info, append=TRUE, sep="\n")

	SA.data[is.na(SA.data)] <- sample(seq.num, length(SA.data[is.na(SA.data)]), replace = TRUE)
      
    }else{}
    
    if ( any(is.na(QC.data)>0) ) {
      #cat(">> NA are present in QC -> convert to between threshold", "\n")
      cat( ">> NA are present in QC -> convert to between threshold",
           file= f.info, append=TRUE, sep="\n")
      
      QC.data[is.na(QC.data)] <- sample(seq.num, length(QC.data[is.na(QC.data)]), replace = TRUE)
      
    }else{}
    
    # PCA with all sample
  res.pca <- prcomp(t(SA.data), 
                    scale = TRUE)
  # predict PCA coordinate
  QC.coord <- predict(res.pca, newdata = t(QC.data))             
  
  ### a) plot PCA
  ###############
  # plot PCA with QC
  
  png(filename =paste(outpath, "/PCA_",tile.list[i],".png", sep=""),
      width = 480, height = 480, units = "px")
  
  p <- fviz_pca_ind(res.pca, 
                        geom="point", 
                        pointsize = 4,
                        title = tile.list[i])
  
  print(p)
      
  print(  fviz_add(p, QC.coord, geom="point", 
               color ="red",
               addlabel = FALSE ) )
  
  dev.off()
  
  
  # B) calculate RSD% for each peak of QC
  sd.QC <- apply(QC.data, 1, sd)
  mean.QC <- apply(QC.data, 1, mean)
  
  RSD <- sd.QC/mean.QC*100
  all.RSD <- cbind(all.RSD,RSD)
  
  ## C) euclidienne distance
  # dist.QC <- stats::dist(as.matrix(QC.data), method = "euclidean",diag = FALSE)
  dist.QC <- vegdist(as.matrix(QC.data), method = "bray")
  dist.QC[lower.tri(dist.QC)] <- NA
  dist.QC[lower.tri(dist.QC, diag = TRUE)] <- NA
  all.dist <- cbind(all.dist, na.omit(reshape2::melt(dist.QC))$value)
  
  ## D) cor matrix
  res <- cor(QC.data)
  res[lower.tri(res)] <- NA
  res[lower.tri(res, diag = TRUE)] <- NA
  all.Rsquare <- cbind(all.Rsquare, na.omit(reshape2::melt(res))$value)
  
}

colnames(all.RSD) <- tile.list
colnames(all.dist) <- tile.list
colnames(all.Rsquare) <- tile.list

#########################
## b) frequency graphe
#########################
# declaring the break points
break_points <-seq(0, 100, by=10)
# transforming the data
myFUN <- function(var1,var2) { cut(abs(var1), var2,right=FALSE) }
                          
data_transform <- apply(all.RSD,2,
                      FUN=myFUN, var2= break_points)

# calculating prop cumulative frequency
cumulative_freq <- apply(data_transform,2,
                         function(x){ c(0, cumsum(prop.table(table(x)))*100) } )  

max_length <- max(sapply(cumulative_freq, length)) # Find max length

# Pad each vector with NA values to make them all the same length
padded_vectors <- lapply(cumulative_freq, function(x) {
  if(length(x) < max_length) {
    c(x, rep(100, max_length - length(x)))
  } else {
    x
  }
})

# Convert the list of padded vectors to a matrix
result_matrix <- do.call(cbind, padded_vectors)

# Convert the matrix to a data frame or table
cumulative_freq <- as.data.frame(result_matrix)

col.plot <- brewer.pal(12, "Paired")

png(filename =paste(outpath, "/cumulative frequency.png", sep=""),
    width = 480, height = 480, units = "px")

# plotting the data
plot(break_points, 
     cumulative_freq[,1],
     las=1, type="b",
     col= col.plot[1],
       bg = col.plot[1],
       pch=21,
     xlab="% RSD",
     ylab="Cumulative Frequency (%)")

s.pch <- c(21,22,23,24,25,21,22,23,24,25,21,22,23)
# creating line graph
for (i in 2:ncol(cumulative_freq))
    {
      points(break_points, 
             cumulative_freq[,i],
             type="b",
             col= col.plot[i],
             bg = col.plot[i],
             pch=s.pch[i]
             )
    }

legend("bottomright",
       legend=colnames(cumulative_freq),
       col= col.plot,
       pt.bg= col.plot,
       pch=s.pch,
       box.lty=0)

dev.off()

################################################################
png(filename =paste(outpath, "/pearson corr.png", sep=""),
    width = 480, height = 480, units = "px")

      boxplot(all.Rsquare,
              las=2)
dev.off()

################################################################
png(filename =paste(outpath, "/bray-curtis distance.png", sep=""),
    width = 480, height = 480, units = "px")

    boxplot(all.dist,
            las=2)
dev.off()

