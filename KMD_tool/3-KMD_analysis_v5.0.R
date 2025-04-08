##################################################
## Title: KMD Kendrick Mass Defect Calculation and Plot
##################################################
## version:4.0 batch sample
## Date: 2024-03-05
## Author: Amergin McDavid and Boris Droz @ Oregon state University 
## Modified by Peter Bright 2024-07-24
###############################################################################
## libraries and setup (See bottom of code for more first time set up instructions)
library(MetaboCoreUtils)
library(data.table)
library(magrittr)
library(purrr)
library(dplyr)
library(glue)
workdir <- getwd()
###Run the KMD function stored in separate R script called "func_KMDv3.1.R"
source(paste(workdir,"/func_KMDv3.1.R",sep=""))
###############################################################################
## Parameters to Change Here
## Input data
# sample.list <- "sample_list_KMD.csv"
sample.list <- "sample_list_KMD.csv"
feat.files <- "featureGroups.txt"# is an input but is on an output folder

# suspect list to match features to
fn.susp.list <- glue("{workdir}/input/KMD_CF2_neg_SuspectList_2024-07-19.csv")
seed = 2 #Set a number here if the labelling on the output graph sucks. Random, but reproducible for a given number
#################################################
# Mass defect filtering based on the suspect list
#################################################
MD.minmax <-c(-0.49,0.50) # based on the merge suspect list (OCDE,NIST,NORMAN,2EPA)
#MD.minmax <-c(-0.25,0.1) # from Zwiener paper based on OECD suspect list

################################################
## Kendrick Mass Defect (KMD) function parameter
################################################

form_unit <- "CF2" #base unit C2F4 CF2O C2F4O
mz_acc <- 4 # mz accurency (ppm) of the HR MS used to set 
             # Mass tolerance for homologous series (Da)
              # capture what is lower than this diff
              # high value increase the number of peak in the same HS
              # 0.005 Da is the default parameter in FindPFAS
              # should be smaller than Condition for suspect match list
int_tresh <- 1000 # intensity or area count threshold 
thr_dmz = 8 # threshold window for mz filtering
              # error associate to each mz form same HS serie
n_min <- 3 # number minimal to consider one HS

##################################
# Condition for suspect match list
##################################
tolerance <- 8 # 'mz' tolerance (ppm)

################################################################################
#                                   Start Script                                  
################################################################################
setwd(workdir)

# creat output
date <- Sys.Date()
output <- creat.subDir(paste(workdir,"/output",sep=""), paste(date,"_KMD_result_", form_unit, sep=""))

################################################################################
# read data
df.samp <- read.csv(paste(workdir,"/input/",sample.list,sep=""), sep=",",header=TRUE)
df <- read.table(paste(workdir,"/output/",feat.files,sep=""), header=TRUE)
df.susp.list <- read.csv(fn.susp.list)

# select data
p.samp <- df.samp$filename[df.samp$sampletype=="SA"]
# unified names
names(df) <- sub("^X", "", names(df))# rename if name start by X -- 
names(df) <- sub("[.]", "-", names(df))
################################################################################
## write the option of the run
################################
f.info <- paste(output,"/AA_INFO_PARA.txt",sep="")

cat("######################################################",file= f.info,append=TRUE, sep="\n")
cat( paste("*** Kendrick Mass defect model  ---", Sys.Date()), file= f.info, sep="\n")
cat("######################################################",file= f.info,append=TRUE, sep="\n")
cat("R-script KMD_analysis_v3.1 - Droz 2024",file= f.info,append=TRUE, sep="\n")
cat("######################################################",file= f.info,append=TRUE, sep="\n")
cat(paste("Mass defect filtering from", MD.minmax[1],"to", MD.minmax[2]),file= f.info,append=TRUE, sep="\n")
cat("Sample list analyzed",file= f.info,append=TRUE, sep="\n")
cat("--------------------",file= f.info,append=TRUE, sep="\n")
cat(p.samp,file= f.info,append=TRUE, sep="\n")
cat("##############",file= f.info,append=TRUE, sep="\n")
cat("KMD PARAMETER",file= f.info,append=TRUE, sep="\n")
cat("##############",file= f.info,append=TRUE, sep="\n")
cat(paste("KMD on base unit:", form_unit),file= f.info,append=TRUE, sep="\n")
cat(paste("Homologue serie (HS) tolerance (Da):", (mz_acc /1000)),file= f.info,append=TRUE, sep="\n")
cat(paste("Intensity count threshold:", int_tresh),file= f.info,append=TRUE, sep="\n")
cat(paste("Error window in mz (ppm) associate to the unit:", thr_dmz),file= f.info,append=TRUE, sep="\n")
cat(paste("Number minimal to consider one HS:", n_min),file= f.info,append=TRUE, sep="\n")
cat("#####################",file= f.info,append=TRUE, sep="\n")
cat("Suspect List Matching",file= f.info,append=TRUE, sep="\n")
cat("#####################",file= f.info,append=TRUE, sep="\n")
cat(paste("Suspect list used:",basename(fn.susp.list)),file= f.info,append=TRUE, sep="\n")
cat(paste("mz tolerance (ppm):", tolerance),file= f.info,append=TRUE, sep="\n")

for (nd in 1: length(p.samp))
    {
################################################################################
## Kendrick mass defect (KMD) analysis 
#######################################

## filter the data using mass defect
MD_range= MD.minmax ## range based on suspect list values
MD <- df$mz - round2(df$mz,0)    # mass defect of each feature
df <- df[MD>=MD_range[1]&MD<=MD_range[2],]   

dmz <- calculateMass(form_unit)
  
df.kmd <- KMD_FUN (df$mz,
                   df$ret,
                   df[,names(df)==p.samp[nd]],
                  form_unit= form_unit,
                  hs_tol = (mz_acc /1000),
                  int_tresh = int_tresh)

write.table(df.kmd, file=paste(output,"/",p.samp[nd],"_KMD_",form_unit,"nofilt.csv",sep=""),sep=",", 
            append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)


  if ((any(table(df.kmd$HS_num) > 3))) {

  df.kmd <- KMD_filter(df.kmd,
                form_unit= form_unit, 
                thr_dmz = thr_dmz, 
                n_min= n_min)

  ###############################################################################
  ## Match with suspect list
  ##################################################################################################original

  # Initialize matched_rows with the correct number of columns and column names
  matched_rows <- cbind(df.kmd, matrix(NA, nrow = nrow(df.kmd), ncol = ncol(df.susp.list)))
  names(matched_rows) <- c(names(df.kmd), names(df.susp.list))

  # Rename specific columns if needed
  names(matched_rows)[c(4, 5)] <- c("exp_modulo", "exp_KMD")

  # Ensure tolerance is defined
  if (missing(tolerance) || is.na(tolerance)) {
    stop("Tolerance value is missing or NA")
  }

  # Loop through 'mz' values from the suspect list
  for (i in 1:nrow(df.susp.list)) {
    mz_val <- df.susp.list$MONOISOTOPIC_MASS[i]
    
    # Calculate score
    cond1 <- df.kmd$mz - mz_val
    df.kmd$score <- sqrt(cond1^2)
    
    # Filter out rows with NA in mz or score
    df.kmd <- df.kmd[!is.na(df.kmd$mz) & !is.na(df.kmd$score), ]
    
    # Calculate condition
    df.kmd$cond <- abs(df.kmd$mz - mz_val) <= tolerance / 1000
    
    # Check if any condition is met
    if (any(!is.na(df.kmd$cond)) && any(df.kmd$cond, na.rm = TRUE)) {
      min_row_pos <- which.min(df.kmd$score[df.kmd$cond])
      min_row_pos <- which(df.kmd$cond)[min_row_pos]
      
      # Ensure the row indices are within bounds
      if (min_row_pos <= nrow(matched_rows)) {
        # Debugging prints to verify indices
        start_col <- ncol(df.kmd) -1 ####Change here if columns names and values are mismatched in output
        end_col <- (ncol(df.kmd) + ncol(df.susp.list)) -2 ####Change here if columns names and values are mismatched in output
        # print(paste("Matching row:", min_row_pos))
        # print(paste("Assigning columns from", start_col, "to", end_col))
        
        # Ensure correct number of columns are assigned
        if ((end_col - start_col + 1) == ncol(df.susp.list)) {
          matched_rows[min_row_pos, start_col:end_col] <- df.susp.list[i, ]
        } else {
          stop("Column index mismatch: check the indices.")
        }
      }
    }
  }
  # Define a function to check if a vector is strictly increasing
  increasing_RT <- function(x) {
    all(diff(x) > 0)
  }
  increasing_MZ <- function(x) {
    all(diff(x) >= 12.0) #Homologous series must lose a carbon, 12
  }

  # Filter dataframe
  matched_rows <- matched_rows %>%
    group_by(HS_num) %>%
    filter(increasing_RT(ret)) %>%
    filter(increasing_MZ(mz)) %>%
    filter(all(is.na(residual) | residual == first(na.omit(residual)))) %>%
    ungroup()%>%
    mutate(HS_num = as.integer(factor(HS_num)))

  Column_order = c( "int", "exp_modulo", "INCHIKEY", "FORMULA",
    "SMILES",	"NumHDonors",	"NumHAcceptors",	"MD",	"KMD", "modulo", "parent",	"rep_unit",	"rep_mass",	"rep_num",
    "residual", "residual_charge",	"residual_mass","MONOISOTOPIC_MASS", "TYPE",  "LIST_ID","HS_num", "Homologues","exp_KMD", "mz", "ret","PREFERRED_NAME", "ACRONYM")
  matched_rows <- matched_rows[ ,Column_order]

if (nrow(matched_rows) < 2) {
  next  # Skip the current iteration if matched_rows is empty
}

  write.table(matched_rows, file=paste(output,"/",p.samp[nd],"_KMD_",form_unit,"matchsus.txt",sep=""),sep="\t", 
              append=FALSE, row.names=FALSE,col.names=TRUE, quote=TRUE)

  # ################################################################################
  # ## Plot result and save result
  # ##############################
  # Initialize and save the file as a png
  png(filename = paste(output,"/",p.samp[nd],"_KMD_",form_unit ,"_plot.png",sep=""),width = 480, height = 480 )
  par(mar=c(5, 4, 2, 5) )

  x_min <- min(matched_rows$mz, na.rm = TRUE)
  x_max <- max(matched_rows$mz, na.rm = TRUE) 
  if (is.na(x_min) || is.na(x_max) || x_min == Inf || x_max == -Inf) {
    x_min <- 0
    x_max <- 1000
    }

    # Set the y-axis limits to capture the data, but still have a reasonable range (0.2)
  y_min <- min(matched_rows$exp_KMD, na.rm = TRUE)
  y_max <- max(matched_rows$exp_KMD, na.rm = TRUE)
  if (is.na(y_min) || is.na(y_max) || y_min == Inf || y_max == -Inf) {
   y_min <- 0  # Default values when all exp_KMD values are NA
   y_max <- 0.2
  } else {
  if (y_max - y_min < 0.2) {
  y_min <- y_min - 0.05
  y_max <- y_max + 0.05
  } }

  
  plot(NULL,NULL,
      xlab="m/z",
      ylab= paste("KMD ",form_unit,sep=""),
      las=1,
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max))

  pal <- colorRamp(c("blue", "green", "orange", "red"))    # 1) choose colors
  col <- rgb(pal((matched_rows$ret - min(matched_rows$ret)) / diff(range(matched_rows$ret))), max=255)  # 2) interpolate numbers for RT color scaling
  col.ramp <- rgb(pal(seq(0, 1, length.out = 20)), maxColorValue = 255)

  pt <-rep(c(21,22,23,24,25),round2(length(unique(matched_rows$HS_num)),0)/5+1) # point

  #### Set the Target, Suspect, and Unknown series to have colored lines and point outlines
  col.l <- ifelse(is.na(matched_rows$TYPE), "red", #marker color
                  ifelse(matched_rows$TYPE == "Target", "blue",
                  ifelse(matched_rows$TYPE == "Suspect", "black", "red")))
          
  lwd.s <- ifelse(is.na(matched_rows$TYPE), 1.5,  # Handle NA cases (linewidth)
                  ifelse(matched_rows$TYPE == "Suspect", 1.5,
                  ifelse(matched_rows$TYPE == "Target", 1.5, 1.5)))


  # Set a random seed for reproducibility of label placement
  set.seed(seed)
  num_unique_hs <- length(unique(matched_rows$HS_num))
  if (num_unique_hs > 0) { label_offsets <- runif(num_unique_hs, min = 0.00, max = (max(matched_rows$exp_KMD, na.rm = TRUE) / 3))  # Random offsets
 } else {
  label_offsets <- rep(0, num_unique_hs)}
 
### A loop for plotting each point, line, and HS number
  for (i in 1:length(unique(matched_rows$HS_num))) {
      # Get Coordinates for the current HS_num
      current_data <- matched_rows[matched_rows$HS_num == i, ]
      
      # If current_data has less than 1 row, skip this iteration
      if (nrow(current_data) < 1) {
          next
      }
      
      # Calculate the position for the label
      last_point <- tail(current_data, 1)  # Tail = last point, usually is as neat as possible
      x_pos <- last_point$mz
      y_pos <- last_point$exp_KMD + 0.008  # Offset label point by this much

      # Label the series with the HS Number
      text(x=x_pos, 
          y=y_pos,
          labels=i, 
          pos=3,  # Position above the point
          font=2, 
          cex=1.0,  # Adjust text size as needed
          col=col.l[matched_rows$HS_num == i], 
          bg="white")  # Background color for readability

      # Draw line between point and the HS number label
      segments(x0=last_point$mz, 
              y0=last_point$exp_KMD, 
              x1=x_pos, 
              y1=y_pos, 
              col="black", 
              lty=1)  # Line type for connection

      # Plot the line between the actual data points in the HS series
      lines(current_data$mz,
            current_data$exp_KMD,
            type="b",
            col=col.l[matched_rows$HS_num == i],
            lwd=lwd.s[matched_rows$HS_num == i],
            bg=col[matched_rows$HS_num == i],
            pch=pt[i],
            lty=2)
  }

  legend.col <- function(col, lev, ylabel=" RT (min)"){
    
    opar <- par
    n <- length(col)
    bx <- par("usr")
    
    box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
                bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
    box.cy <- c(bx[3], bx[3])
    box.sy <- (bx[4] - bx[3]) / n
    
    xx <- rep(box.cx, each = 2)
    
    par(xpd = TRUE)
    for(i in 1:n){
      
      yy <- c(box.cy[1] + (box.sy * (i - 1)),
              box.cy[1] + (box.sy * (i)),
              box.cy[1] + (box.sy * (i)),
              box.cy[1] + (box.sy * (i - 1)))
      polygon(xx, yy, col = col[i], border = col[i])
      
    }
    par(new = TRUE)
    plot(0, 0, type = "n",
        ylim = c(min(lev), max(lev)),
        yaxt = "n", ylab = "",
        xaxt = "n", xlab = "",
        frame.plot = FALSE)
    axis(side = 4, las = 2, tick = TRUE)
    mtext(ylabel, side=4, line=3, cex.lab=1,las=0)
    
    par <- opar
  }
  legend.col(col = col.ramp, lev = matched_rows$ret/60, ylabel=" RT (min)")
  box()

  dev.off() }


else {
    write.table(df.kmd, file=paste(output,"/",p.samp[nd],"_KMD_",form_unit,"_no_H_series.txt",sep=""), 
                append = FALSE, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = TRUE )
              
    print(glue("\n\nThe {p.samp[nd]} kmd dataframe has no homologous series detected in the sample.\nView the df.kmd dataframe through the terminal by typing 'df.kmd' directly or view the no_H_series.csv for more information."))
  }
    }
print('===============Code Finished===============')
## Description:
###############
## Performed KMD from a feature list
## modified version of the KMD function python code https://github.com/JonZwe/PFAScreen
# ---- first need to create two folder named: "input" and "output" in your working directory

## Reference
## Zweigle, J.; et al. Anal. Bioanal. Chem. 2023. https://doi.org/10.1007/s00216-023-05070-2.

###############################################################################
### https://rdrr.io/github/rformassspectrometry/MetaboCoreUtils/
############ KMD - 
# install.packages("remotes")
# remotes::install_github("rformassspectrometry/MetaboCoreUtils", force = TRUE)