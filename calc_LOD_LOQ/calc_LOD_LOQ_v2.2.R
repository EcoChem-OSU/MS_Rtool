###############################################################################
## Limit of Detection (LOD) and Limit of quantitfication (LOQ) method 
###############################################################################
# v1.0 Trever
# v2.2 mod B.Droz - 2024-01-22 @ Oregon State University

## Description:
## Input should be a csv file with special character (e.g., - % ? ! /)
## Exemple of input and output ... be sure to have similar structure file

####################################################################
## Experimental aspect : (---- IMPORTANT TO SET YOUR EXPERIMENT ---)
####################################################################
## - Area count (or intensity) should be evaluate in function of the concentration
## - Vial(1999) recommend to have at least 5 data points to make the calculation
## - data should be spread within 200 fold (e.g., 1 to 200 or 200 to 40000) 
##    to follow paper experimental setting, data closer to LOD are better.
## - Ideally one or two point are lower than the expected LOD
## - Do not include conc 0 and area 0 into the calculation
## - Data point do necessarly need to be equal separate (weight normalization)

############
# Reference
## Vial, J.; Jardy, A. Experimental Comparison of the Different Approaches 
##              To Estimate LOD and LOQ of an HPLC Method. 
##              Anal. Chem. 1999, 71 (14), 2672-2677. DOI: 10.1021/ac981179n.
################################################################################

##############################
## PARAMETER --NEED TO BE MODIFIED IF NEEDED
#############
workdir <- "C:/Users/drozditb/Documents/OSU_data_analysis/"

setwd(workdir)

df <- read.table(file.choose(),sep= ",", # select data
                 header=TRUE,na.strings = "NA" ) # not number should be NA
head(df)

# choose the column for concentration and for species to measured
# not necessary if you have similar input than the example
nx <- 1 ## concentration 
ny <- 2:ncol(df)  ##species list

plot.opt <- "YES" # do you want a plot plotted in you screen? YES or NO
set.thres <- 200  # set the threshold on how many fold the range of data is accepted.
################################################################################
###############################################################################

## DO NOT MODIFIED SOMETHING BELOW JUST RUN THE CODE ####

################################################################################
###############################################################################
### LOQ LOD FUNCTION
####################

fun_LODQ <- function(data, x.label=x,y.label=y, plot.check="YES") 
  ## x refer to conc
  ## y refer to area or intensity
  {
  #x <- (c(0.5, 1, 5, 10, 25, 50, 100))
  x <-data[,colnames(data)==x.label]
  y <-data[,colnames(data)==y.label]
  
  #fit simple linear regression model
  model <- lm(y~x)
  # #define weights 
  wt <- 1 / lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2
  
  #wt <- length(x)*(1/x)/(sum(1/x))
  
  #perform weighted least squares regression
  wls_model <- lm(y ~ x, weights=wt)
  
  if (plot.check=="YES") {
      plot(y ~ x, data=data,# creat a scater plot
         xlab= x.label,
         ylab= paste("area", y.label) ) 
      abline(wls_model) #plot the model
  }
  
  RSS <- sqrt( sum((wls_model$residuals)^2) / (nrow(data)-2) )
  
  # inverse 
  LOD <- abs(RSS-wls_model$coefficients[1]) /wls_model$coefficients[2]
  LOQ <- 3*LOD
  
  df.out <- c(LOD = LOD,LOQ = LOQ,
    rsquare = summary(wls_model)$r.squared,
   RMSE= sqrt(mean(wls_model$residuals^2)))
  
  return(df.out)
}

###############################################################################
##############################################################################
######## CODE START HERE
#########################
# rename if name start by X -- 
names(df) <- sub("^X", "", names(df))

d.out <- NULL
names.out <- NULL

f.info <- "consol_log.txt"

for (i in ny)
    {
      # select data to follow Vial methodological plan .....
      cat("RUN LOD LOQ CALC FOR", names(df)[i],"...", "\n")  
      cat("RUN LOD LOQ CALC FOR", names(df)[i],"...", "\n",file= f.info,append=TRUE, sep="\n")
  
      df[,i] <- replace(df[,i], which(df[,i] < 0), NA) # replace neg values by NA
      
      df.in <- na.omit(df[,c(nx,i)])
      df.in <- df.in[!df.in[,1]<=0 | !df.in[,2]<=0,] # remove line with 0 0 data
      
      # select data with measurement not egal to zero keep only one
      if (any(df.in[,2]==0)){
           pos <-max(which(df.in[,2]==0))
           df.in <- df.in[pos:nrow(df.in) ,]
      }else{  }
      
      h.r <- set.thres* min(df.in[!df.in[,2]==0,1]) # higher setting
      
      df.in <- df.in[df.in[,1]<=h.r ,]
      
      if ( nrow(df.in)<5 ) {
        cat("Less then 5 point with pos data ... not enough to run LOD...", "\n")
        cat("Less then 5 point with pos data ... not enough to run LOD...", "\n", 
              names(df)[i],"...", "\n",file= f.info,append=TRUE, sep="\n")
      }else{
        names.out <- c(names.out, names(df)[i])
        d.out<-rbind(d.out, fun_LODQ ( df.in, 
                                       x.label=names(df)[nx],
                                       y.label=names(df)[i],
                                       plot.check=plot.opt ) )
      }
    }
 
row.names(d.out) <- names.out

write.csv(d.out, "output_LOD_LOQ.csv")
