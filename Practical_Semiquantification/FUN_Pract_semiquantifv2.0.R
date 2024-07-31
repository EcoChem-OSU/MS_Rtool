#########################################################
## FUNCTION FOR Practical Semi-quantification calibration
########################################################
## version:2.0
## Date: 2024-03-18
## Author: Dunping Cao @ Oregon state University modified by B.Droz 

## Reference to cite
## Cao, D.; et al. J. Am. Soc. Mass Spectrom. 2023, 34 (5), 939–947. 
## https://doi.org/10.1021/jasms.3c00019.

####################Function we use for uncertainty####################################
####################Please do not change these#########################################
pred_int_upper <- function(x_new,alpha,data_set,response_new,index){
  n<- dim(data_set)[1]
  response <- data_set[,2]
  x <- data_set[,1]
  X <- matrix(x,ncol=1)
  ywt <- 1/(x^index)
  oumiga_inverse <- diag(ywt)
  H <- X%*%solve(t(X)%*%oumiga_inverse%*%X)%*%t(X)%*%oumiga_inverse 
  sigma_hat2 <- (matrix(response,nrow=1)%*%oumiga_inverse%*%(diag(1,n)-H)%*%matrix(response,ncol=1))/(n-1)
  X_new <- matrix(x_new,nrow=1)
  cal <- lm(formula = response ~0+ x, weights = ywt)
  upper_bound <- X_new%*%matrix(as.vector(cal[[1]]),ncol=1)+qt(1-alpha/2,n-1)*sqrt(sigma_hat2[1,1]*( (X_new[1,1])^index + X_new%*%solve(t(X)%*%oumiga_inverse%*%X)%*%t(X_new)))
  return(upper_bound-response_new)
}

pred_int_lower <- function(x_new,alpha,data_set,response_new,index){
  n<- dim(data_set)[1]
  response <- data_set[,2]
  x <- data_set[,1]
  X <- matrix(x,ncol=1)
  ywt <- 1/(x^index)
  oumiga_inverse <- diag(ywt)
  H <- X%*%solve(t(X)%*%oumiga_inverse%*%X)%*%t(X)%*%oumiga_inverse 
  sigma_hat2 <- (matrix(response,nrow=1)%*%oumiga_inverse%*%(diag(1,n)-H)%*%matrix(response,ncol=1))/(n-1)
  X_new <- matrix(x_new,nrow=1)
  cal <- lm(formula = response ~0+ x, weights = ywt)
  lower_bound <- X_new%*%matrix(as.vector(cal[[1]]),ncol=1)-qt(1-alpha/2,n-1)*sqrt(sigma_hat2[1,1]*( (X_new[1,1])^index + X_new%*%solve(t(X)%*%oumiga_inverse%*%X)%*%t(X_new)))
  return(lower_bound-response_new)
}

Uncertainty_lower <- function(x_new,alpha,data_set,index){
  n<- dim(data_set)[1]
  response <- data_set[,2]
  x <- data_set[,1]
  ywt <- 1/(x^index)
  cal <- lm(formula = response ~ 0+x, weights = ywt)
  X_new <- matrix(x_new,nrow=1)
  response_new <- (X_new%*%matrix(as.vector(cal[[1]]),ncol=1))[1,1]
  mtry<- try(uniroot(pred_int_upper,c(0,10*max(x)),alpha=alpha,data_set=data_set,response_new=response_new,index=index)$root,silent = TRUE)
  if (class(mtry) != "try-error"){
    uncertainty_lower <- uniroot(pred_int_upper,c(0,10*max(x)),alpha=alpha,data_set=data_set,response_new=response_new,index=index)$root
  }else{
    uncertainty_lower <- 0
  }
  return(uncertainty_lower)
}

Uncertainty_upper <- function(x_new,alpha,data_set,index){
  n<- dim(data_set)[1]
  response <- data_set[,2]
  x <- data_set[,1]
  ywt <- 1/(x^index)
  cal <- lm(formula = response ~ 0+x, weights = ywt)
  X_new <- matrix(x_new,nrow=1)
  response_new <- (X_new%*%matrix(as.vector(cal[[1]]),ncol=1))[1,1]
  mtry<- try(uniroot(pred_int_lower,c(0,10*max(x)),alpha=alpha,data_set=data_set,response_new=response_new,index=index)$root,silent = FALSE)
  if (class(mtry) != "try-error"){
    uncertainty_upper <- uniroot(pred_int_lower,c(0,10*max(x)),alpha=alpha,data_set=data_set,response_new=response_new,index=index)$root
  }else{
    uncertainty_upper <- NaN
  }  
  return(uncertainty_upper)
}

table_uncertainty <- function(x_series,alpha,data_set,index){
  lower_bound <- map_dbl(x_series, ~ Uncertainty_lower(.x, alpha=alpha,data_set=data_set,index=index))
  upper_bound <- map_dbl(x_series, ~ Uncertainty_upper(.x, alpha=alpha,data_set=data_set,index=index))
  return(data.frame(x_series,lower_bound,upper_bound))
}

####################Function we use for uncertainty###################################
####################Please do not change these########################################
###EXTRA FUNCTION
###########################################
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

# creat output
date <- Sys.Date()
output <- creat.subDir(workdir,paste(date,output.folder.name,sep="") )
################################################################################
#########################################
## Cheak and download or releaod package
## v1.0 Emmanuel Rey 2013 EAWAG
########################################
# check library version and load library
check.lib <- function(package.list) {
  
  for(pkg in package.list){
    
    libTest <- try(library(pkg,character.only=TRUE),silent=TRUE)
    
    if(class(libTest)=='try-error'){
      
      updteTest <- try(install.packages(pkg))   
      
      if(class(updteTest)=='try-error'){update.packages(pkg)}
      
      else{install.packages(pkg)}
      
      libTest <- try(library(pkg,character.only=TRUE),silent=TRUE)
      
      print(pkg)
      
      }
    }
  }


