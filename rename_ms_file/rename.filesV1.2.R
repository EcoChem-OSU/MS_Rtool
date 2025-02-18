###################################################
### rename files names from wiff and mZXML by folder
####################################################
## DUMMY TOOL CODE - 
## version:1.2
## Date: 2024-01-12
## Author: Boris Droz @ Oregon State University
###############################################################################
## Description:

## require two folder wiff and mzxml in the workdir folder with the file 
## you want to rename (same number of files)
## the code directly rename it on the folder
## files should be in the same order to rename it
## The code came with two options:

### Opt. 1: rename wiff similar as the mzxml --assuming file are same order
### Opt2. first get name on the wiff and mzxml into a table
## ..., and then rename it on a table in a"new_names" column 

###############################################################################
#### RUN SECTION OF THE CODE STEP BY BLOCK 
###############################################################################

###############################################################################
###SET PARAMETER
#################
wifftype <- ".wiff" # could set ".wiff" or ".wiff2"
outtype <- ".mzXML" #set mzML or mzXML
workdir <- "C:/Users/drozditb/Documents/OSU_data_analysis/name_test"

# workdir <- "D:/Patroon_NTS/RIL"
setwd(workdir)
###############################################################################
## START SCRIPT ## To run anyway
#################
f.fns.wiff <- paste(workdir,"/wiff/",sep="")
f.fns.mzxml <- paste(workdir,"/mzxml/",sep="")
  
fns.wiff <- list.files(path= f.fns.wiff,
                      pattern= paste(wifftype,"$",sep="") ,full.names = FALSE)

fns.mzxml <-list.files(path=f.fns.mzxml,
                       pattern=outtype ,full.names = FALSE)

## to check if number of files are similar
print(fns.wiff); length(fns.wiff)
print(fns.mzxml); length(fns.mzxml)

################################################################################
### Opt. 1: rename wiff similar as the mzxml --assuming file are same order
################################################################################
################################################################################

old.fns <- paste(f.fns.wiff,fns.wiff,sep="")
new.fns <- paste(f.fns.wiff,
                 sapply(strsplit(fns.mzxml,"[.]"),"[[",1),wifftype,sep="")

file.rename( from= old.fns, to= new.fns)


fns.wiff <-list.files(path= f.fns.wiff,
                    pattern=".timeseries.data" ,full.names = FALSE)
if (length(fns.wiff)==0)
    { }else{
      old.fns <- paste(f.fns.wiff,fns.wiff,sep="")
      new.fns <- paste(f.fns.wiff,
                      sapply(strsplit(fns.mzxml,"[.]"),"[[",1),".timeseries.data",sep="")
    
      file.rename( from= old.fns, to= new.fns)
      }



fns.wiff <-list.files(path= f.fns.wiff,
                     pattern=".wiff.scan" ,full.names = FALSE)
if (length(fns.wiff)==0)
    { }else{
    
    old.fns <- paste(f.fns.wiff,fns.wiff,sep="")
    new.fns <- paste(f.fns.wiff,
                     sapply(strsplit(fns.mzxml,"[.]"),"[[",1),".wiff.scan",sep="")
    
    file.rename( from= old.fns, to= new.fns)
    }

################################################################################
## Opt2. first get name on the wiff and mzxml into a table
## ..., and then rename it on a table in a"new_names" column
################################################################################
################################################################################

df.fn.names <- data.frame(cbind(wiff_name= sapply(strsplit(fns.wiff,"[.]"),"[",1),
                            mzxml_name= sapply(strsplit(fns.mzxml,"[.]"),"[",1) ))

                          
# save it and change what you whant in excel or other
write.csv(df.fn.names,"new_sample_namelist.csv",row.names = FALSE)
###############################################################################
###############################################################################
# modify all filename according to a created "new_names" column
###############################################################
#### DONT put any extension file on the new name list #########
###############################################################
df.filename <- read.csv("new_sample_namelist.csv",header=TRUE); head(df.filename)

## rename wiff
old.fns <- paste(f.fns.wiff, paste(df.filename$wiff_name,wifftype, sep=""), sep="" )
new.fns <- paste(f.fns.wiff, paste(df.filename$new_names,wifftype, sep=""), sep="" )
                 
file.rename( from= old.fns, to= new.fns)

# ## rename timesseries
# fns.wiff <-list.files(path= f.fns.wiff,
#                       pattern=".timeseries.data" ,full.names = FALSE)
fns.wiff <- paste(f.fns.wiff, paste(df.filename$wiff_name,".timeseries.data" , sep=""), sep="" )

if (length(fns.wiff)==0)
    { }else{
    old.fns <- fns.wiff #paste(f.fns.wiff,fns.wiff,sep="")
    new.fns <- paste(f.fns.wiff, paste(df.filename$new_names,
                                       ".timeseries.data", sep=""), sep="" )
    
    file.rename( from= old.fns, to= new.fns)
    }
## rename wifff.scan
# fns.wiff <-list.files(path= f.fns.wiff,
#                       pattern=".wiff.scan" ,full.names = FALSE)
fns.wiff <- paste(f.fns.wiff, paste(df.filename$wiff_name,".wiff.scan" , sep=""), sep="" )


if (length(fns.wiff)==0)
    { }else{
    old.fns <- fns.wiff# paste(f.fns.wiff,fns.wiff,sep="")
    new.fns <- paste(f.fns.wiff, paste(df.filename$new_names,
                    ".wiff.scan",sep=""), sep="" ) 
    
    file.rename( from= old.fns, to= new.fns)
    }

## rename the mzXML files
old.fns <- paste(f.fns.mzxml, paste(df.filename$mzxml_name,outtype, sep=""), sep="" )
new.fns <- paste(f.fns.mzxml, paste(df.filename$new_names,
                                    outtype,sep=""), sep="" ) 

file.rename( from= old.fns, to= new.fns)
  
