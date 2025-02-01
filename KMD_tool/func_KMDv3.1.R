# #################################################################################################################################################
# FUNCTION check and produced subDir folder
###########################################
#February 2017 -- mod on the 2024-11-06 
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

################################################################################
#                                   Round FUNCTION                                  #
################################################################################
round2 <- function(x, n) { # check out https://janajarecki.com/blog/r-does-not-round-2-5-to-3/
  x <- as.numeric(x)
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z * posneg
}

################################################################################
#                                   KMD FUNCTION                                  #
################################################################################
KMD_FUN <-function(mz,ret,int,
                   form_unit= "CF2",
                   hs_tol = 0.005, 
                   int_tresh = 500)
{
  # dependency
  require(MetaboCoreUtils)
  require(data.table)
  
  # Calculate KMD
  dmz <- calculateMass(form_unit)
  km <- mz *round2(dmz,0)/dmz
  kmd <- round2(km,0)- km
  
  # Calculate modulo: Compounds from the same 
  ## homologous series bear an identical modulo
  modulo <- mz %% dmz
  
  df.kmd <- data.frame(mz = mz, ret = ret, int=int, mod = modulo, kmd = kmd)
  df.kmd <- df.kmd[df.kmd$int>=int_tresh,]
  df.kmd <- na.omit (df.kmd)
  df.kmd <- df.kmd[order(df.kmd$mod),]
  
  ## assign all feature to one homologue serie
  HS_num <- rep(0, length(df.kmd$mz) )
  exit_loop <- FALSE
  
  for (n in seq_along(df.kmd$mz) ) {
    if (HS_num[n] == 0) {
      
      HS_num[n] <- n
      i <- n
      
      if (i < length(df.kmd$mz) ) {
        while (abs(df.kmd$mod[i] - df.kmd$mod[n]) < hs_tol && !exit_loop ) {
          HS_num[i] <- n 
          i <- i + 1
          if (i >= length(df.kmd$mz) ){ exit_loop <-TRUE }  } 
      } else {HS_num[i] <- n
      exit_loop <-TRUE}
    }
  }
  
  df.kmd$HS_num <- HS_num
  
  return (df.kmd)
}  

################################################################################
#                                   KMD Filtering                                  #
################################################################################
KMD_filter <-function(df.kmd, 
                      form_unit= "CF2", 
                      thr_dmz = 10,
                      n_min=3)
{
  # Calculate mass unit
  dmz <- calculateMass(form_unit)
  
  ## Filtering by mz and RT 
  ## Iterate through each data point 
  # asum we keep the first point
  df.kmd -> df.kmd.start
  df.kmd.out <-NULL
  
  for (l in 1:length(unique(df.kmd$HS_num)) ) # treat each HS separately
  {
    df.kmd <- df.kmd.start
    
    df.kmd <- df.kmd[df.kmd$HS_num == unique(df.kmd$HS_num)[l], ]# select the HS
    
    if (nrow(df.kmd) < n_min) {}else{
      
      df.kmd <- df.kmd[order(df.kmd$mz),] # order to make the job easier
      
      out.list <- list() # keep all computing
      
      seq <- 1:(nrow(df.kmd)-n_min) # keep pos seq to test
      seq <- seq[seq>0]
      df.kmd -> df.kmd.temp # temp sequence
      
      for (k in seq )
      {
        df.kmd <- df.kmd.temp[k:nrow(df.kmd.temp),]
        
        i=2
        n=2
        exit_loop0 <- FALSE
        while (!exit_loop0 ) { 
          if (i <= n) {
            exit_loop <- FALSE
            i <- n
            j <- 1
            while (!exit_loop ) {
              # Check if the diff. between the current the previous values meet  dy criteria
              if ( (df.kmd$mz[i] - df.kmd$mz[i-1]) < ((j*dmz)-thr_dmz) ) {
                ## if diff lower than threshold --> false positif
                df.kmd <- df.kmd[-i, ] # remove the false positive
                exit_loop <-TRUE
                i <- i-1 # restart at one lower index
                if(i<2){i=2} #in case of the 2sd don't meat criteria
              }else { }
              if ( (df.kmd$mz[i] - df.kmd$mz[i-1]) > ((j*dmz)+thr_dmz) ){ 
                #if diff higher than threshold -->  maybe false positif 
                # check next dy
                j=j+1} else{ 
                  exit_loop <-TRUE
                  i <- i+1 # case of empty data
                }
              # If the point passes the check, add its index to indices_to_keep
              if (n>=nrow(df.kmd)) {exit_loop <-TRUE} else {} # that the end
            } #end while
          } else {n <- n+1 } #end if restart i position 
          if (n>=nrow(df.kmd)) {exit_loop0 <-TRUE} # end of the end
        }
        
        # do th same for remove dRT not pos.
        i=2
        exit_loop <- FALSE
        while (!exit_loop ) { 
          if ( (df.kmd$ret[i] - df.kmd$ret[i-1]) < 0  )
          {
            df.kmd <- df.kmd[-i, ] # remove the false positive
            exit_loop <-TRUE
            i <- i-1 # restart at one lower index
            if(i<2){i=2}
          } else {
            i <- i+1
          }
          
          if (i>=nrow(df.kmd)) {exit_loop <-TRUE} else {}
        }
        
        out.list <-c(out.list,list(df.kmd) )
        
      } #end loop alll option of start
      
      # keep longest serie
      long.list <- unlist(lapply(out.list, function (x) {nrow(x)}))
      
      out.list <- out.list[[which.max(long.list)]]
      
      # check if unique values or duplicate within HS
      out.list <- out.list[order(out.list$mz),] 
      
      unique_values <- numeric(0)# Initialize an empty vector to store unique values
      out.kmd.HS <- NULL
        # Iterate through each value
        for (i in seq_along(out.list$mz)) {
          # Check if the current value is within the buffer range of any existing unique value
          if (!any(abs(unique_values - out.list$mz[i]) <= 0.0005)) {
            # If not, add it to the unique values vector
            unique_values <- c(unique_values, out.list$mz[i])
            out.kmd.HS <- rbind(out.kmd.HS,out.list[i,] )                
          }
        }
        
      df.kmd.out <- rbind(df.kmd.out, out.kmd.HS)
      
    } # end loop HS
    
  }# end loop if enough n HS
  
  # Calculate number of members in each HS
  HS_num_temp <-as.matrix(table(df.kmd.out$HS_num))
  hsnumber <- rep(0, length(df.kmd.out$HS_num) )
  for (i in 1:nrow(HS_num_temp) )
  {
    hsnumber[df.kmd.out$HS_num %in% as.numeric(row.names(HS_num_temp))[i] ]<- HS_num_temp[i,1]
  }
  
  df.kmd.out$Homologues <- hsnumber
  
  ## filter if number of homologues is not greater than specified minimum value
  # once again
  df.kmd.out <- df.kmd.out[df.kmd.out$Homologues>=n_min,] 
  
  # re-organized 
  # Find unique values in the original vector
  un_val <- unique(df.kmd.out$HS_num)
  conv_val <- integer(length(df.kmd.out$HS_num))
  
  # Map each unique value to a consecutive integers
  for (i in seq_along(un_val)) {
    conv_val[df.kmd.out$HS_num == un_val[i]] <- i
  }
  df.kmd.out$HS_num <- conv_val
  
  return (df.kmd.out)
}
