#' Save maximum likelihood output

#' @description Saves the output of maximum likelihood functions - asym, od50, scal, sigma and lnLik.

#' @inheritParams maxLik

#' @return A dataframe "projectName_ML.df" is saved to the global environment and a .csv file "projectName_ML.csv" is exported to the "parameter_files" directory. 

#' @export

#' @author Aleeza C. Gerstein

saveML2Param <- function(projectName){
  #Sys.Date() : asks system about current time and date
  fileFolder <- paste(Sys.Date(), projectName, sep="_")
  newdir <- file.path(getwd(), "parameter_files")
  newdir2 <- file.path(getwd(), "parameter_files", projectName)
  if (!file.exists(newdir)){		
    dir.create(newdir, showWarnings = FALSE)
    cat(paste("\nCreating new directory: ", newdir), sep="")
  }
  if (!file.exists(newdir2)){		
    dir.create(newdir2, showWarnings = FALSE)
    cat(paste("\nCreating new directory: ", newdir2), sep="")
  }
  
  ML2.df <- .ML2param(projectName)
  
  ML2df <- paste(projectName, "_ML2.df", sep="")
  
  filename2 <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ML2.csv", sep=""))	
  
  cat(paste("\n", ML2df, " has been written to the global environment", sep=""))
  assign(ML2df, ML2.df, inherits=TRUE)
  cat(paste("\nSaving files: ", filename2, sep=""))
  
  write.csv(ML2.df, file=filename2, row.names=FALSE)	
}

.ML2param <- function(projectName){
  data <- eval(parse(text=projectName))
  ML2 <- eval(parse(text=paste(projectName, ".ML2", sep="")))
  #lapply: extracts the parametrs of ML2 list
  asymA <- round(unlist(lapply(ML2, function(x) x$par[1])), 2)
  od50A <- round(unlist(lapply(ML2, function(x) x$par[2])), 2)
  scalA <- round(unlist(lapply(ML2, function(x) x$par[3])), 2)
  sigma <- round(unlist(lapply(ML2, function(x) x$par[4])), 2)
  asymB <- round(unlist(lapply(ML2, function(x) x$par[5])), 2)
  od50B <- round(unlist(lapply(ML2, function(x) x$par[6])), 2)
  scalB <- round(unlist(lapply(ML2, function(x) x$par[7])), 2)
  lnLik <- round(unlist(lapply(ML2, function(x) x$lnLik)), 2)
  type <- round(unlist(lapply(ML2, function(x) x[[9]][["type"]])), 2)
  ML2.df <- data.frame(line = names(data), asymA, od50A, scalA, sigma, asymB, od50B, scalB, lnLik)
  return(ML2.df)
}

saveML2Param("nazli3")
