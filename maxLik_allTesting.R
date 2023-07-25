maxLikTest <- function(projectName, standType ="one", clearHalo, diskDiam = 6, standardLoc = 2.5, maxDist=25, ymax=200, xplots = 4, height = 8,  width = 8, needML = TRUE, popUp = TRUE, nameVector = TRUE, overwrite = TRUE, plotParam = TRUE, savePDF= TRUE, plotSub = NA, plotCompon=FALSE){
  options(warn=-1)
  
  data <- eval(parse(text=projectName))
  
  dotedge <- diskDiam/2+0.7
  

    dotMax <- max(sapply(data, function(x) {x[which(x[,1] > standardLoc)[1], 2]}))
    standard <-c( sapply(data, function(x) {dotMax-x[which(x[,1] > standardLoc)[1], 2]}))
    
    if(needML){
      ML2 <- lapply(c(1:length(data)), .getstats2Log, data=data, dotedge=dotedge, maxDist=maxDist, stand=standard, maxSlope=20)
      assign(paste(projectName, ".ML2", sep=""), ML2, inherits=TRUE)
      cat(paste("\n", projectName, ".ML2 has been written to the global environment\n", sep=""))
    }
  
  alarm()
}

.curveLog <- function(x, asym, od50, scal, asymB, od50B, scalB) { asym*exp(scal*(x-od50))/(1+exp(scal*(x-od50)))+asymB*exp(scalB*(x-od50B))/(1+exp(scalB*(x-od50B)))}

.curveNegLog <- function(x, asym, od50, scal, asymB, od50B, scalB) { asym*exp(-1 * scal*(x-od50))/(1+exp(-1 * scal*(x-od50))) - asymB*exp(-1 * scalB*(x-od50B))/(1+exp(-1 * scalB*(x-od50B)))}

.curveParadox <-  function(x, slope, height, shift, drop, asym, midpoint, scal) {
  
  a <- exp(slope * (x - shift))
  b <- 1 - 4 * exp(slope * (x - shift))
  c <- exp(2 * slope * (x - shift))
  d <- (1 + exp(slope * (x - shift))) ^ 4
  
  drop * a * (b + c) / d + height + asym*exp(scal*(x-midpoint))/(1+exp(scal*(x-midpoint)))
}

.getstats2Log <- function(i, data, stand, dotedge=dotedge, maxDist=maxDist, maxSlope=100){
  cat(".")
  
  # data <- speed
  # stand <- standard
  # dotedge <- 3.7
  # maxDist <- 40
  # maxSlope <- 100
  # i <- 5
  
  startX <- which(data[[i]][,1] > dotedge+0.5)[1]
  stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  data[[i]] <- subset(data[[i]], data[[i]]$x != "NA")
  data[[i]]$x <- data[[i]]$x+ stand[i] -min(data[[i]]$x+stand[i])  #the model only fits when it goes down to 0
  data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
  data[[i]]$distance <- log(data[[i]]$distance)
  
  sumsquares.fit <- function(theta){
    asym<-theta[[1]]
    od50<-theta[[2]]
    scal<-theta[[3]]
    sigma<-theta[[4]]
    asymB<-theta[[5]]
    od50B<-theta[[6]]
    scalB<-theta[[7]]
    y<-data[[i]]$x
    x<-data[[i]]$distance
    res <- dnorm(y, (asym*exp(scal*(x-od50))/(1+exp(scal*(x-od50)))+asymB*exp(scalB*(x-od50B))/(1+exp(scalB*(x-od50B)))), sigma, log= T)
    sum(res)
  }
  
  sumsquares.fit2 <- function(theta){
    asym<-theta[[1]]
    od50<-theta[[2]]
    scal<-theta[[3]]
    sigma<-theta[[4]]
    asymB<-theta[[5]]
    od50B<-theta[[6]]
    scalB<-theta[[7]]
    y<-data[[i]]$x
    x<-data[[i]]$distance
    res <- dnorm(y, (asym*exp(-1 * scal*(x-od50))/(1+exp(-1 * scal*(x-od50))) - asymB*exp(-1 * scalB*(x-od50B))/(1+exp(-1 * scalB*(x-od50B)))), sigma, log= T)
    sum(res)
  }
  
  sumsquares.fit3 <- function(theta){
    slope<-theta[[1]]
    height<-theta[[2]]
    shift<-theta[[3]]
    drop<-theta[[4]]
    sigma<-theta[[5]]
    asym<-theta[[6]]
    midpoint<-theta[[7]]
    scal <- theta[[8]]
    y<-data[[i]]$x
    x<-data[[i]]$distance
    
    a <- exp(slope * (x - shift))
    b <- 1 - 4 * exp(slope * (x - shift))
    c <- exp(2 * slope * (x - shift))
    d <- (1 + exp(slope * (x - shift))) ^ 4
    
    res <- dnorm(y, drop * a * (b + c) / d + height + asym*exp(scal*(x-midpoint))/(1+exp(scal*(x-midpoint))), sigma, log= T)
    sum(res)
  }
  
  lowOD <- min(data[[i]]$x)
  highOD <- quantile(data[[i]]$x, 0.99, na.rm = TRUE)
  lowerLog <- c(0, 0, 0,0, 0, 0, 0)
  upperLog <- c(highOD, log(maxDist), maxSlope, 10, highOD,  log(maxDist), maxSlope)
  lowerParadox <- c(0, 0, 0, 0, 0, 0, 0, 0)
  upperParadox <- c(maxSlope, highOD, maxDist, 1000, 10, highOD, maxDist, maxSlope)

  par.tryLog <-c(asym = 0.7*highOD, od50 = log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2, asymB = 0.7*highOD, od50B = log(maxDist)/4, scalB = maxSlope*0.01)
  par.tryNegLog <- c(asym = 0.7*highOD, od50 = log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2, asymB = 0.7*highOD, od50B = log(maxDist)/4, scalB = maxSlope*0.01)
  par.tryParadox <- c(slope = 0.1, height = 0.8*highOD, shift = 0, drop = 200, sigma = 0.2, asym = 0.8*highOD, midpoint = 0, scal = 0)
  
  mlpoint<-c()
  mlpointLog<-find.mle(sumsquares.fit,par.tryLog, method="subplex",upper=upperLog,lower=lowerLog, control=list(maxit=50000))
  mlpointNegLog<-find.mle(sumsquares.fit2,par.tryNegLog, method="subplex",upper=upperLog,lower=lowerLog, control=list(maxit=50000))
  mlpointParadox<-find.mle(sumsquares.fit3,par.tryParadox, method="subplex",upper=upperParadox,lower=lowerParadox, control=list(maxit=50000))
  
  mlpoint <- if (mlpointLog$lnLik > mlpointNegLog$lnLik) mlpointLog else mlpointNegLog
  mlpoint <- if (mlpointParadox$lnLik > mlpoint$lnLik) mlpointParadox else mlpoint
  
  df1 <- c(log = mlpointLog$lnLik, negLog = mlpointNegLog$lnLik, para = mlpointParadox$lnLik)
  df2 <- c('Log', 'NegLog', 'Para')
  maxi <- df2[which.max(df1)]
  
  mlpoint <- append(mlpoint, maxi)
  names(mlpoint[[length(mlpoint)]]) <- 'type'
  
  cat(i, ': ', maxi, '\n')

  cat('Log: ')
  cat(mlpointLog$lnLik, '\n')
  cat('Neg Log: ')
  cat(mlpointNegLog$lnLik, '\n')
  cat('Paradox: ')
  cat(mlpointParadox$lnLik, '\n\n')
  
  mlpoint
}
