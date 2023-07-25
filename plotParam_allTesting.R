.singleFoG2 <- function(data, ML2, stand, clearHaloStand, dotedge = 3.4, maxDist = 40, ymax = 200, i, label, plotCompon=FALSE){
  
  startX <- which(data[[i]][,1] > dotedge+0.5)[1]
  stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand
  data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
  
  xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  if (ML2[[i]][length(ML2[[i]])] == 'Log') {
    yy<- .curveLog(xx, ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7])
  } else if (ML2[[i]][length(ML2[[i]])] == 'NegLog') {
    yy <- .curveNegLog2(xx, ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3]) 
  } else {
    yy <- .curveParadox(xx, ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[4], ML2[[i]]$par[6], ML2[[i]]$par[7], ML2[[i]]$par[8])
  }

  #RAD
  ploty <- data[[i]]$x
  ploty[ploty < 0] <-0

  plot(data[[i]]$distance, ploty, cex=0.7, col=grey(0.7), type="p", ylim=c(0, ymax), xlim=c(0, maxDist -dotedge), xaxt="n", yaxt="n", xlab="", ylab="")
  axis(2, labels=FALSE)
  yyplot <- (yy+min(data[[i]]$x))
  yyplot[yyplot < 0] <- 0
  points(exp(xx), yyplot, type="l", col="black", lwd=3)
  
  if(length(xx)<1){
    xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  }
  
  # yy<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], log(xx))
  # yy <- (yy+min(data[[i]]$x))
  # yy[yy < 0] <- 0
  # if (slope >1){
  #   
  #   # points(xx, yy, type="l", col="black", lwd=2)
  #   
  #   if(plotCompon){
  #     xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  #     yy2.1<- .curve(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3],xx)
  #     yy2.2<- .curve(ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7],xx)
  #     yy1plot <- (yy2.1 +min(data[[i]]$x))
  #     yy1plot[yy1plot <0] <-0
  #     yy2plot <- (yy2.2 +min(data[[i]]$x))
  #     yy2plot[yy2plot <0] <-0
  #     points(exp(xx), yy1plot , type="l", col="orange", lwd=2, lty=2)
  #     points(exp(xx), yy2plot, type="l", col="orange", lwd=2, lty=2)
  #   }
  # }
  mtext(label, side=3, cex=0.6)
}

plotParamTest <- function(projectName, ML2, stand,  clearHaloStand, standardLoc = 2.5, ymax=200, dotedge = 3.4, maxDist= 40, xplots = 4, height = 10, width=7, overwrite = TRUE, popUp = TRUE, label=label, savePDF = TRUE, plotSub = plotSub, plotCompon=plotCompon){
  data <- eval(parse(text=projectName))
  if(is.na(plotSub[1])){
    plotSub <- 1:length(data)
  }
  fileFolder <- projectName
  dir.create(file.path(getwd(), "figures"), showWarnings= FALSE)
  dir.create(file.path(getwd(), "figures", fileFolder), showWarnings= FALSE)
  t <- file.path("figures", projectName , paste(projectName, "_test.pdf", sep=""))
  if (!overwrite){
    if (file.exists(t)){
      t <- file.path("figures", projectName , paste(projectName, "_FoG_2_FoG", FoG, "_RAD", RAD, ".pdf", sep=""))
      if (file.exists(t)){
        k <- 2
        while(file.exists(t)){
          k <- k+1
          t <- file.path("figures", projectName, paste(projectName, "_FoG_", k, "_FoG", FoG, "_RAD", RAD, ".pdf", sep=""))
        }
      }
    }
  }
  
  if(xplots > length(plotSub)){
    xplots <- length(plotSub)
  }
  if (ceiling(length(plotSub)/xplots) < 6) {
    yplots<- ceiling(length(plotSub)/xplots)}
  else {yplots<- 6}
  numpages <- ceiling(length(plotSub)/(xplots*yplots))
  if(savePDF){
    pdf(t, width=width, height=height)
  }
  # if(!savePDF){
  # quartz(width=width, height=height)
  # }
  par(mfrow=c(yplots , xplots), mar=c(1,1,1,1), oma=c(4,5,1,1))
  for (k in plotSub){
    .singleFoG2(data = data, ML2 = ML2, dotedge = dotedge, maxDist = maxDist, ymax = ymax, stand = stand, i = k,  clearHaloStand = clearHaloStand, label=label[k], plotCompon=plotCompon)
    if(numpages == 1){
      if (k >= xplots*yplots-xplots+1){
        axis(1, cex.axis=1)
      }
      else {axis(1, cex.axis=1, labels= FALSE)}
    }
    if(numpages == 2){
      if (k >= xplots*yplots-xplots+1 & k < xplots*yplots+1){
        axis(1, cex.axis=1)
      }
      if (k >= 2*xplots*yplots-xplots+1){
        axis(1, cex.axis=1)
      }
      else {axis(1, cex.axis=1, labels= FALSE)}	
    }				
    if(numpages == 3){
      if (k >= xplots*yplots-xplots+1 & k < xplots*yplots+1 | k >= 2*xplots*yplots-xplots+1 & k < 2*xplots*yplots+1 | k >= 3*xplots*yplots-xplots+1){
        axis(1, cex.axis=1)
      }
      else{axis(1, labels=FALSE)}
    }				
    axis(1, labels=FALSE)
    j <- 1
    while (j <= numpages){
      if (k %in% seq(1, j*yplots*xplots, by=xplots)) {axis(2, cex.axis=1, las=2)}
      j <- j+1
    }
  }
  
  mtext("Distance (mm)", outer=TRUE, side=1, line=2, cex=1.2)
  mtext("Pixel intensity", outer=TRUE, side=2, line=2, cex=1.2)
  
  if(savePDF){
    dev.off()	
    cat(paste("\nFigure saved: ", t, sep=""))
    
    if(popUp){
      tt <- paste("open", t)
      system(tt)
    }
  }
}

plotParamTest('speed', speed.ML2, standard, cHS, label = names(speed), plotSub = NA, plotCompon = FALSE)
