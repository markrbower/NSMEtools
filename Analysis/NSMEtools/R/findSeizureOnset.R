findSeizureOnset <- function(subjectDateDir) {
  library(tidyverse)
  library(dplyr)
  library(bcp)
  library(ggplot2)
  library(cowplot)
  
  setwd(paste0('/Users/markrbower/Dropbox/Documents/Concepts/2022_07_01_NSMEtools/NSMEtools/Presentation/Figures/', subjectDateDir, '/b_stim') )
  
  filenames <- list.files('.','spectrum')
  # load a file to get the frequency information
  sumPower <- 0
  if ( length(filenames) > 0 ) {
    load( filenames[1] )
    idx <- which( spectrum.f > 15 ) # All indices above 15 Hz.
    N <- dim( spectrum.m )[1] # Number of spectrogram slices
    M <- dim( spectrum.m )[2] # Number of frequency bins
    sumPower <- vector( mode='numeric', length=N )
    sumSpec <- spectrum.m * 0
    cnt = 0
    validThetaCount <- 0
    for ( filename in filenames ) {
      cnt = cnt + 1
      print( paste0( cnt, " of ", length(filenames) ) )
      load( filename )
      ss <- sum(sum(spectrum.m))
      theta_filename <- str_replace(filename,'spectrum','1_theta')
      load(theta_filename)
      thetaSum <- sum( thetaResult$valid )
      if ( ss > 500E3 & thetaSum > 10 ) {
        validThetaCount <- validThetaCount + 1
        sumSpec <- sumSpec + spectrum.m
        sumPower <- sumPower + sapply( 1:N, function(j) {sum(spectrum.m[j,idx])} )
        
        idx <- which( spectrum.f > 10 & spectrum.f < 50 )
        sumPowerMid <- sapply(1:N,function(j) {sum(spectrum.m[j,idx])})
        idx <- which( spectrum.f > 7 & spectrum.f < 9 )
        sumPowerTheta <- sapply(1:N,function(j) {sum(spectrum.m[j,idx])})
        ratio <- sumPowerTheta / sumPowerMid
        # image( spectrum.m )
        # plot( sumPowerTheta, main=paste0("theta") )
        # plot( ratio, main=paste0("theta / broadband : ", thetaSum, ' : ', ss) )
      }
    }
#    plot( sumPower )
    par( mfrow= c(3,1) )
    
    print( validThetaCount )
    # Convert matrix to a dataframe
    sumSpec.df <- reshape2::melt(sumSpec, c("x", "y"), value.name = "z")
    # sumSpec.m <- matrix(data=sumSpec,nrow=M,ncol=N,dimnames=list(1:M,1:N))
    #sumSpec.df <- data.frame(x=rep(as.numeric(dimnames(sumSpec.m)[[1]]),each=ncol(sumSpec.m)),y=rep(as.numeric(dimnames(sumSpec.m)[[2]]),times=nrow(sumSpec.m),z=as.vector(sumSpec.m))
    
    p1 <- ggplot(data=sumSpec.df,aes(x=x,y=y,fill=z))+geom_tile()+guides(fill = guide_colourbar(label = FALSE,ticks = FALSE))+theme(legend.position = "none")
    # image( sumSpec )
    
    idx <- which( spectrum.f > 10 & spectrum.f < 50 )
    sumPowerMid <- sapply(1:N,function(j) {sum(sumSpec[j,idx])})
    bcp.sp <- bcp( sumPower, p0=0.01 )
#    plot(bcp.sp$posterior.prob)
    index_prob = cbind(1:N, bcp.sp$posterior.prob)
    colnames(index_prob) = c("index", "prob")
    index_prob = as.data.frame(index_prob) %>% arrange(desc(prob))
    good_idx <- which( index_prob$prob > 0.95 )
    #  head(index_prob)
    seizureStartIdx <- index_prob$index[ min( good_idx ) ]
    sumPowerMid <- sapply(1:N,function(j) {sum(sumSpec[j,idx])})
    idx <- which( spectrum.f > 7.5 & spectrum.f < 8.5 )
    sumPowerTheta <- sapply(1:N,function(j) {sum(sumSpec[j,idx])})
    ratio <- sumPowerTheta / sumPowerMid
    df <- data.frame(index=1:length(ratio),theta=sumPowerTheta,ratio=ratio)
    p2 <- ggplot( df, aes(x=index,y=theta) ) + geom_line()
    p3 <- ggplot( df, aes(x=index,y=ratio) ) + geom_line()
    
    # velocity
    fid <- file( './video_position_smoothed.txt', 'r' )
    t = as.numeric( strsplit(readLines(fid,1),",")[[1]] )
    x = as.numeric( strsplit(readLines(fid,1),",")[[1]] )
    y = as.numeric( strsplit(readLines(fid,1),",")[[1]] )
    close(fid)

    p_out <- plot_grid(p1, p2, p3, align = "v", ncol = 1, axis="l", rel_heights = c(.4, .3, .3))
    print(p_out)
    
  } else {
    print( paste0( "No files found in ", getwd() ) )
  }

}


