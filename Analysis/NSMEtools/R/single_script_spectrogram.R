single_script_spectrogram <- function() {
library(foreach)
library(here)
library(parallel)
library(doParallel)
library(multitaper)

stepSeconds <- 5
widthSeconds <- 15
idiv_overlap <- function( start=1, stop=100, step=5, width=10 ) {
  stop_limit <- stop - width + 1
  starts  <- seq(from=start, to=stop_limit, by=step)
  windows <- sapply(1:length(starts),function(x) {starts[x] + 1:width - 1})
}

loess.filter <- function(x,span,rn) loess(formula=paste(x,"rn",sep="~"),data=df,span=span)$fitted


output_dir <-  "ResultsAndFigures"
sps = 32000
filedir = "/home/bm662/project/Data/Stim/Stim5/10.19.20/Test9aShamStim/2020-10-19_12-39-36/"
ncs_files <- list.files(filedir, pattern='ncs$')
mef_filedir = ( "/home/bm662/Documents/Concepts/NSMEtools/Analysis/NSMEtools/R/mef2/" )
spectrumMaxFreq <- 50  # In Hertz
password <- ''
mef_files <- list.files(mef_filedir, pattern='mef$')
mef_file  = mef_files[1]
iterationCount <- 0
mef_fileroot <- tools::file_path_sans_ext( mef_file )
mef_fullpath <- paste0( mef_filedir, mef_file )
info <- meftools::mef_info( c( mef_fullpath,password) )
iter_cont <- MEFcont( filename=mef_fullpath, password=password)
mi <- iter_cont$nextElem()
md <- mi$hasNext()
mef_data <- mi$nextElem()
windows <- idiv_overlap(start=1,stop=length(mef_data),step=stepSeconds*sps,width=widthSeconds*sps)
M <- ncol(windows)
tmp <- spec.mtm( ts( mef_data[windows[,1]], frequency = 32000 ), nw=4.0, k=7, jackknife = TRUE, plot=FALSE )
spectrumLimit <- max(which(tmp$freq<spectrumMaxFreq))
spectrum.f <- tmp$freq[1:spectrumLimit]
rn <- as.numeric(rownames(df))
windowPeaks <- which(tmp$freq>=4 & tmp$freq<=10)
windowValleys <- setdiff(which(tmp$freq>=2 & tmp$freq<=20),windowPeaks)
wide <- length(which(tmp$freq<=0.5))
nearbyWindow <- -wide:wide
valleyU <- c(0,0)
thetaResult <- data.frame(start=vector(mode='double',length=M),stop=vector(mode='double',length=M),valid=vector(mode='double',length=M))
df=log10(data.frame(m=tmp$spec[1:spectrumLimit],u=tmp$mtm$jk$upperCI[1:spectrumLimit],l=tmp$mtm$jk$lowerCI[1:spectrumLimit]))
N <- spectrumLimit
spectrum.m <- matrix(nrow=M,ncol=N)
spectrum.l <- matrix(nrow=M,ncol=N)
spectrum.u <- matrix(nrow=M,ncol=N)
start_time = Sys.time()
S <- foreach(j = 1:5, .export = c("mef_data","windows"),.packages=c('multitaper'))  %dopar% {
    spec.mtm( ts( mef_data[windows[,j]], frequency=32000 ), nw=4.0, k=7, jackknife=TRUE, plot=FALSE )
}
end_time = Sys.time()
print( paste0( "Runtime: ", end_time - start_time ) )

j = 1
      tmp <- S[[j]]
      print( length(tmp$spec) )
      df=log10(data.frame(m=tmp$spec[1:spectrumLimit],u=tmp$mtm$jk$upperCI[1:spectrumLimit],l=tmp$mtm$jk$lowerCI[1:spectrumLimit]))
      spectrum.m[j,] <- df$m
      spectrum.l[j,] <- df$l
      spectrum.u[j,] <- df$u
#      p <- ggplot(df,aes(x=seq(1,spectrumLimit))) + geom_line(aes(y=m),color='black') + geom_line(aes(y=u),color='green') + geom_line(aes(y=l),color='green')
      vars <- colnames(df)
      sdf <- as.data.frame( lapply(vars,loess.filter,span=0.15,rn=rn), col.names=vars )
#      p <- ggplot(df,aes(x=seq(1,spectrumLimit))) + geom_line(aes(y=m),color='black') + geom_line(aes(y=u),color='green') + geom_line(aes(y=l),color='green')
#      print(p)

      dsdx <- diff(sign(diff(sdf$l)))
      peaks <- which( dsdx[windowPeaks] < 0 ) + min(windowPeaks)
      peakLidx <- which( sdf$l == max(sdf$l[peaks]) )
      peakL <- max( df$l[peakLidx + nearbyWindow] )


      dsdx <- diff(sign(diff(sdf$u)))
      valleysIdx <- which( dsdx[windowValleys] > 0 )
      valleys <- windowValleys[valleysIdx]
      tmp1 <- which(valleys < peakLidx)
      tmp2 <- which(valleys > peakLidx)
      if ( length(tmp1)>0 & length(tmp2)>0 ) {
        idx <- valleys[length(tmp1)]
        valleyU[1] <- min( df$u[idx + nearbyWindow])
        idx <- valleys[tmp2[1]]
        valleyU[2] <- min( df$u[idx + nearbyWindow])
        result <- peakL>valleyU[1] & peakL>valleyU[2]
      } else {
        result <- FALSE
      }
      thetaResult$start[j] <- min(windows[,j])
      thetaResult$stop[j] <- max(windows[,j])
      thetaResult$valid[j] <- result

    steps <- paste0( 'step_', stepSeconds, '_width_', widthSeconds )
    # Spectrogram
    out_filename <- paste0( here(), '/', output_dir, '/', mef_fileroot, '_', iterationCount, '_', steps, '_', iterationCount, '.pdf' )
    print( out_filename )
    pdf( out_filename )
    image(spectrum.m[,1:spectrumLimit], col=grey(seq(0,1,length=256)))
    dev.off()
    # Spectrum
    out_filename <- paste0( here(), '/', output_dir, '/', mef_fileroot, '_', iterationCount, '_', steps, '_spectrum.RDdata' )
    save( spectrum.f, spectrum.m, spectrum.l, spectrum.u, file=out_filename )
    # Theta result
    out_filename <- paste0( here(), '/', output_dir, '/', mef_fileroot, '_', iterationCount, '_', steps, '_', iterationCount, '_theta.RDdata' )
    save( thetaResult, file=out_filename )
}

