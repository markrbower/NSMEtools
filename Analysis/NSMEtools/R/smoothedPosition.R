smoothedPosition  <- function() {
# smoothedPosition.R
#
# Mark R. Bower
# October 22, 2022
# Yale
library( stringr )

fillFirstLast <- function(x) {
  idx <- which( !is.na(x) )
  first <- min(idx)
  for (j in 1:first) {x[j] <- x[first]}
  last <- max(idx)
  for (j in last:length(x)) {x[j] <- x[last]}
  return(x)
}

# Load data
rootDir <- '/Users/markrbower/Dropbox/Documents/Concepts/2022_07_01_NSMEtools/NSMEtools/Presentation/Figures/Stim5'
filename <- 'b_stim/video_position.txt'
out_filename <- str_replace( filename, 'position', 'position_smoothed')
dirs <- list.files(path=rootDir,pattern="^\\d+_\\d+$")
for ( dir in dirs ) {
  dataDir <- paste0( rootDir, '/',  dir )
  setwd( dataDir )
  fid <- file( filename, 'r' )
  t = as.numeric( strsplit(readLines(fid,1),",")[[1]] )
  x = as.numeric( strsplit(readLines(fid,1),",")[[1]] )
  y = as.numeric( strsplit(readLines(fid,1),",")[[1]] )
  close(fid)

  # Interpolate zeros
  df = data.frame( t=t, x=x, y=y )
  
  zx <- which( df$x == 0 )
  df$x[zx] <- NA
  df$x <- fillFirstLast( df$x )
  df <- df %>%  mutate( x = na.approx(x), maxgap=Inf, na.rm=FALSE )
  zy <- which( df$y == 0 )
  df$y[zy] <- NA
  df$y <- fillFirstLast( df$y )
  df <- df %>%  mutate( y = na.approx(y), maxgap=Inf, na.rm=FALSE )
  
  # Interpolate big differences
  dx <- diff(df$x)
  zx <- which( abs(dx) > 30 )
  df$x[zx+1] <- NA
  df$x <- fillFirstLast( df$x )
  df <- df %>%  mutate( x = na.approx(x), maxgap=Inf, na.rm=FALSE )
  dy <- diff(df$y)
  zy <- which( abs(dy) > 30 )
  df$y[zy+1] <- NA
  df$y <- fillFirstLast( df$y )
  df <- df %>%  mutate( y = na.approx(y), maxgap=Inf, na.rm=FALSE )

  # Save edited position file.
  sink( file=out_filename, type='output' )
  cat( paste0( df$t, collapse=',') )
  cat( '\n' )
  cat( paste0( df$x, collapse=',') )
  cat( '\n' )
  cat( paste0( df$y, collapse=',') )
  cat( '\n' )
  sink()
}
}

