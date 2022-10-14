dataDir <- "/Users/markrbower/Dropbox/Documents/Data/rat_test"
setwd( dataDir )
filename <- file.path( dataDir, "Ipsi_Hippoc_TT9/CSC34.ncs", fsep = .Platform$file.sep )
file.exists(filename)
meftools::ncs2mef( c(filename) )
  