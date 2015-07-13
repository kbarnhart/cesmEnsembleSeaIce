# Bootstrap code
library(reshape2)
library(boot)
library(ncdf)

# Read in data
setwd("/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/")
SIfile <- open.ncdf("./justNSIF_ensembleAndBG.nc") # ncdf


print(paste("File",SIfile$filename,"contains",SIfile$nvars,"variables"))
for( i in 1:SIfile$nvars ) {
  v <- SIfile$var[[i]]
  print(paste("Here is information on variable number",i))
  print(paste(" Name: ",v$name))
  print(paste(" Units:",v$units))
  print(paste(" Missing value:",v$missval))
  print(paste(" # dimensions :",v$ndims))
  print(paste(" Variable size:",v$varsize))
}

data<- get.var.ncdf(SIfile,varid='nSIF_ensemble')
data_yr<- get.var.ncdf(SIfile,varid='time_ensemble')

close.ncdf(SIfile)
rm(SIfile)


# Sample mean function for bootstrapping
samplemean <- function(x, d) {
  return(mean(x[d]))
}

# Sample standard deviation function for bootstrapping
samplesd <- function(x, d) {
  return(sd(x[d]))
}

ny<-80
njj<-60
nii<-156

vals=data[nii,njj,ny,]

boot_result <- boot(vals,statistic=samplemean,R=1000)
boot_result
mean(vals,na.rm=T)
boot.ci(boot_result,conf = 0.95)

# Stdev

std_df <- read.csv("./stdsVals.txt",header=F)
transp_std <- t(std_df)
samplesd <- function(x, d) {
  return(sd(x[d]))
}


boot_result_stdev <- boot(transp_mean[1:30,1],statistic=samplesd,R=1000)
boot_result_stdev
sd(transp_mean[,1],na.rm=T)

boot.ci(boot_result_stdev)
