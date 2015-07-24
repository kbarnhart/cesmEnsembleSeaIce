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


# create ouptu matrices
mean_ub <- array(data = NA, dim=c()

# eventually this will be a loop

ny<-80
njj<-60
nii<-156

vals=data[nii,njj,ny,]

boot_result_mean <- boot(vals,statistic=samplemean,R=1000)
mean <- mean(vals,na.rm=T)
mean_ci <- boot.ci(boot_result_mean,conf = 0.95,type = "basic")
mean_ub <- mean_ci$basic[4]
mean_lb <- mean_ci$basic[5]

boot_result_stdev <- boot(vals,statistic=samplesd,R=1000)
sd <- sd(vals,na.rm=T)
sd_ci <- boot.ci(boot_result_stdev,conf = 0.95, type = "basic")
sd_ub <- sd_ci$basic[4]
sd_lb <- sd_ci$basic[5]
