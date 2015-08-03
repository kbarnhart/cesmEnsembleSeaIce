# Bootstrap code
library(reshape2)
library(boot)
library(ncdf)

ns=c("nh","sh")
rcp=c("RCP85", "RCP45")

# start profiler to look for efficiency improvements
Rprof("profileBootstrapCI.out")

for (n in 1:2)
{
  for (r in 1:2)
  {
      fileName<-paste("./justNSIF_ensembleAndBG.", ns[n],".", rcp[r], ".nc",collapse ='',sep='')
      print(fileName)

      # Read in data
      setwd("/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/")
      SIfile <- open.ncdf(fileName) # ncdf
      
      # put this in a loop for all 4 similar NCDFs 
      # copy the file and add the extra variable names. 
      
      
      print(paste("File",SIfile$filename,"contains",SIfile$nvars,"variables"))
      
      
      data<- get.var.ncdf(SIfile,varid='nSIF_ensemble')
      data_yr<- get.var.ncdf(SIfile,varid='time_ensemble')
      
      
      
      ni<-dim(data)[1]
      nj<-dim(data)[2]
      ny<-dim(data)[3]
      
      
      # Sample mean function for bootstrapping
      samplemean <- function(x, d) {
        return(mean(x[d]))
      }
      
      # Sample standard deviation function for bootstrapping
      samplesd <- function(x, d) {
        return(sd(x[d]))
      }
      
      
      # create ouptu matrices
      mean_ub <- array(data = NA, dim=c(ni, nj, ny))
      mean_lb <- array(data = NA, dim=c(ni, nj, ny))
      sd_ub <- array(data = NA, dim=c(ni, nj, ny))
      sd_lb <- array(data = NA, dim=c(ni, nj, ny))
      
      sd_var <- array(data = NA, dim=c(ni, nj, ny))
      mean_var <- array(data = NA, dim=c(ni, nj, ny))
      
      for (nii in 1:ni)
      {
        print(nii)
        for (njj in 1:nj)
        {
        
          for (nyy in 1:ny)
          {
            vals=data[nii,njj,nyy,]
            naTest<-is.na(mean(vals))
            if (naTest==FALSE)
            {
              waterTest<-(mean(vals)<360)
              
              if (waterTest==TRUE)
              {
            
              boot_result_mean <- boot(vals,statistic=samplemean,R=1000)
              mean_var[nii, njj, nyy] <- mean(vals,na.rm=T)
              mean_ci <- boot.ci(boot_result_mean,conf = 0.95,type = "basic")
              mean_ub[nii, njj, nyy] <- mean_ci$basic[5]
              mean_lb[nii, njj, nyy] <- mean_ci$basic[4]
              
              boot_result_stdev <- boot(vals,statistic=samplesd,R=1000)
              sd_var[nii, njj, nyy] <- sd(vals,na.rm=T)
              sd_ci <- boot.ci(boot_result_stdev,conf = 0.95, type = "basic")
              sd_ub[nii, njj, nyy] <- sd_ci$basic[5]
              sd_lb[nii, njj, nyy] <- sd_ci$basic[4] 

              }
            }
          }
        }
      }
      
      
      
      
      niDIM<-dim.def.ncdf( SIfile$dim$ni$name, SIfile$dim$ni$units, SIfile$dim$ni$vals, unlim=FALSE, create_dimvar=FALSE)
      njDIM<-dim.def.ncdf( SIfile$dim$nj$name, SIfile$dim$nj$units, SIfile$dim$nj$vals, unlim=FALSE, create_dimvar=FALSE)
      nyDIM<-dim.def.ncdf( SIfile$dim$time_ensemble$name, SIfile$dim$time_ensemble$units, SIfile$dim$time_ensemble$vals, unlim=FALSE, create_dimvar=TRUE)
      nverticesDIM<-dim.def.ncdf( SIfile$dim$nvertices$name, SIfile$dim$nvertices$units, SIfile$dim$nvertices$vals, create_dimvar=FALSE)
      
      newDims=list(niDIM, njDIM, nyDIM)
      
      # Make new output variables
      var_meanUB <-var.def.ncdf('meanUB', 'days', dim=newDims, 1e+30, longname='Upper bound for the mean number of sea ice free days', prec="single")
      var_meanLB <-var.def.ncdf('meanLB', 'days', dim=newDims, 1e+30, longname='Lower bound for the mean number of sea ice free days', prec="single")
      var_mean <-var.def.ncdf('mean', 'days', dim=newDims, 1e+30, longname='Mean number of sea ice free days', prec="single")
      
      var_sdUB <-var.def.ncdf('sdUB', 'days', dim=newDims, 1e+30, longname='Upper bound for standard deviation of the number of sea ice free days', prec="single")
      var_sdLB <-var.def.ncdf('sdLB', 'days', dim=newDims, 1e+30, longname='Lower bound for standard deviation of the number of sea ice free days', prec="single")
      var_sd <-var.def.ncdf('sd', 'days', dim=newDims, 1e+30, longname='Standard deviation of the number of sea ice free days', prec="single")
      
      # get lat/lon variables
      
      TLAT<-SIfile$var$TLAT
      TLON<-SIfile$var$TLON
      latt_bounds<-SIfile$var$latt_bounds
      lont_bounds<-SIfile$var$lont_bounds
      
      TLAT_vals<-get.var.ncdf(SIfile,varid='TLAT')
      TLON_vals<-get.var.ncdf(SIfile,varid='TLON')
      latt_bounds_vals<-get.var.ncdf(SIfile,varid='latt_bounds')
      lont_bounds_vals<-get.var.ncdf(SIfile,varid='lont_bounds')
      
      newFN<-paste("./justNSIF_ensembleAndBG.", ns[n],".", rcp[r], ".BootCI.TEST.nc",collapse ='',sep='')
      # Create the test file
      nc <- create.ncdf( newFN, list(var_meanUB,var_meanLB,var_sdUB,var_sdLB, var_sd, var_mean, TLAT, TLON, latt_bounds, lont_bounds) )
      # Write some data to the file
      put.var.ncdf(nc, var_meanUB, mean_ub)
      put.var.ncdf(nc, var_meanLB, mean_lb)
      put.var.ncdf(nc, var_sdUB, sd_ub)
      put.var.ncdf(nc, var_sdLB, sd_lb)
      put.var.ncdf(nc, var_sd, sd_var)
      put.var.ncdf(nc, var_mean, mean_var)
      
      put.var.ncdf(nc, TLAT, TLAT_vals)
      put.var.ncdf(nc, TLON, TLON_vals)
      put.var.ncdf(nc, latt_bounds, latt_bounds_vals)
      put.var.ncdf(nc, lont_bounds, lont_bounds_vals)
      close.ncdf(nc)
      rm(nc)
      
      close.ncdf(SIfile)
      rm(SIfile)

  }
}

# profile code at end
Rprof(NULL)
summaryRprof("profileBootstrapCI.out")
