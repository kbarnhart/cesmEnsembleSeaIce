# Respones to meotreviewer questions
## Setup ####
library(ncdf4)
library(ggplot2)
library(stringr)
#setwd("~/Dropbox/seaIceDataForChris/")
setwd("/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/")
getwd()

outF<-"~/git/cesmEnsembleSeaIce/assumptionCheckFigures"

SIfile <- nc_open("./DLM__MCMC_results_final.nh.RCP85.nc")

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


head(SIfile$var)


# Variables necessary to get residuals
#nSIF
#[1] " # dimensions : 4"
#[1] " Variable size: 320" " Variable size: 104" " Variable size: 250" " Variable size: 30" 

#dlm_levelmean
#[1] " # dimensions : 3"
#[1] " Variable size: 320" " Variable size: 104" " Variable size: 250"

#So, to get residuals, take mean of nSIF across 30, compare against dlm_levelmean.

## Create residual dataframe ####

## Extract observations, apply mean across 30 samples ##
observations <- ncvar_get(SIfile,varid='nSIF')
table(is.nan(observations))
observations[is.nan(observations)] <- NA
observations_mean <- apply(X=observations,MARGIN=c(1,2,3),FUN=mean,na.rm=T)

## Compare against model results ## 
model_pred <- ncvar_get(SIfile,varid='dlm_levelmean') 
mod_pred_ub<- ncvar_get(SIfile,varid='dlm_upperBound')
mod_pred_lb<- ncvar_get(SIfile,varid='dlm_lowerBound')

shift<-ncvar_get(SIfile,varid='shift_year')
shift[is.nan(shift)] <- NA

emerge<-ncvar_get(SIfile,varid='emergence_year')
emerge[is.nan(emerge)] <- NA


mod_pred_broadcast=array(0,dim=dim(observations))
for (idx in 1:dim(observations)[4]){
mod_pred_broadcast[,,,idx]=model_pred
}
residuals <- model_pred - observations_mean
residuals_all <- mod_pred_broadcast - observations

# Sanity check existence of .
table(is.na(model_pred))
table(is.na(observations_mean))
table(is.na(residuals))

summary(model_pred)
summary(observations_mean)
summary(residuals)

## Look at individual points
# Drew Point: {'ni': 209, 'nj': 64}
model_pred[209,64,]
observations_mean[209,64,]
drewpoint_df <- residuals[209,64,]
drewpoint_noNA_df <- drewpoint_df[!is.na(drewpoint_df)]
qplot(drewpoint_df)
qqnorm(drewpoint_df)
qqline(drewpoint_df)

# Partial auto-correlation plots on residuals
acf(drewpoint_noNA_df)
pacf(drewpoint_noNA_df) 
# Very close to independent.  
# Minor correlation in the residuals, likely because of smoothing process

# Heteroskedacity 
qplot(1:length(drewpoint_noNA_df),
      drewpoint_noNA_df) + geom_smooth()
# Constant variance throughout entire range

# Effect of truncation at end
# Model is not calculated at end for drew point.  What about some other point?
qplot(1:length(observations_mean[209,64,]),observations_mean[209,64,])
qplot(1:length(model_pred[209,64,]),model_pred[209,64,])

# Variable transformation
# Sigmoidal transforms of relative freq does little.  
# Other suggested transforms (sqrt) does not make sense, 
# as this is a different type of count data for which 
drewpoint_logit <- log(observations_mean[209,64,]/365) - 
  log(1 - observations_mean[209,64,]/365)
qplot(1:length(drewpoint_logit),drewpoint_logit)


# Bering Strait: {'ni': 200, 'nj': 44} Entirely open water at end
model_pred[200,44,]
observations_mean[200,44,]
beringstrait_df <- residuals[200,44,]
beringstrait_noNA_df <- beringstrait_df[!is.na(beringstrait_df)]
qplot(beringstrait_df)
qqnorm(beringstrait_df)
qqline(beringstrait_df)

# Partial auto-correlation plots on residuals
acf(beringstrait_noNA_df)
pacf(beringstrait_noNA_df) 
# Very close to independent.  
# Minor correlation in the residuals, likely because of smoothing process


# Heteroskedacity 
qplot(1:length(beringstrait_noNA_df),
      beringstrait_noNA_df) + geom_smooth()
# Constant variance throughout entire range

# Effect of truncation at end
# Model is not calculated at end for drew point.  What about some other point?
qplot(1:length(observations_mean[200,44,]),observations_mean[200,44,])
qplot(1:length(model_pred[200,44,]),model_pred[200,44,])

# Variable transformation
# Bering Strait transforms 
beringstrait_logit <- log(observations_mean[200,44,]/365) - 
  log(1 - observations_mean[200,44,]/365)
qplot(1:length(beringstrait_logit),beringstrait_logit)
qplot(1:length(beringstrait_logit),tanh(beringstrait_logit))



# Central A - CAA: {'ni': 159, 'nj': 93} Entirely open water at end
model_pred[159,93,]
observations_mean[159,93,]
CAA_df <- residuals[159,93,]
CAA_noNA_df <- CAA_df[!is.na(CAA_df)]

pdf("CAA_QQ_overall_then_separated.pdf")
#qplot(CAA_df)
qqnorm(CAA_df)
qqline(CAA_df)

# CAA is decidedly not normal overall.  
# However, if one breaks the group into two pieces, .
# Piece one: Truncated at zero
#qplot(CAA_noNA_df[1:100])
qqnorm(CAA_noNA_df[1:100])
qqline(CAA_noNA_df[1:100])

# Piece two: Non-truncated
#qplot(CAA_noNA_df[101:length(CAA_noNA_df)])
qqnorm(CAA_noNA_df[101:length(CAA_noNA_df)])
qqline(CAA_noNA_df[101:length(CAA_noNA_df)])
dev.off()


# Partial auto-correlation plots on residuals
acf(CAA_noNA_df)
pacf(CAA_noNA_df) 
# Very close to independent.  
# Minor correlation in the residuals, likely because of smoothing process


# Heteroskedacity 
qplot(1:length(CAA_noNA_df),
      CAA_noNA_df) + geom_smooth()
# Constant variance throughout entire range

# Effect of truncation at end
# Model is not calculated at end for drew point.  What about some other point?
qplot(1:length(observations_mean[159,93,]),observations_mean[159,93,])
qplot(1:length(model_pred[159,93,]),model_pred[159,93,])

# Variable transformation
# Bering Strait transforms 
CAA_logit <- log(observations_mean[159,93,]/365) - 
  log(1 - observations_mean[159,93,]/365)
qplot(1:length(CAA_logit),CAA_logit)


## Loop through points and checks ####
startYear<-1850
require(reshape)
fnt_sz<-24
## Read in list of points:
#pointlist <- read.csv("~/git/cesmEnsembleSeaIce/PointList.csv")
pointlist <- read.csv("~/git/cesmEnsembleSeaIce/assumptionCheckPoints.csv")

for (idx in 1:nrow(pointlist)){
  ## dataframe setup for each site
  print(pointlist$sitename[idx])
  
  site_all_obs_df <-observations[pointlist$ni[idx],pointlist$nj[idx],,]
  site_obs_df <- observations_mean[pointlist$ni[idx],pointlist$nj[idx],]
  site_model_df <- model_pred[pointlist$ni[idx],pointlist$nj[idx],]
  
  site_model_ub_df <- mod_pred_ub[pointlist$ni[idx],pointlist$nj[idx],]
  site_model_lb_df <- mod_pred_lb[pointlist$ni[idx],pointlist$nj[idx],]
  
  site_resid_df <- residuals[pointlist$ni[idx],pointlist$nj[idx],]
  site_resid_df_all <-residuals_all[pointlist$ni[idx],pointlist$nj[idx],,]
  
  site_resid_noNA_df <- site_resid_df[!is.na(site_resid_df)]
  timeAxis<-startYear:(startYear+length(site_model_df)-1)

  sy<-shift[pointlist$ni[idx],pointlist$nj[idx]]
  ey<-emerge[pointlist$ni[idx],pointlist$nj[idx]]

  
  ## If there is a truncated vs non-truncated section, add limits here ##
  # If there is a point < 3 and a point >=3, set idx as first occurence of this in bound 1
  # If there is a point <= 362 and a point >362, set idx as first occurence of this in bound 2
  # Result is NA if this does not exist
  pointlist$bound1[idx] <- which(any(site_model_df<3) & site_model_df >= 3)[1]
  pointlist$bound2[idx] <- which(any(site_model_df<362) & site_model_df >= 362)[1] 
  
  
  # Normality checks. Note the impact of truncated distributions.
  # Also, no automated way to separate truncated section from linearly changing sections
  filename <- str_c(outF, '/',pointlist$sitename[idx],"_QQ.pdf")
  pdf(filename)
  qqnorm(site_resid_noNA_df, main='',cex=1.725)
  print(qqline(site_resid_noNA_df))
  dev.off()
  
  qq_df<-data.frame(site_resid_noNA_df)
  
  y <- quantile(qq_df$site_resid_noNA_df, c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  
  filename <- str_c(outF, '/',pointlist$sitename[idx],"_QQ_v2.pdf")
  pdf(filename)
  qplot_out<-ggplot(qq_df, aes(sample=site_resid_noNA_df)) + geom_point(stat = "qq",size=4, pch=1)+
    geom_abline(slope = slope, intercept = int, colour='blue', linetype=3)+
    xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + 
    theme(panel.background = element_blank())+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
    theme(axis.text = element_text(size = fnt_sz),
          axis.title = element_text(size=fnt_sz))
  
  print(qplot_out)
  dev.off()

  # If bound1 and bound 2 exist
  # 4 cases: 
  # Neither bound exists
  if(is.na(pointlist$bound1[idx]) & is.na(pointlist$bound2[idx]) ){}
  # Bound 1 exists, bound 2 does not
  else if(!is.na(pointlist$bound1[idx])  & is.na(pointlist$bound2[idx]) ){
    qq_dat_l<-site_resid_noNA_df[1:(pointlist$bound1[idx]-70)]
    qq_l_df<-data.frame(qq_dat_l)
    filename <- str_c(outF, '/',pointlist$sitename[idx],"_QQlowertruncated.pdf")
    pdf(filename)
    y <- quantile(qq_l_df$qq_dat_l, c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    qplot_out<-ggplot(qq_l_df, aes(sample=qq_dat_l)) + geom_point(stat = "qq",size=4, pch=1)+
      geom_abline(slope = slope, intercept = int, colour='blue', linetype=3)+
      xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + 
      theme(panel.background = element_blank())+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
      theme(axis.text = element_text(size = fnt_sz),
            axis.title = element_text(size=fnt_sz))
    print(qplot_out)
    dev.off()
    
    qq_dat_u<-site_resid_noNA_df[(pointlist$bound1[idx]-70 + 1):length(site_resid_noNA_df)]
    qq_u_df<-data.frame(qq_dat_u)
    
    filename <- str_c(outF, '/',pointlist$sitename[idx],"_QQuntruncated.pdf")
    pdf(filename)
    y <- quantile(qq_u_df$qq_dat_u, c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    qplot_out<-ggplot(qq_u_df, aes(sample=qq_dat_u)) + geom_point(stat = "qq",size=4, pch=1)+
      geom_abline(slope = slope, intercept = int, colour='blue', linetype=3)+
      xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + 
      theme(panel.background = element_blank())+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
      theme(axis.text = element_text(size = fnt_sz),
            axis.title = element_text(size=fnt_sz))
    print(qplot_out)
    dev.off()
  }
  # Bound 2 exists, bound 1 does not
  else if (is.na(pointlist$bound1[idx])  & !is.na(pointlist$bound2[idx]) ){
    qq_dat_un<-site_resid_noNA_df[1:(pointlist$bound2[idx]-70 + 1)]
    qq_u_df<-data.frame(qq_dat_un)
    filename <- str_c(outF, '/',pointlist$sitename[idx],"_QQuntruncated.pdf")
    pdf(filename)
    y <- quantile(qq_u_df$qq_dat_un, c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    qplot_out<-ggplot(qq_u_df, aes(sample=qq_dat_un)) + geom_point(stat = "qq",size=4, pch=1)+
      geom_abline(slope = slope, intercept = int, colour='blue', linetype=3)+
      xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + 
      theme(panel.background = element_blank())+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
      theme(axis.text = element_text(size = fnt_sz),
            axis.title = element_text(size=fnt_sz))
    print(qplot_out)
    dev.off()
    
    qq_dat_l<-site_resid_noNA_df[(pointlist$bound2[idx]-70 + 1):length(site_resid_noNA_df)]
    qq_up_df<-data.frame(qq_dat_l)
    
    filename <- str_c(outF, '/',pointlist$sitename[idx],"_QQuppertruncated.pdf")
    pdf(filename)
    y <- quantile(qq_up_df$qq_dat_l, c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    qplot_out<-ggplot(qq_up_df, aes(sample=qq_dat_l)) + geom_point(stat = "qq",size=4, pch=1)+
      geom_abline(slope = slope, intercept = int, colour='blue', linetype=3)+
      xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + 
      theme(panel.background = element_blank())+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
      theme(axis.text = element_text(size = fnt_sz),
            axis.title = element_text(size=fnt_sz))
    print(qplot_out)
    dev.off()
  }
  # Both bound 1 and bound 2 exist
  else if (!is.na(pointlist$bound1[idx]) & !is.na(pointlist$bound2[idx]) ){
    qq_dat_l<-site_resid_noNA_df[1:(pointlist$bound1[idx]-70)]
    qq_l_df<-data.frame(qq_dat_l)
    filename <- str_c(outF, '/',pointlist$sitename[idx],"_QQlowertruncated.pdf")
    pdf(filename)
    y <- quantile(qq_l_df$qq_dat_l, c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    qplot_out<-ggplot(qq_l_df, aes(sample=qq_dat_l)) + geom_point(stat = "qq",size=4, pch=1)+
      geom_abline(slope = slope, intercept = int, colour='blue', linetype=3)+
      xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + 
      theme(panel.background = element_blank())+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
      theme(axis.text = element_text(size = fnt_sz),
            axis.title = element_text(size=fnt_sz))
    print(qplot_out)
    dev.off()
    
    qq_dat_un<-site_resid_noNA_df[(pointlist$bound1[idx]-70 + 1):(pointlist$bound2[idx]-70)]
    qq_un_df<-data.frame(qq_dat_u)
    
    filename <- str_c(outF, '/',pointlist$sitename[idx],"_QQuntruncated.pdf")
    pdf(filename)
    y <- quantile(qq_un_df$qq_dat_u, c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    qplot_out<-ggplot(qq_un_df, aes(sample=qq_dat_un)) + geom_point(stat = "qq",size=4, pch=1)+
      geom_abline(slope = slope, intercept = int, colour='blue', linetype=3)+
      xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + 
      theme(panel.background = element_blank())+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
      theme(axis.text = element_text(size = fnt_sz),
            axis.title = element_text(size=fnt_sz))
    print(qplot_out)
    dev.off()
    
    qq_dat_up<-site_resid_noNA_df[(pointlist$bound2[idx]-70):length(site_resid_noNA_df)]
    qq_up_df<-data.frame(qq_dat_up)
    
    filename <- str_c(outF, '/',pointlist$sitename[idx],"_QQuppertruncated.pdf")
    pdf(filename)
    y <- quantile(qq_up_df$qq_dat_up, c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    qplot_out<-ggplot(qq_up_df, aes(sample=qq_dat_up)) + geom_point(stat = "qq",size=4, pch=1)+
      geom_abline(slope = slope, intercept = int, colour='blue', linetype=3)+
      xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + 
      theme(panel.background = element_blank())+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
      theme(axis.text = element_text(size = fnt_sz),
            axis.title = element_text(size=fnt_sz))
    print(qplot_out)
    dev.off()
    
  }

  
  # plot the model results (include upper and lower bounds, and truncation locations)
  filename <- str_c(outF, '/', pointlist$sitename[idx],"_DLMresults.pdf")
  pdf(filename)
  dlm_df <- data.frame(tA = timeAxis,SM = site_model_df, SM_ub=site_model_ub_df, SM_lb=site_model_lb_df) # Separate df for mean of simulation results
  sim_results_mean_dlm <- ggplot(dlm_df, aes(timeAxis,site_model_df))+geom_line()+
        xlim(c(1920,2100))+ ylim(c(-1,366))+
        xlab('Year')+ylab('Number of open water days')+
        geom_line(aes(x=timeAxis, y = site_model_ub_df, colour = "FF9999"))+
        geom_line(aes(x=timeAxis, y = site_model_lb_df, colour = "FF9999"))+
        geom_vline(xintercept=pointlist$bound1[idx]+startYear, linetype=2)+
        geom_vline(xintercept=pointlist$bound2[idx]+startYear)+
        geom_vline(xintercept=sy, linetype=3) +
        geom_vline(xintercept=ey, linetype=4) +
        theme(legend.position="none")+
        theme(panel.background = element_blank())+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
        theme(axis.text = element_text(size = fnt_sz),
          axis.title = element_text(size=fnt_sz))
  print(sim_results_mean_dlm)
  dev.off()
  
  obs_df_m <- melt(site_all_obs_df, id.vars=c("mod"))
  obs_df_m$timeAxis <- timeAxis
  
  # plot the original observations (includ truncation locations)
   filename <- str_c(outF, '/', pointlist$sitename[idx],"_origObs.pdf")
   pdf(filename)
   mean_obs_df <- data.frame(tA = timeAxis,SO = site_obs_df) # Separate df for mean of simulation results
  
   sim_results_mean <- ggplot(obs_df_m,aes(timeAxis,value)) + geom_line(aes(alpha=0.5)) + geom_point(aes(alpha=0.55)) + #aes(alpha=0.5, colour="black")) + 
    xlim(c(1920,2100)) + ylim(c(0,365)) + 
    xlab("Year") + ylab("Number of open water days") + 
    geom_line(data=mean_obs_df,aes(x=tA, y=SO, alpha=1)) + 
    geom_vline(xintercept=pointlist$bound1[idx]+startYear, linetype=2) + 
    geom_vline(xintercept=pointlist$bound2[idx]+startYear) +
    geom_vline(xintercept=sy, linetype=3) +
    geom_vline(xintercept=ey, linetype=4) +
    theme(legend.position="none")+
    theme(panel.background = element_blank())+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
    theme(axis.text = element_text(size = fnt_sz),
          axis.title = element_text(size=fnt_sz))
          
  print(sim_results_mean)
  dev.off()
  
  # DLM results and model output
  filename <- str_c(outF, '/', pointlist$sitename[idx],"_origObs+DLM.pdf")
  pdf(filename)#, width=4.7244, height=4.7244)
  sim_results_mean <- ggplot(obs_df_m,aes(timeAxis,value)) + geom_line(aes(alpha=0.5)) + geom_point(aes(alpha=0.55)) + #aes(alpha=0.5, colour="black")) + 
    xlim(c(1920,2100)) + ylim(c(0,365)) + 
    xlab("Year") + ylab("Number of open water days") + 
    geom_line(data=mean_obs_df,aes(x=tA, y=SO, alpha=1)) + 
    geom_line(data=dlm_df ,aes(x=tA, y=SM, alpha=1)) + 
    geom_line(data=dlm_df,aes(x=tA, y=SM_ub, alpha=1,colour = "FF9999")) +
    geom_line(data=dlm_df,aes(x=tA, y=SM_lb, alpha=1,colour = "FF9999")) +
    geom_vline(xintercept=pointlist$bound1[idx]+startYear, linetype=2) +
    geom_vline(xintercept=pointlist$bound2[idx]+startYear) +
    geom_vline(xintercept=sy, linetype=3) +
    geom_vline(xintercept=ey, linetype=4) +
    theme(legend.position="none")+
    theme(panel.background = element_blank())+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
    theme(axis.text = element_text(size = fnt_sz),
          axis.title = element_text(size=fnt_sz))
  print(sim_results_mean)
  dev.off()
  
  # Partial Auto-correlation checks
  filename <- str_c(outF, '/',pointlist$sitename[idx],"_PACF.pdf")
  pdf(filename)
  print(pacf(site_resid_noNA_df, main='', cex=1.725))
  dev.off()
  
  filename <- str_c(outF, '/',pointlist$sitename[idx],"_PACF_v2.pdf")
  bacf <- pacf(site_resid_noNA_df, main='', plot = FALSE)
  CI<- qnorm(c(0.025, 0.975))/sqrt(length(site_resid_noNA_df))
  bacfdf <- with(bacf, data.frame(lag, acf))
  pdf(filename)
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    ylim(c(-0.7,0.7))+
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0))+
    geom_hline(aes(yintercept = CI[1], linestyle=3, colour='blue')) +
    geom_hline(aes(yintercept = CI[2], linestyle=3, colour='blue')) +
    xlab("Lag") + ylab("Partial ACF") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
    theme(panel.background = element_blank())+
    theme(axis.text = element_text(size = fnt_sz),
          axis.title = element_text(size=fnt_sz))
  print(q)
  dev.off()
  
  # Residuals of mean check
  
  
  resid_df<-data.frame(tA=timeAxis[72:250], resid=site_resid_noNA_df)
  
  filename <- str_c(outF, '/',pointlist$sitename[idx],"_Residuals_mean.pdf")
  pdf(filename)
  print(ggplot(data=resid_df, aes(x=tA, y=resid)) + geom_point() + 
              xlim(c(1920,2100)) +  xlab('Year') + ylab('Residual (days)') +
              geom_smooth()+theme(panel.background = element_blank() )+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
          theme(axis.text = element_text(size = fnt_sz),
          axis.title = element_text(size=fnt_sz)))
  dev.off()

# residuals of all points  
 resid_df_all<-data.frame(tA=timeAxis, resid=site_resid_df_all)
 resid_df_m <- melt(resid_df_all, id.vars=c("tA"))
 filename <- str_c(outF, '/',pointlist$sitename[idx],"_Residuals_all.pdf")
 pdf(filename)
 print(ggplot(data=resid_df_m, aes(x=tA, y=value)) + geom_point(alpha=0.5) + 
        xlim(c(1920,2100)) +  xlab('Year') + ylab('Residual (days)') +
        geom_smooth()+theme(panel.background = element_blank())+
         theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
         theme(axis.text = element_text(size = fnt_sz),
               axis.title = element_text(size=fnt_sz))) 
 dev.off()

# Logit variable transformation. Not useful
  site_logit <- log(site_obs_df/365) - 
    log(1 - site_obs_df/365)
  filename <- str_c(outF, '/',pointlist$sitename[idx],"_LogitTransform.pdf")
  pdf(filename)
  print(qplot(timeAxis,site_logit,xlim=c(1920,2100), ylim=c(-2,7),
              xlab='Year',ylab='Logit of the number of open water days')+
          theme(axis.text = element_text(size = fnt_sz),
                axis.title = element_text(size=fnt_sz)))
  dev.off()
  
  rm(resid_df)
  rm(dlm_df)
  rm(mean_obs_df)
}

pointlist
 

