print 'importing modules'
import nio
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from datetime import timedelta
import os
import glob

## to add:
#1. compatability with sea ice obs times (2 day obs, real years)
#2. distance to nearest ice edge analysis
#3. freeze up after Jan 1 compatability
#3. make first day of open water 1, and last 365 as default
#. make nsif 365 as default


import cProfile, pstats, StringIO
pr = cProfile.Profile()
pr.enable()

owThresh=15 # sea ice concentration for "open Water"
midSIF=250  # day of year for "middle of SIF season" 215~= August 1st
midIce=70

fn='/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/satelliteObs.timeseries.nc'
daysPerMonth=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
cumDaysInYear=np.cumsum(daysPerMonth)

daysPerMonthLeap=[31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
cumDaysInYearLeap=np.cumsum(daysPerMonthLeap)

startTime = datetime.now()
print 'Starting averaging of', fn

fn_monthOut=fn[:-2]+'monthlyAvg.nc'
#
f=nio.open_file(fn, 'r')
numyears=35
t=f.variables['time'][:]

reftime=datetime(1979, 1, 1, 0, 0, 0)

ni=f.dimensions['ni']
nj=f.dimensions['nj'] 

# for each cell do the following:
# 1. Create monthly averages of SIC

#** add montjly sea ice area estimate
nummonths=numyears*12

# do the montly averaging
monthlyTimestep=np.zeros(numyears*12)
monthlyTimeBounds=np.zeros((numyears*12,2))

monthlyAvgSIC=np.zeros((numyears*12, nj, ni))
arrayItter=0
satKey='satelliteSIC'
startYear=1979
for i in range(numyears): 
    #print 'year = ',i+int(t[0]/365)
    for j in range(12):
        startDate=datetime(startYear+i, j+1, 1, 0, 0, 0)
        try:
            endDate=datetime(startYear+i, j+2, 1, 0, 0, 0)
        except: # this will break for december
            endDate=datetime(startYear+i+1, 1, 1, 0, 0, 0)
        startOrd=startDate-reftime
        endOrd=endDate-reftime
        
                
        inRange=(t>=startOrd.days)*(t<endOrd.days)
        indStart=np.where(inRange==True)[0][0]
        indStop=np.where(inRange==True)[0][-1]+1
        
        
        #indStart=int(i*daysPerYear+np.remainder(cumDaysInYear[j-1], daysPerYear))
        #indStop=int(i*daysPerYear+cumDaysInYear[j]-1)
        #
        print i, j, indStart, indStop, (datetime.now()-startTime)
        selt=f.variables['time'][indStart:indStop]
        selIce=f.variables[satKey][indStart:indStop,:,:]
        meanIce=selIce.mean(axis=0)
        meanIce[selIce[0,:,:]==f.variables[satKey].__dict__['_FillValue']]=f.variables[satKey].__dict__['_FillValue']
        monthlyAvgSIC[arrayItter,:,:]=meanIce
        monthlyTimestep[arrayItter]=round(selt.mean(axis=0))
        monthlyTimeBounds[arrayItter, 0]=f.variables['time_bounds'][indStart,0]
        monthlyTimeBounds[arrayItter, 1]=f.variables['time_bounds'][indStop,1]
        arrayItter+=1
        
## Create this monthly averaged file as a new netcdf
fn_monthOut=fn[:-2]+'monthlyAvg.nc'
fMonth=Dataset(fn_monthOut, 'w',format='NETCDF4')

# create all the dimentions, set time to unlimited 
for k in f.dimensions.keys():
    if f.unlimited(k)==True:
        fMonth.createDimension(k, None)
    else:
        fMonth.createDimension(k, f.dimensions[k])
        
# use the netCDF4 instead of pyNIO since it seems to work much better with unlimited variables      
fMonthVars={}
for key in {'TLAT', 'TLON','time_bounds', 'time'}:
    #print 'creating ', key
    # the netCDF4 module requires that if a fill value exists, it must be declared when the variable is created.
    try:
        fMonthVars[key]=fMonth.createVariable(key, f.variables[key].typecode(), f.variables[key].dimensions, fill_value=f.variables[key].__dict__['_FillValue'])
    except:
        fMonthVars[key]=fMonth.createVariable(key, f.variables[key].typecode(), f.variables[key].dimensions)
    # sett all the attribute keys.
    atts = f.variables[key].__dict__
    for attKey in atts.keys():
        if attKey != '_FillValue':
            setattr(fMonth.variables[key],attKey,atts[attKey])         
            
# create the montly averaged sea ice variable 
monthAvgKey=satKey+'_monthAvg'
fMonthVars[monthAvgKey]=fMonth.createVariable(monthAvgKey, f.variables[satKey].typecode(), f.variables[satKey].dimensions,fill_value=f.variables[satKey].__dict__['_FillValue'])
#print 'creating aice_d_monthAvg'
atts = f.variables[satKey].__dict__
for attKey in atts.keys():
    if attKey is not '_FillValue':
        setattr(fMonth.variables[monthAvgKey],attKey,atts[attKey])
setattr(fMonth.variables[monthAvgKey],'long_name','mean monthly ice area  (aggregate)')   

# put data into variables, first the ones we are copying over. 
#print 'putting data into standard variables'
for key in {'TLAT', 'TLON'}:
    fMonthVars[key][:,:]=f.variables[key][:]
    
# now, the ones we have created (those with time as a dimention)
fMonthVars['time'][:]=monthlyTimestep
fMonthVars['time_bounds'][:,:]=monthlyTimeBounds
fMonthVars[satKey+'_monthAvg'][:,:,:]=monthlyAvgSIC

# close and delete the output netCDF variable, retain f, so we can do the next part of the analysis.
fMonth.close() 
del fMonth

print'monthly averaging finished, moving onto analysis', (datetime.now()-startTime)
# 2. First day, last day, and duration of open water
# Save these things to a new netcdf 

fMonth=Dataset(fn_monthOut, 'r',format='NETCDF4')
monthlyAvgSIC=fMonth.variables[satKey+'_monthAvg'][:,:,:]
fMonth.close()

iceMax=monthlyAvgSIC.max(axis=0)
fillVal=f.variables[satKey].__dict__['_FillValue']
landMask=f.variables[satKey][0,:,:]==fillVal

numAnalyisYears=numyears#-1 

iceVal=1000.
firstDOY=iceVal*np.ones((numAnalyisYears, nj, ni))
lastDOY=iceVal*np.ones((numAnalyisYears, nj, ni))
firstDOYcontinous=iceVal*np.ones((numAnalyisYears, nj, ni))
lastDOYcontinous=iceVal*np.ones((numAnalyisYears, nj, ni))
duration=iceVal*np.ones((numAnalyisYears, nj, ni))
durationContinous=iceVal*np.ones((numAnalyisYears, nj, ni))
numSIF=365.*np.ones((numAnalyisYears, nj, ni))
springMix=iceVal*np.ones((numAnalyisYears, nj, ni))
fallMix=iceVal*np.ones((numAnalyisYears, nj, ni))

yearlyTimestep=np.zeros(numyears)
yearlyTimeBounds=np.zeros((numyears,2))

iceAffected=(iceMax>0.1)*(iceMax<101)
land=iceMax>100.5
#totalAffected=iceAffected.sum(axis=None)
notIceAffected=iceAffected==False

for i in range(numyears):
    startDate=datetime(startYear+i, 1, 1, 0, 0, 0)+timedelta(days=midIce)
    endDate=datetime(startYear+i+1, 1, 1, 0, 0, 0)+timedelta(days=midIce)
    
    numdays=endDate-startDate
    
    startOrd=startDate-reftime
    endOrd=endDate-reftime
 
    inRange=(t>=startOrd.days)*(t<endOrd.days)
    indStart=np.where(inRange==True)[0][0]
    indStop=np.where(inRange==True)[0][-1]+1
        
    selt=f.variables['time'][indStart:indStop]
    yearlyTimestep[i]=round(selt.mean(axis=0))
    yearlyTimeBounds[i, 0]=f.variables['time_bounds'][indStart,0]
    yearlyTimeBounds[i, 1]=f.variables['time_bounds'][indStop,1]
    
    yearConc=f.variables[satKey][indStart:indStop,:,:]
    
    percDone=round(i*100.0/numyears)
    print 'year = ',i,'  : ', str(indStart), str(indStop) ,str(percDone), '% of this model ', (datetime.now()-startTime)

    for nii in range(ni):
        for nji in range(nj):         
            if iceAffected[nji, nii]==True: # not land    
                selConc=yearConc[:,nji, nii] 
    
       	#we only want this analysis if there are fully sea ice covered days
                fullSIBool=any((selConc>80)*(selConc<101))
                if fullSIBool==True:
                    numFullSI=sum((selConc>80)*(selConc<101))
                else: 
                    numFullSI=0 
                                      
                selTimeDOYs=np.asarray(range(1, numdays.days+1))+midIce
                owInds=selConc<=owThresh                    
                owDOYs=selTimeDOYs[owInds]                    
                difs=np.diff(owDOYs)
                owInds=(difs<2)
                
                # we only want the analysis of first and last day if there are
                # more than a certain number of SI days each year. like 60
                if numFullSI>=0:#60: 
                    if sum(owInds)>0:
                        tempInd=list(np.where(owInds==True)[0]+1)	
                        owDOYs=owDOYs[tempInd]
    
                    if len(owDOYs)>0:
                        firstOW=owDOYs[0]
                    
                        lastOW=owDOYs[-1]-1   
                        
                                             			
                        difs=np.diff(owDOYs)                        
                        earlySeason=((difs>1)*(owDOYs[1:]<midSIF))
    
                        if sum(earlySeason)>0:
                            eSInds=np.where(earlySeason==True)
                            eSInds=eSInds[0]+1	
                            if eSInds[-1]<len(owDOYs):		
                                firstOWc=owDOYs[eSInds[-1]]
                            else:
                                firstOWc=firstOW
                        else:
                            firstOWc=firstOW				
                            lateSeason=((difs>1)*(owDOYs[1:]>midSIF))                        
                        if sum(lateSeason)>0:
                            lSInds=np.where(lateSeason==True)
                            lSInds=lSInds[0]-1
                            if lSInds[0]<len(owDOYs):
                                lastOWc=owDOYs[lSInds[0]]
                            else:
                                lastOWc=lastOW
                        else:
                            lastOWc=lastOW   
                            
                        # save the values:                                              
                        firstDOY[i, nji, nii]=firstOW
                        lastDOY[i, nji, nii]=lastOW
                        firstDOYcontinous[i, nji, nii]=firstOWc
                        lastDOYcontinous[i, nji, nii]=lastOWc
                        springMix[i, nji, nii]=firstOWc-firstOW
                        fallMix[i, nji, nii]=lastOW-lastOWc
                        durationContinous[i, nji, nii]=lastOW-firstOW
                        
                numSIF[i, nji, nii]=len(owDOYs) 
    
        # after looping through each year, mask
        numSIF[i, landMask]=fillVal
        firstDOY[i,landMask]=fillVal
        lastDOY[i, landMask]=fillVal
        firstDOYcontinous[i, landMask]=fillVal
        lastDOYcontinous[i,landMask]=fillVal
        springMix[i, landMask]=fillVal
        fallMix[i, landMask]=fillVal
        durationContinous[i, landMask]=fillVal
        
        numSIF[i, notIceAffected]=fillVal
        firstDOY[i,notIceAffected]=fillVal
        lastDOY[i, notIceAffected]=fillVal
        firstDOYcontinous[i, notIceAffected]=fillVal
        lastDOYcontinous[i,notIceAffected]=fillVal
        springMix[i, notIceAffected]=fillVal
        fallMix[i, notIceAffected]=fillVal
        durationContinous[i, notIceAffected]=fillVal  
                    
from scipy import stats
# calculate trends in NSIF
slope=fillVal*np.ones((nj, ni))
R=fillVal*np.ones((nj, ni))
pval=fillVal*np.ones((nj, ni))
slope_masked=fillVal*np.ones((nj, ni))
for nii in range(ni):
    print nii, '/', ni
    for nji in range(nj):
        if iceAffected[nji, nii]==True: # not land
            s, icpts, r, p, ster = stats.linregress(yearlyTimestep/365.,numSIF[:, nji, nii])               
            slope[nji, nii]=s
            R[nji, nii]=r
            pval[nji, nii]=p
            if p<0.05:
                slope_masked[nji, nii]=s

# 2.5 Save this to a netCDF

fnAnalysis=fn[:-2]+'Analysis.nc'
fAn=Dataset(fnAnalysis, 'w',format='NETCDF4')

# create all the dimentions, set time to unlimited 
for k in f.dimensions.keys():
    if f.unlimited(k)==True:
        fAn.createDimension(k, None)
    else:
        fAn.createDimension(k, f.dimensions[k])
        
# use the netCDF4 instead of pyNIO since it seems to work much better with unlimited variables      
fAnVars={}
for key in {'TLAT', 'TLON','time_bounds', 'time'}:
    print 'creating ', key
    # the netCDF4 module requires that if a fill value exists, it must be declared when the variable is created.
    try:
        fAnVars[key]=fAn.createVariable(key, f.variables[key].typecode(), f.variables[key].dimensions, fill_value=f.variables[key].__dict__['_FillValue'])
    except:
        fAnVars[key]=fAn.createVariable(key, f.variables[key].typecode(), f.variables[key].dimensions)
    # sett all the attribute keys.
    atts = f.variables[key].__dict__
    for attKey in atts.keys():
        if attKey != '_FillValue':
            setattr(fAn.variables[key],attKey,atts[attKey])  

# put data into variables, first the ones we are copying over. 
print 'putting data into standard variables'
for key in {'TLAT', 'TLON'}:
    fAnVars[key][:,:]=f.variables[key][:]
# now, the ones we have created (those with time as a dimention)
fAnVars['time'][:]=yearlyTimestep
fAnVars['time_bounds'][:,:]=yearlyTimeBounds

dataKey=['numSIF','firstDOY','lastDOY','firstDOYcontinous','lastDOYcontinous','springMix','fallMix','durationContinous']
data=[numSIF,firstDOY,lastDOY,firstDOYcontinous,lastDOYcontinous,springMix,fallMix,durationContinous]
units=['number of days', 'day of year', 'day of year', 'day of year', 'day of year','number of days','number of days','number of days']
longName=['Total number of open water days per year', 
'First day of open water',
'Last day of open water',
'First day of continous open water',
'Last day of continuous open water',
'Spring mixed period',
'Fall mixed period',
'Number of days of continous open water']

stdAttributes={'_FillValue': np.array([  1.00000002e+30], dtype=float),
'cell_measures': 'area: tarea',
'cell_methods': 'time: mean',
'comment': 'none',
'coordinates': 'TLON TLAT time',
'missing_value': np.array([  1.00000002e+30], dtype=float),
'time_rep': 'averaged'}

for i in range(len(dataKey)):
    key=dataKey[i]
    print 'creating ', key
    # the netCDF4 module requires that if a fill value exists, it must be declared when the variable is created.
    fAnVars[key]=fAn.createVariable(key, 'f', f.variables[satKey].dimensions, fill_value=stdAttributes['_FillValue'])
    
    # set all the attribute keys.
    for attKey in stdAttributes.keys():
        if attKey != '_FillValue':
            setattr(fAn.variables[key],attKey,stdAttributes[attKey])
    setattr(fAn.variables[key],'long_name',longName[i])
    setattr(fAn.variables[key],'units',units[i])
    fAnVars[key][:,:,:]=data[i] 
     
cKey='pval'
fAnVars[cKey]=fAn.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fAn.variables[cKey],'long_name','Trend pValue') 
setattr(fAn.variables[cKey],'units','')   
setattr(fAn.variables[cKey], 'coordinates', 'nm TLON TLAT')

cKey='R'
fAnVars[cKey]=fAn.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fAn.variables[cKey],'long_name','R2') 
setattr(fAn.variables[cKey],'units','')
setattr(fAn.variables[cKey], 'coordinates', 'nm TLON TLAT')

cKey='slope'
fAnVars[cKey]=fAn.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fAn.variables[cKey],'long_name','Trend Slope') 
setattr(fAn.variables[cKey],'units','Days per Year')
setattr(fAn.variables[cKey], 'coordinates', 'nm TLON TLAT')

fAnVars['slope'][:,:]=slope
fAnVars['pval'][:,:]=pval
fAnVars['R'][:,:]=R   

cKey='slope_masked'
fAnVars[cKey]=fAn.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fAn.variables[cKey],'long_name','Trend Slope (masked by p value)') 
setattr(fAn.variables[cKey],'units','Days per Year')
setattr(fAn.variables[cKey], 'coordinates', 'nm TLON TLAT')
fAnVars['slope_masked'][:,:]=slope_masked



fAn.close()  

  
    
      
        
            
f.close()
del f
del fAn
print'finished analysis of ',fn , (datetime.now()-startTime)
#     
### There are some other things to consider. Like prediction of open water.  
## A basic trendline might be a super lamo way to show things once persistant open water occurs
## by what year is persistant open water predicted
#
# where do models agree and where do they disagree?
	# how does SD of metrics change through time. 
	
	# plot dur vs time for drew point for each of the happy ensembles. 
	
# add something with closest sea ice to each part of the coast (implications for polar bears)
# SST
# implications for biogeochem and permafrost. 

pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print s.getvalue()	
	