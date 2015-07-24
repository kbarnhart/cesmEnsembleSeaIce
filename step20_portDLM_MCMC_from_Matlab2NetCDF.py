
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec
from netCDF4 import Dataset
import nio
from scipy import stats
from datetime import datetime

#import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['font.sans-serif']=['Arial']

startTime = datetime.now()

import glob


fn_monthOut=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/DLM__MCMC_results_final.nc'


## get the Sea ice days info 
path=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/'
#path=u'/Users/katherinebarnhart/Desktop/RESEARCH/seaIcePPF/'
#os.chdir(path)
dirList=glob.glob(path+'b.e11.B20TRC5CNBDRD-BRCP85C5CNBDRD.f09_g16.*nh.Analysis.nc')
dirList2=glob.glob(path+'b.e11.B1850C5CN.f09_g16.005.cice.h1.aice_d_nh.*.Analysis.nc')#path=u'/Users/katherinebarnhart/Desktop/RESEARCH/seaIcePPF/'

analysisFiles=[]    
for fn in dirList:
    analysisFiles.append(Dataset(fn, 'r'))
    print fn, (datetime.now()-startTime)
    
analysisFiles2=[]    
for fn in dirList2:
    analysisFiles2.append(Dataset(fn, 'r'))
    print fn, (datetime.now()-startTime)    
    
numModels=len(dirList)
numYears=(2100-1850)

NI=320
NJ=104
NM=30

key='aice_d'
f=nio.open_file(u'/Volumes/Pitcairn/seaicePPF/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h1.aice_d_nh.19200101-20051231.nc','r')
landMask=f.variables[key][0,:,:]>120
fillVal=f.variables[key].__dict__['_FillValue']
del f


numSIF=np.nan*np.ones((numModels,numYears,NJ,NI))
first=np.nan*np.ones((numModels,numYears,NJ,NI))
last=np.nan*np.ones((numModels,numYears,NJ,NI))

i=0
for fn in dirList:
    f=Dataset(fn, 'r')
    
    ns=f.variables['numSIF'][:,:,:]
    numSIF[i,-ns.shape[0]:,:,:]=ns
    
    fs=f.variables['firstDOY'][:,:,:]
    first[i,-ns.shape[0]:,:,:]=fs
    
    ls=f.variables['lastDOY'][:,:,:]
    last[i,-ns.shape[0]:,:,:]=ls
    
    print fn, (datetime.now()-startTime)
    i+=1
    f.close()

nSIF=np.array(numSIF)
first=np.array(first)
last=np.array(last)


                                                                                                                                                                                                                          
numSIF1850=[]
fir1850=[]
las1850=[]

for fn in dirList2:
    f=Dataset(fn, 'r')
    numSIF1850.extend(f.variables['numSIF'][:,:,:])
    fir1850.extend(f.variables['firstDOY'][:,:,:])
    las1850.extend(f.variables['lastDOY'][:,:,:])
    print fn, (datetime.now()-startTime) 
    f.close()   
    
nSIF1850=np.array(numSIF1850)
first1850=np.array(fir1850)
last1850=np.array(las1850)
del numSIF1850
del fir1850
del las1850


#### now get the output from the DLM
path=u'/Users/katherinebarnhart/MATLABwork/seaIceEmergence/output_v3/'
dirListYear=glob.glob(path+'year*.csv')


dirListLevelMean=glob.glob(path+'levelmean*.csv')
dirListLevelStd=glob.glob(path+'levelstd*.csv')
dirListSlopeMean=glob.glob(path+'slopemean*.csv')
dirListSlopeStd=glob.glob(path+'slopestd*.csv')

NJ=104
NI=320
NM=30

e_year_mmean=np.nan*np.ones((NJ, NI))
e_year=np.nan*np.ones((NJ, NI))
s_year=np.nan*np.ones((NJ, NI))

e_year_bgmean=np.nan*np.ones((NJ, NI))
slope=np.nan*np.ones((NM, NJ, NI))

levelmean=np.nan*np.ones((numYears,NJ,NI))
levelstd=np.nan*np.ones((numYears,NJ,NI))
slopestd=np.nan*np.ones((numYears,NJ,NI))
slopemean=np.nan*np.ones((numYears,NJ,NI))

# get values from dlm output
for rFN in dirListYear:

    rFile = open(rFN, 'r')
    print rFN
    rData=rFile.readlines()


    for i in range(1,len(rData)):
        #print i
        line=rData[i].strip().split(',')
        ni=int(line[0])-1# adjust from 1 based index to 0 based
        nj=int(line[1])-1
        
        s_year[nj, ni]=int(line[2])
        e_year[nj, ni]=int(line[4])
        e_year_mmean[nj, ni]=int(line[3])
    rFile.close()
    del rFile


for rFN in dirListLevelMean:

    rFile = open(rFN, 'r')
    print rFN
    rData=rFile.readlines()


    for i in range(1,len(rData)):
        #print i
        line=rData[i].strip().split(',')
        ni=int(line[0])-1# adjust from 1 based index to 0 based
        nj=int(line[1])-1
        
        rest=line[2:]
        levelmean[-len(rest):, nj, ni]=rest
    rFile.close()
    del rFile


for rFN in dirListLevelStd:

    rFile = open(rFN, 'r')
    print rFN
    rData=rFile.readlines()


    for i in range(1,len(rData)):
        #print i
        line=rData[i].strip().split(',')
        ni=int(line[0])-1# adjust from 1 based index to 0 based
        nj=int(line[1])-1
        
        rest=line[2:]
        levelstd[-len(rest):, nj, ni]=rest
    rFile.close()
    del rFile

for rFN in dirListSlopeMean:

    rFile = open(rFN, 'r')
    print rFN
    rData=rFile.readlines()


    for i in range(1,len(rData)):
        #print i
        line=rData[i].strip().split(',')
        ni=int(line[0])-1# adjust from 1 based index to 0 based
        nj=int(line[1])-1
        
        rest=line[2:]
        slopemean[-len(rest):, nj, ni]=rest
    rFile.close()
    del rFile


for rFN in dirListSlopeStd:

    rFile = open(rFN, 'r')
    print rFN
    rData=rFile.readlines()


    for i in range(1,len(rData)):
        #print i
        line=rData[i].strip().split(',')
        ni=int(line[0])-1# adjust from 1 based index to 0 based
        nj=int(line[1])-1
        
        rest=line[2:]
        slopestd[-len(rest):, nj, ni]=rest
    rFile.close()
    del rFile

lowerBound=levelmean-2.*levelstd
upperBound=levelmean+2.*levelstd



isnan=np.isnan(levelmean)
levelmean[isnan]=fillVal 

isnan=np.isnan(levelstd)
levelstd[isnan]=fillVal 

isnan=np.isnan(slopestd)
slopestd[isnan]=fillVal 

isnan=np.isnan(slopemean)
slopemean[isnan]=fillVal 

isnan=np.isnan(lowerBound)
lowerBound[isnan]=fillVal

isnan=np.isnan(upperBound)
upperBound[isnan]=fillVal

# now calculate the post-emergence slope

# use the year based on 
emerge=e_year

slope=np.nan*np.ones((NM, NJ, NI))
slope_masked=np.nan*np.ones((NM, NJ, NI))
R=np.nan*np.ones((NM, NJ, NI))
pval=np.nan*np.ones((NM, NJ, NI)) 
icemask=np.nanmean(nSIF1850, axis=0)<350

time=np.arange(1850,2101)

for nii in range(NI):
    print 'calculating trendlines', str(nii)+'/'+str(NI), (datetime.now()-startTime)
    for njj in range(NJ):
        if  (icemask[njj, nii]==True): 
            for nmm in range(NM):             
                                
                #### choose last year that there is less than 350 days per year of open water
                
                fyind=np.where(nSIF[nmm,:, njj, nii]<350)[0]
                if len(fyind)==0:
                    stopInd=250
                else:
                    stopInd=fyind[-1]
                postEmergeInd=np.where(time>emerge[njj, nii])[0]
                
                
                if len(postEmergeInd)>0:
                    startInd=postEmergeInd[0]
                    
                    if stopInd>startInd+4:
                        s, icpts, r, p, ster = stats.linregress(time[startInd:stopInd],nSIF[nmm,startInd:stopInd, njj, nii])               
                        slope[nmm, njj, nii]=s
                        R[nmm, njj, nii]=r
                        pval[nmm, njj, nii]=p 
                        if p<0.05:
                            slope_masked[nmm, njj, nii]=s
                                                                                                                                                                           
meanSlope=np.nanmean(slope, axis=0)
stdSlope=np.nanstd(slope, axis=0)                                           

meanSlope_masked=np.nanmean(slope_masked, axis=0)
stdSlope_masked=np.nanstd(slope_masked, axis=0)                                                    
        
isnan=np.isnan(slope_masked)
slope_masked[isnan]=fillVal  
isinf=np.isinf(slope_masked)
slope_masked[isinf]=fillVal                                                                                                                                                      
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
isnan=np.isnan(slope)
slope[isnan]=fillVal  
isinf=np.isinf(slope)
slope[isinf]=fillVal       
    
print 'making the netcdf'

modelFile= '/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/allModels.Summary.nh.RCP85.nc'
f=nio.open_file(modelFile, 'r')

## Create this monthly averaged file as a new netcdf
fMonth=Dataset(fn_monthOut, 'w',format='NETCDF4')
fillVal=f.variables['nSIFMean'].__dict__['_FillValue'][0]

# create all the dimentions, set time to unlimited 
for k in f.dimensions.keys():
    if f.unlimited(k)==True:
        fMonth.createDimension(k, None)
    else:
        fMonth.createDimension(k, f.dimensions[k])

                
# use the netCDF4 instead of pyNIO since it seems to work much better with unlimited variables      
fMonthVars={}
for key in {'TLAT', 'TLON','latt_bounds','lont_bounds', 'time', 'time_bounds'}:
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
for key in {'TLAT', 'TLON','latt_bounds','lont_bounds', 'time_bounds'}:       
    fMonthVars[key][:,:]=f.variables[key][:]
fMonthVars['time'][:]=f.variables['time'][:]    


satKey='numSIF'     

#fMonth.createDimension('nm', 30)
   
cKey='nSIF'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nm', 'time','nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Number of Sea Ice Free Days (original model output)') 
setattr(fMonth.variables[cKey],'units','days')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['nSIF'][:,:,:,:]=numSIF   
   
   # dlm results
   
   
cKey='dlm_levelmean'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('time','nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Mean level of nSIF in DLM') 
setattr(fMonth.variables[cKey],'units','number of sea ice free days per year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['dlm_levelmean'][:,:,:]=levelmean  
   
cKey='dlm_levelstd'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('time','nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Standard deviation of level of nSIF in DLM') 
setattr(fMonth.variables[cKey],'units','number of sea ice free days per year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['dlm_levelstd'][:,:,:]=levelstd  

cKey='dlm_upperBound'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('time','nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Upper bound of predicted nSIF in DLM') 
setattr(fMonth.variables[cKey],'units','number of sea ice free days per year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['dlm_upperBound'][:,:,:]=upperBound

cKey='dlm_lowerBound'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('time','nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Lower bound of predicted nSIF in DLM') 
setattr(fMonth.variables[cKey],'units','number of sea ice free days per year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['dlm_lowerBound'][:,:,:]=lowerBound  

   
cKey='dlm_slopemean'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('time','nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Mean slope of nSIF in DLM') 
setattr(fMonth.variables[cKey],'units','number of sea ice free days per year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['dlm_slopemean'][:,:,:]=slopemean  
   
cKey='dlm_slopestd'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('time','nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Standard deviation of slope of nSIF in DLM') 
setattr(fMonth.variables[cKey],'units','number of sea ice free days per year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['dlm_slopestd'][:,:,:]=slopestd     
   
   
# post emergence slope
cKey='meanPostEmergenceExpansion'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Mean of Post Emergence Expansion of Open Water') 
setattr(fMonth.variables[cKey],'units','days per year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['meanPostEmergenceExpansion'][:,:]=meanSlope


cKey='stdPostEmergenceExpansion'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Standard Deviation of Slope') 
setattr(fMonth.variables[cKey],'units','days per year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['stdPostEmergenceExpansion'][:,:]=stdSlope
 
                                                                                             
cKey='meanPostEmergenceExpansion_masked'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Mean of Slope (masked by Pval)') 
setattr(fMonth.variables[cKey],'units','days per year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['meanPostEmergenceExpansion_masked'][:,:]=meanSlope_masked

cKey='stdPostEmergenceExpansion_masked'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Standard Deviation of Slope (masked by Pval)') 
setattr(fMonth.variables[cKey],'units','days per year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
fMonthVars['stdPostEmergenceExpansion_masked'][:,:]=stdSlope_masked
                                                                                                                                                                                                                                                  
cKey='pval_PostEmergenceExpansion'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Trend pValue for post emergence expansion of open water') 
setattr(fMonth.variables[cKey],'units','')   
setattr(fMonth.variables[cKey], 'coordinates', 'nm TLON TLAT')
fMonthVars['pval_PostEmergenceExpansion'][:,:,:]=pval

cKey='R_PostEmergenceExpansion'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Trend R2 for post emergence expansion of open water') 
setattr(fMonth.variables[cKey],'units','')
setattr(fMonth.variables[cKey], 'coordinates', 'nm TLON TLAT')
fMonthVars['R_PostEmergenceExpansion'][:,:,:]=R   

cKey='PostEmergenceExpansionRate'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Trend Slope for post emergence expansion of open water (by model)') 
setattr(fMonth.variables[cKey],'units','Days per Year')
setattr(fMonth.variables[cKey], 'coordinates', 'nm TLON TLAT')
fMonthVars['PostEmergenceExpansionRate'][:,:,:]=slope

cKey='PostEmergenceExpansion_masked'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Trend Slope for post emergence expansion of open water (by model) masked by Pval <0.05') 
setattr(fMonth.variables[cKey],'units','Days per Year')
setattr(fMonth.variables[cKey], 'coordinates', 'nm TLON TLAT')
fMonthVars['PostEmergenceExpansion_masked'][:,:,:]=slope_masked

## shift years

cKey='shift_year'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Shift Year') 
setattr(fMonth.variables[cKey],'units','year')  
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT') 
fMonthVars['shift_year'][:,:]=s_year        
                                                                     
cKey='lag_time'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Lag Time (emergence (95% based) minus shift year)') 
setattr(fMonth.variables[cKey],'units','years')  
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT') 
fMonthVars['lag_time'][:,:]=e_year-s_year        
                                                                                                                                                                                                                                                                                                                                                                                                              
cKey='emergence_year'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Emergence Year(based on BG 95th range)') 
setattr(fMonth.variables[cKey],'units','year')  
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT') 
fMonthVars['emergence_year'][:,:]=e_year        


cKey='emergence_year_bgmean'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Emergence Year (based on BG mean)') 
setattr(fMonth.variables[cKey],'units','year')  
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT') 
fMonthVars['emergence_year_bgmean'][:,:]=e_year_bgmean        

cKey='emergence_year_mmean'
fMonthVars[cKey]=fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Emergence Year (model mean)') 
setattr(fMonth.variables[cKey],'units','year')  
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT') 
fMonthVars['emergence_year_mmean'][:,:]=e_year_mmean        
    
                                                     
fMonth.close()   
f.close()