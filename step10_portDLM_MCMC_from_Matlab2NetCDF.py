
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

## get the Sea ice days info 
path=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/'
#path=u'/Users/katherinebarnhart/Desktop/RESEARCH/seaIcePPF/'
#os.chdir(path)
dirList=glob.glob(path+'b.e11.f09_g16.001.*.timeseries.Analysis.nc')
dirList2=glob.glob(path+'b.e11.B1850C5CN.*.timeseries.Analysis.nc')#path=u'/Users/katherinebarnhart/Desktop/RESEARCH/seaIcePPF/'

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
f=nio.open_file(u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/b.e11.f09_g16.001.cice.h1.aice_d_nh.002.timeseries.nc','r')
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
dirList=glob.glob(path+'year*.csv')


NJ=104
NI=320
NM=30

e_year_mmean=np.nan*np.ones((NJ, NI))
e_year=np.nan*np.ones((NJ, NI))
s_year=np.nan*np.ones((NJ, NI))

e_year_bgmean=np.nan*np.ones((NJ, NI))
slope=np.nan*np.ones((NM, NJ, NI))


for rFN in dirList:

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

modelFile= u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/b.e11.f09_g16.001.cice.h1.aice_d_nh.001.timeseries.Analysis.nc'
f=nio.open_file(modelFile, 'r')

## Create this monthly averaged file as a new netcdf
fn_monthOut=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/DLM__MCMC_results_v2.nc'
fMonth=Dataset(fn_monthOut, 'w',format='NETCDF4')
fillVal=f.variables['numSIF'].__dict__['_FillValue'][0]

# create all the dimentions, set time to unlimited 
for k in f.dimensions.keys():
    if f.unlimited(k)==True:
        fMonth.createDimension(k, None)
    else:
        fMonth.createDimension(k, f.dimensions[k])

                
# use the netCDF4 instead of pyNIO since it seems to work much better with unlimited variables      
fMonthVars={}
for key in {'TLAT', 'TLON','latt_bounds','lont_bounds'}:
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
for key in {'TLAT', 'TLON','latt_bounds','lont_bounds'}:       
    fMonthVars[key][:,:]=f.variables[key][:]
    
satKey='numSIF'     

fMonth.createDimension('nm', 30)
   
cKey='meanSlope'
fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Mean of Slope') 
setattr(fMonth.variables[cKey],'units','year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')

cKey='stdSlope'
fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Standard Deviation of Slope') 
setattr(fMonth.variables[cKey],'units','year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
 
                                                                                             
cKey='meanSlope_masked'
fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Mean of Slope (masked by Pval)') 
setattr(fMonth.variables[cKey],'units','year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')

cKey='stdSlope_masked'
fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Standard Deviation of Slope (masked by Pval)') 
setattr(fMonth.variables[cKey],'units','year')   
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT')
                                                                                                                                                                                                                                                  

cKey='pval'
fMonth.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Trend pValue') 
setattr(fMonth.variables[cKey],'units','')   
setattr(fMonth.variables[cKey], 'coordinates', 'nm TLON TLAT')

cKey='R'
fMonth.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','R2') 
setattr(fMonth.variables[cKey],'units','')
setattr(fMonth.variables[cKey], 'coordinates', 'nm TLON TLAT')

cKey='slope'
fMonth.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Trend Slope') 
setattr(fMonth.variables[cKey],'units','Days per Year')
setattr(fMonth.variables[cKey], 'coordinates', 'nm TLON TLAT')

cKey='slope_masked'
fMonth.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Trend Slope Masked by Pval <0.05') 
setattr(fMonth.variables[cKey],'units','Days per Year')
setattr(fMonth.variables[cKey], 'coordinates', 'nm TLON TLAT')


fMonth.variables['slope_masked'][:,:,:]=slope_masked
fMonth.variables['slope'][:,:,:]=slope
fMonth.variables['meanSlope'][:,:]=meanSlope
fMonth.variables['stdSlope'][:,:]=stdSlope

fMonth.variables['meanSlope_masked'][:,:]=meanSlope_masked
fMonth.variables['stdSlope_masked'][:,:]=stdSlope_masked

fMonth.variables['pval'][:,:,:]=pval
fMonth.variables['R'][:,:,:]=R   

## shift years

cKey='shift_year'
fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Shift Year') 
setattr(fMonth.variables[cKey],'units','year')  
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT') 
                                                                     
cKey='lag_time'
fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Lag Time (emergence-shift)') 
setattr(fMonth.variables[cKey],'units','year')  
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT') 
                                                                                                                                                                                                                                                                                                                                                                                                              
cKey='emergence_year'
fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Emergence Year(based on BG 95th range)') 
setattr(fMonth.variables[cKey],'units','year')  
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT') 

cKey='emergence_year_bgmean'
fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Emergence Year (based on BG mean)') 
setattr(fMonth.variables[cKey],'units','year')  
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT') 

cKey='emergence_year_mmean'
fMonth.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
setattr(fMonth.variables[cKey],'long_name','Emergence Year (model mean)') 
setattr(fMonth.variables[cKey],'units','year')  
setattr(fMonth.variables[cKey], 'coordinates', 'TLON TLAT') 

fMonth.variables['emergence_year'][:,:]=e_year        
fMonth.variables['shift_year'][:,:]=s_year        

fMonth.variables['lag_time'][:,:]=e_year-s_year        
fMonth.variables['emergence_year_bgmean'][:,:]=e_year_bgmean        
fMonth.variables['emergence_year_mmean'][:,:]=e_year_mmean        
    
                                                     
fMonth.close()   
f.close()