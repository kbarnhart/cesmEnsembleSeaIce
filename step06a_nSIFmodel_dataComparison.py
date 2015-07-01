#model obs compariosn

#for each year. for each obs cell. get the correspoding large cell and find what percentile the observation is in
#(we hope for consistent fifties)

# do all metrics:
print 'importing modules'
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import glob
import haversine as h

path=u'/Volumes/Svalbard/seaIcePPF/'

dirList=glob.glob(path+'b.e11.f09_g16.001.*.timeseries.Analysis.nc')
dirListOrig=glob.glob(path+'b.e11.f09_g16.001.*.timeseries.nc')
obs=u'/Volumes/Svalbard/seaIcePPF/satelliteObs.timeseries.nc'
obsA=u'/Volumes/Svalbard/seaIcePPF/satelliteObs.timeseries.Analysis.nc'

startTime = datetime.now()

fMod=Dataset(dirList[0], 'r')
fModOrig=Dataset(dirListOrig[0], 'r')
fObs=Dataset(obs, 'r')
fObsA=Dataset(obsA, 'r')

files=[]
for fn in dirList:
    files.append(Dataset(fn, 'r'))
    print fn

    
numModels=len(dirList)
numYears=(2100-1850)

ni=320
nj=104

numSIF=np.nan*np.ones((numModels,numYears,nj,ni))
first=np.nan*np.ones((numModels,numYears,nj,ni))
last=np.nan*np.ones((numModels,numYears,nj,ni))

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

obsValnSIF=fObsA.variables['numSIF'][:,:,:].data
obsValfirst=fObsA.variables['firstDOY'][:,:,:].data
obsVallast=fObsA.variables['lastDOY'][:,:,:].data


aice=fModOrig.variables['aice_d'][0,:,:]
obice=fObs.variables['satelliteSIC'][0,:,:]

lats=fMod.variables['TLAT'][:,:]
lons=fMod.variables['TLON'][:,:]

latsOB=fObs.variables['TLAT'][:,:]
lonsOB=fObs.variables['TLON'][:,:]

latSelect=lats
latSelect[aice>100]=np.nan

lonSelect=lons
lonSelect[aice>100]=np.nan

latSelectOB=latsOB
latSelectOB[obice>100]=np.nan

lonSelectOB=lonsOB
lonSelectOB[obice>100]=np.nan

njo=448
nio=304

fillVal=fModOrig.variables['aice_d'].__dict__['_FillValue']

numObsYears=fObsA.variables['numSIF'][:,0,0].size

nSIFperc=fillVal*np.ones(fObsA.variables['numSIF'][:,:,:].shape)
firstPerc=fillVal*np.ones(fObsA.variables['numSIF'][:,:,:].shape)
lastPerc=fillVal*np.ones(fObsA.variables['numSIF'][:,:,:].shape)

time=np.zeros(numObsYears)
timeBounds=np.zeros((numObsYears,2))
offs=1979-1850
modelTime=np.zeros(numObsYears)
for i in range(numObsYears):
    time[i]=fObsA.variables['time'][i]   
    timeBounds[i,0]=fObsA.variables['time_bounds'][i,0]
    timeBounds[i,1]=fObsA.variables['time_bounds'][i,1]
    modelTime[i]=fMod.variables['time'][i+offs]/365.
    
print 'starting analysis'
ocean=obice.mask==False

dataset={}
for njj in range(njo):
    print njj,'/', njo, (datetime.now()-startTime)
    for nii in range(nio):
        
        if ocean[njj, nii]==True:
            
            sLat=latSelectOB[njj, nii]
            sLon=lonSelectOB[njj, nii]

            dists=h.haversine(sLon, sLat, lonSelect, latSelect)
            ind=np.where(dists==np.nanmin(dists))            
            nj=ind[0][0]
            ni=ind[1][0]
            for i in range(numObsYears):
                
                # mask appropriately. 
                if obsValnSIF[i,njj,nii]<800:
                    modelValsSIF = numSIF[:,i+offs,nj,ni]
                    modelValsSIF=modelValsSIF[modelValsSIF<800]  
                           
                    if modelValsSIF.size>25:
                        numOver=np.sum(modelValsSIF>obsValnSIF[i,njj,nii])
                        numUnder=np.sum(modelValsSIF<=obsValnSIF[i,njj,nii])
                        
                    
                        perc=(float(numUnder)/modelValsSIF.size)*0.5+((modelValsSIF.size-float(numOver))/modelValsSIF.size)*0.5 # linear interp between percentiles
                        
                        if numOver==modelValsSIF.size or numUnder==modelValsSIF.size:
                            perc=0.5
                        
                    nSIFperc[i, njj, nii]=perc
                
                if obsValfirst[i,njj,nii]<800:
                    modelValsFirst = first[:,i+offs,nj,ni]
                    modelValsFirst=modelValsFirst[modelValsFirst<800]
                    
                    if modelValsFirst.size>25:
                        numOver=np.sum(modelValsFirst>obsValfirst[i,njj,nii])
                        numUnder=np.sum(modelValsFirst<=obsValfirst[i,njj,nii])
                        percFirst=(float(numUnder)/modelValsFirst.size)*0.5+((modelValsFirst.size-float(numOver))/modelValsFirst.size)*0.5
                        if numOver==modelValsFirst.size or numUnder==modelValsFirst.size:
                            percFirst=0.5
                        firstPerc[i, njj, nii]=percFirst
                    
                if obsVallast[i,njj,nii]<800:        
                    modelValsLast = last[:,i+offs,nj,ni]
                    modelValsLast=modelValsLast[modelValsLast<800]
                    
                    if modelValsLast.size>25:
                        numOver=np.sum(modelValsLast>obsVallast[i,njj,nii])
                        numUnder=np.sum(modelValsLast<=obsVallast[i,njj,nii])
                    
                        percLast=(float(numUnder)/modelValsLast.size)*0.5+((modelValsLast.size-float(numOver))/modelValsLast.size)*0.5
                        if numOver==modelValsLast.size or numUnder==modelValsLast.size:
                            percLast=0.5
                        lastPerc[i, njj, nii]=percLast        
    
dataset['nSIFperc']={'data':nSIFperc}
dataset['nSIFperc']['units']='percentile'
dataset['nSIFperc']['longName']='Percentile of number of sea ice free days observed'

dataset['firstPerc']={'data':firstPerc}
dataset['firstPerc']['units']='percentile'
dataset['firstPerc']['longName']='Percentile of first day of open water observed'

dataset['lastPerc']={'data':lastPerc}
dataset['lastPerc']['units']='percentile'
dataset['lastPerc']['longName']='Percentile of last day of open water observed'

# save as a netcdf

import nio
f=nio.open_file(obs, 'r')

fnAnalysis='/Volumes/Svalbard/seaIcePPF/obsModelComparison.nc'
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
fAnVars['time'][:]=time
fAnVars['time_bounds'][:,:]=timeBounds

dataKey=dataset.keys()

stdAttributes={'_FillValue': np.array([  1.00000002e+30], dtype=float),
'cell_measures': 'area: tarea',
'cell_methods': 'time: mean',
'comment': 'none',
'coordinates': 'TLON TLAT time',
'missing_value': np.array([  1.00000002e+30], dtype=float),
'time_rep': 'averaged'}

for i in range(len(dataKey)):
    key=dataKey[i]
    d=dataset[key]
    data=d['data']
    print 'creating ', key
    # the netCDF4 module requires that if a fill value exists, it must be declared when the variable is created.
    if len(data.shape)==3:
        fAnVars[key]=fAn.createVariable(key, 'f',  ('time', 'nj', 'ni'), fill_value=stdAttributes['_FillValue'])
        fAnVars[key][:,:,:]=data
        ('time', 'nj', 'ni')
    elif len(data.shape)==2:
        fAnVars[key]=fAn.createVariable(key, 'f',  ('nj', 'ni'), fill_value=stdAttributes['_FillValue'])
        fAnVars[key][:,:]=data
            
    # set all the attribute keys.
    for attKey in stdAttributes.keys():
        if attKey != '_FillValue':
            setattr(fAn.variables[key],attKey,stdAttributes[attKey])
    setattr(fAn.variables[key],'long_name',d['longName'])
    setattr(fAn.variables[key],'units',d['units'])

fAn.close()    
f.close()
del f
del fAn