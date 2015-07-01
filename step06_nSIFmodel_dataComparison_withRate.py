#model obs compariosn

#for each year. for each obs cell. get the correspoding large cell and find what percentile the observation is in
#(we hope for consistent fifties)

# do all metrics:
print 'importing modules'
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import glob
import sys
from scipy import stats
sys.path.append(u'/Volumes/Pitcairn/seaicePPF/code/')

import haversine as h


path=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/'

pathIn=u'/Volumes/Pitcairn/seaicePPF/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/'
pathIn45=u'/Volumes/Pitcairn/seaicePPF/p/cesm0005/CESM-CAM5-BGC-ME/ice/proc/tseries/daily/aice_d/'
pathOut=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/'


# hard code in the names for the run parts
bgKey=u'B1850C5CN'
runPart1key=u'B20TRC5CNBDRD'
runParts23keyRCP85=u'BRCP85C5CNBDRD'
runParts23keyRCP45=u'BRCP45C5CNBDRD'

obs=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/satelliteObs.timeseries.nc'
obsA=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/satelliteObs.timeseries.Analysis.nc'


fn_monthOut=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/allModels.Summary.nh.RCP85.nc'

dirList=glob.glob(path+'*BRCP85C5CNBDRD*nh*.Analysis.nc')
dirListOrig=glob.glob(pathIn+'*B1850C5CN*nh*')
fMonth=Dataset(fn_monthOut, 'r',format='NETCDF4')
slope_model=fMonth.variables['slope'][:,:,:]
slope_model_masked=fMonth.variables['slope_masked'][:,:,:]

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

obsRate=fObsA.variables['slope'][:,:].data

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


nSIFdiff=fillVal*np.ones(fObsA.variables['numSIF'][:,:,:].shape)
firstdiff=fillVal*np.ones(fObsA.variables['numSIF'][:,:,:].shape)
lastdiff=fillVal*np.ones(fObsA.variables['numSIF'][:,:,:].shape)

slope_comparison=fillVal*np.ones((numModels, njo, nio))

slope_percentile=fillVal*np.ones((njo, nio))
var_percentile=fillVal*np.ones((njo, nio))
var_diff=fillVal*np.ones((njo, nio))
slope_diff=fillVal*np.ones((njo, nio))

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

# also calculate the absolute difference in day between the model mean and the observations


ocean=obice.mask==False
trendTime=np.array(range(1979,2014))
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

						# slope percentiles. 
						slope_comparison[:, njj,nii]=obsRate[njj,nii]/slope_model_masked[:, nj,ni]
						numValidModles=numModels-sum(np.isnan(slope_model_masked[:, nj,ni]))
						numOverSlope=np.sum(slope_model_masked[:, nj,ni]>obsRate[njj,nii])
						numUnderSlope=np.sum(slope_model_masked[:, nj,ni]<=obsRate[njj,nii])
						percSlope=(float(numUnderSlope)/numValidModles)*0.5+((numValidModles-float(numOverSlope))/numValidModles)*0.5 # linear interp between percentiles
						if numOverSlope==numValidModles or numUnderSlope==numValidModles:
								percSlope=0.5
						slope_percentile[njj, nii]=percSlope

						obsN=obsValnSIF[:,njj,nii] # determine the percentile of detrended variance
						sN, iN, rtemp, ptemp, setemp = stats.linregress(trendTime,obsN)
						residObs=obsN-(iN+(sN*trendTime))

						obsVar=np.std(residObs)
						modVars=[]
						modN = numSIF[:,i+offs:i+offs+numObsYears,nj,ni]
						for varItter in range(30):
							sN, iN, rtemp, ptemp, setemp = stats.linregress(trendTime,modN[varItter,:])
							residMod=modN[varItter,:]-(iN+(sN*trendTime))
							modVars.append(np.std(residMod))
						modelVars=np.array(modVars)
						numOverVar=np.sum(modelVars>obsVar)
						numUnderVar=np.sum(modelVars<=obsVar)
						percVar=(float(numUnderVar)/30.)*0.5+((30.-float(numOverVar))/30.)*0.5 # linear interp between percentiles
						if numOverVar==30 or numUnderVar==30:
								percVar=0.5
						var_percentile[njj, nii]=percVar

						var_diff[njj, nii]=np.abs(obsVar-np.nanmean(modelVars))
						slope_diff[njj, nii]=np.abs(obsRate[njj,nii]-np.nanmean(slope_model_masked[:, nj,ni]))
						
						for i in range(numObsYears):    
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
										nSIFdiff[i, njj, nii]=np.abs(obsValnSIF[i,njj,nii]-np.nanmean(modelValsSIF))
										
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
												firstdiff[i, njj, nii]=np.abs(obsValfirst[i,njj,nii]-np.nanmean(modelValsFirst))
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
    											lastdiff[i, njj, nii]=np.abs(obsVallast[i,njj,nii]-np.nanmean(modelValsLast))
dataset['nSIFperc']={'data':nSIFperc}
dataset['nSIFperc']['units']='percentile'
dataset['nSIFperc']['longName']='Percentile of number of sea ice free days observed'

dataset['firstPerc']={'data':firstPerc}
dataset['firstPerc']['units']='percentile'
dataset['firstPerc']['longName']='Percentile of first day of open water observed'

dataset['lastPerc']={'data':lastPerc}
dataset['lastPerc']['units']='percentile'
dataset['lastPerc']['longName']='Percentile of last day of open water observed'

dataset['nSIFdiff']={'data':nSIFperc}
dataset['nSIFdiff']['units']='days'
dataset['nSIFdiff']['longName']='Difference between model mean and observed of number of sea ice free days'

dataset['firstdiff']={'data':firstPerc}
dataset['firstdiff']['units']='days'
dataset['firstdiff']['longName']='Difference between model mean and observed first day of open water'

dataset['lastdiff']={'data':lastPerc}
dataset['lastdiff']['units']='days'
dataset['lastdiff']['longName']='Difference between model mean and observed last day of open water'

dataset['var_percentile']={'data':var_percentile}
dataset['var_percentile']['longName']='Comparison of variance'
dataset['var_percentile']['units']='percentile of 1979-2014 detrended variance'

dataset['slope_percentile']={'data':slope_percentile}
dataset['slope_percentile']['longName']='Comparison of slope'
dataset['slope_percentile']['units']='percentile of 1979-2014 slope'

dataset['var_diff']={'data':var_diff}
dataset['var_diff']['longName']='Difference of model mean and observed 1979-2014 detrended variance'
dataset['var_diff']['units']='Days'

dataset['slope_diff']={'data':slope_diff}
dataset['slope_diff']['longName']='Difference of model mean and observed 1979-2014 trend'
dataset['slope_diff']['units']='Days per year'


# save as a netcdf

import nio
f=nio.open_file(obs, 'r')

fnAnalysis='/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/obsModelComparison.nc'
fAn=Dataset(fnAnalysis, 'w',format='NETCDF4')

# create all the dimentions, set time to unlimited 
for k in f.dimensions.keys():
    if f.unlimited(k)==True:
        fAn.createDimension(k, None)
    else:
        fAn.createDimension(k, f.dimensions[k])
fAn.createDimension('nm', numModels)        
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

			
		isnan=np.isnan(data)
		data[isnan]=fillVal  
		isinf=np.isinf(data)
		data[isinf]=fillVal 

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
    
    
cKey='slope_comparison'
fAnVars[cKey]=fAn.createVariable(cKey, 'f', ('nm','nj','ni'),fill_value=fillVal)
setattr(fAn.variables[cKey],'long_name','Comparison of slopes') 
setattr(fAn.variables[cKey],'units','observed 1979-2014 slope divided by modeled 1979-2014 slope')
setattr(fAn.variables[cKey], 'coordinates', 'TLON TLAT')
fAnVars['slope_comparison'][:,:,:]=slope_comparison

fAn.close()    
f.close()
del f
del fAn