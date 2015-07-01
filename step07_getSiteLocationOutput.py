print 'importing modules'
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import glob
import pickle
import haversine as h


path=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/'
#path=u'/Users/katherinebarnhart/Desktop/RESEARCH/seaIcePPF/'
#os.chdir(path)
dirList=glob.glob(path+'b.e11.f09_g16.001.*.timeseries.Analysis.nc')
dirListOrig=glob.glob(path+'b.e11.f09_g16.001.*.timeseries.nc')

dirList1850=glob.glob(path+'b.e11.B1850C5CN.*.timeseries.Analysis.nc')
dirListOrig1850=glob.glob(path+'b.e11.B1850C5CN.*.timeseries.nc')

obs=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/satelliteObs.timeseries.nc'
obsA=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/satelliteObs.timeseries.Analysis.nc'

startTime = datetime.now()

# drew point is at (70.8, -157)
# this corresponds to 
f=Dataset(dirListOrig[0], 'r')
fObs=Dataset(obs, 'r')
fObsA=Dataset(obsA, 'r')

aice=f.variables['aice_d'][0,:,:]
obice=fObs.variables['satelliteSIC'][0,:,:]

lats=f.variables['TLAT'][:,:]
lons=f.variables['TLON'][:,:]

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

siteOutput={}
siteOutput['Barrow']={'loc':(71.2956, -156.7664)}
siteOutput['Kamchatka']={'loc':(61.840594, 173.970393)}
siteOutput['Kotzebue']={'loc':(66.875461, -162.568523)}
siteOutput['Sea of Okhotsk']={'loc':(59.575740, 147.571006)}
siteOutput['North Sea of Okhotsk']={'loc':(60.849945, 159.809775)}
siteOutput['Chukchi']={'loc':(70.327353, -160.855718)}
siteOutput['Drew Point']={'loc':(70.877725, -153.867593)}
siteOutput['Barrow']={'loc':(71.2956, -156.7664)}
siteOutput['Eastern Beaufort Coast']={'loc':(70.030979, -143.224450)}
siteOutput['Churchill']={'loc':(58.744589, -94.098072)}
siteOutput['Western Canadian Archipelago']={'loc':(73.228622, -124.288321)}
siteOutput['Central Canadian Archipelago']={'loc':(78.518153, -112.774650)}
siteOutput['Canadian Archipelago (Interior A)']={'loc':(73.960780, -110.401603)}
siteOutput['Canadian Archipelago (Interior B)']={'loc':(74.225746, -83.331292)}
siteOutput['Canadian Archipelago (Interior C)']={'loc':(77.766749, -100.492169)}
siteOutput['Northern Ellesmere']={'loc':(83.331922, -78.183598)}
siteOutput['North Greenland']={'loc':(83.372621, -33.958968)}
siteOutput['West Greenland (North)']={'loc':(75.026648, -57.293929)}
siteOutput['West Greenland (South)']={'loc':(64.977696, -52.064437)}
siteOutput['Baffin Bay']={'loc':(72.827891, -65.335921)}
siteOutput['East Greenland (North)']={'loc':(79.581088, -20.775375)}
siteOutput['East Greenland (South)']={'loc':(73.807454, -21.830063)}
siteOutput['Fram Strait']={'loc':(77.968355, 3.538632)}
siteOutput['East Svalbard']={'loc':(77.681053, 24.632383)}
siteOutput['West Svalbard']={'loc':(78.355968, 11.009335)}
siteOutput['North Svalbard']={'loc':(80.198051, 19.654310)}
siteOutput['Barents Sea']={'loc':(77.434891, 41.683165)}
siteOutput['Marre Sale']={'loc':(69.764959, 65.343310)}
siteOutput['Kara Sea']={'loc':(75.813956, 91.229132)}
siteOutput['Cape Bykovsky']={'loc':(71.900216, 130.012575)}
siteOutput['Laptev Sea']={'loc':(75.931939, 114.168586)}
siteOutput['Eastern Siberian Sea (West)']={'loc':(71.466926, 158.326560)}
siteOutput['Eastern Siberian Sea (East)']={'loc':(70.126904, 175.296515)}
siteOutput['Bering Strait']={'loc':(67.156341, -168.994759)}

# add deltas
siteOutput['Ob River Delta']={'loc':(66.75, 69.25)}
siteOutput['Lena River Delta']={'loc':(73.25, 127.25)}
siteOutput['Mackenzie River Delta']={'loc':(69.25, -134.75)}
siteOutput['Yukon River Delta']={'loc':(62.75, -164.75)}
siteOutput['Severnaya-Dvina River Delta']={'loc':(64.75, 40.25)}
siteOutput['Pechora River Delta']={'loc':(68.25, 54.25)}
siteOutput['Yana River Delta']={'loc':(71.25, 135.75)}
siteOutput['Hayes River Delta']={'loc':(56.43, -94.25)}
siteOutput['Back River Delta']={'loc':(67.25, -95.75)}
siteOutput['Winisk River Delta']={'loc':(55.25, -85.25)}
siteOutput['Onega River Delta']={'loc':(63.75, 38.25)}
siteOutput['Burnside River Delta']={'loc':(67.25, -108.25)}
siteOutput['Ellice River Delta']={'loc':(67.75, -104.25)}

# add open water points
siteOutput['Arctic Ocean A']={'loc':(84.891602, -61.583258)}
siteOutput['Arctic Ocean B']={'loc':(81.001454, 143.201901)}
siteOutput['Arctic Ocean C']={'loc':(76.128028, -149.298103)}
siteOutput['Arctic Ocean D']={'loc':(83.433205, -142.618416)}
siteOutput['Arctic Ocean E']={'loc':(84.118592, 56.190179)}
siteOutput['Arctic Ocean F']={'loc':(84.118592, 56.190179)}
siteOutput['Arctic Ocean G']={'loc':(84.118592, 56.190179)}
siteOutput['Arctic Ocean H']={'loc':(80.862938, 164.647214)}
siteOutput['Arctic Ocean I']={'loc':(79.644474, -167.579350)}
siteOutput['Arctic Ocean J']={'loc':(73.410644, -177.071537)}
# make sure all coast on coast

for k in siteOutput.keys():  
    sLat=siteOutput[k]['loc'][0]
    sLon=np.remainder(siteOutput[k]['loc'][1]+360., 360.)
    
    dists=h.haversine(sLon, sLat, lonSelect, latSelect)
    ind=np.where(dists==np.nanmin(dists))
    loc={}
    loc['nj']=ind[0][0]
    loc['ni']=ind[1][0]
    siteOutput[k]['inds']=loc 
    siteOutput[k]['cellLL']=(latSelect[loc['nj'],loc['ni']], lonSelect[loc['nj'],loc['ni']])
    
    dists=h.haversine(sLon, sLat, lonSelectOB, latSelectOB)
    ind=np.where(dists==np.nanmin(dists))
    loc={}
    loc['nj']=ind[0][0]
    loc['ni']=ind[1][0]
    siteOutput[k]['indsOB']=loc 
    siteOutput[k]['obsCellLL']=(latSelectOB[loc['nj'],loc['ni']], lonSelectOB[loc['nj'],loc['ni']])
pickle.dump(siteOutput, open("/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/siteOutput_justLocations.p", "wb" ) )                                    
    
    
print 'opening all netcdfs'
files=[]
for fn in dirListOrig:
    files.append(Dataset(fn, 'r'))
    print fn, (datetime.now()-startTime)
analysisFiles=[]    
for fn in dirList:
    analysisFiles.append(Dataset(fn, 'r'))
    print fn, (datetime.now()-startTime)
    
files1850=[]
for fn in dirListOrig1850:
    files1850.append(Dataset(fn, 'r'))
    print fn, (datetime.now()-startTime)
analysisFiles1850=[]    
for fn in dirList1850:
    analysisFiles1850.append(Dataset(fn, 'r'))
    print fn, (datetime.now()-startTime)
        
print  'finished opening netcdfs'
numYears1850=(2200-401)
numDays1850=numYears1850*365

numDays=(2100-1849)*365
numModels=len(dirListOrig)
numSites=len(siteOutput)
sic=np.nan*np.ones((len(siteOutput),numDays,numModels))
thresh=15.
it=0
day=np.arange(0,numDays)
time=day/365.+1850

numObsYears=fObsA.variables['time'].size

for k in siteOutput.keys():
    siteOutput[k]['sic']=np.nan*np.ones((numDays,numModels))     
    siteOutput[k]['nSIF']=np.nan*np.ones((251,numModels))           
    siteOutput[k]['nSIF1850']=[]
    siteOutput[k]['sic1850']=[]
    siteOutput[k]['nsifobs']=fObsA.variables['numSIF'][:,siteOutput[k]['indsOB']['nj'],siteOutput[k]['indsOB']['ni']]      

    # this is where to add the additional values that irina wants. 
    # irina wants (first day, last day), (full and continuous), (model and observation) 

    # also add a query to the regridded observations.
    

# analysis of all the 1850 runs
# this one starts in MY 402, so not exactly a century. need a running total...
runningInd=0                           
for fitt in range(len(analysisFiles1850)):
    f=analysisFiles1850[fitt]
    numYears=len(f.variables['time'][:])
    for k in siteOutput.keys():
        siteOutput[k]['nSIF1850'].extend(f.variables['numSIF'][:,siteOutput[k]['inds']['nj'],siteOutput[k]['inds']['ni']])
        print str(runningInd), str(len(siteOutput[k]['nSIF1850'])), 'NumSIF1850', k, str(numYears),(datetime.now()-startTime)
    runningInd=runningInd+1

runningInd=0                                     
for fitt in range(len(files1850)):
    f=files1850[fitt]
    numTime=len(f.dimensions['time'])
    numYears=numTime/365   
    #by years, not days
    for t in range(numYears):            
        startInd=t*365
        stopInd= (t+1)*365
        siSlice=f.variables['aice_d'][startInd:stopInd,:,:]
        
        for k in siteOutput.keys():
            siteOutput[k]['sic1850'].extend(siSlice[:,siteOutput[k]['inds']['nj'],siteOutput[k]['inds']['ni']])                                                
            print str(runningInd), str(t),str(len(siteOutput[k]['sic1850'])), 'daily sic 1850',k,(datetime.now()-startTime) 
            
    runningInd=runningInd+1
    
for k in siteOutput.keys():
    for sk in ['nSIF1850','sic1850']:
        temp=np.asarray(siteOutput[k][sk])
        siteOutput[k][sk]=temp
                                                         
# 1850-present model output                                                                                                                                                                                                                       numDays=(2100-1849)*365
startTime=datetime.now()                                              
for fitt in range(len(files)):
    f=files[fitt]
    numTime=len(f.dimensions['time'])
    numYears=numTime/365   
    #by years, not days
    for t in range(numYears):            
        startInd=t*365
        stopInd= (t+1)*365
        siSlice=f.variables['aice_d'][startInd:stopInd,:,:]       
        startInd2=(numDays-numTime)+(t*365)
        stopInd2=(numDays-numTime)+((t+1)*365)
        
        for k in siteOutput.keys():
            siteOutput[k]['sic'][startInd2:stopInd2,fitt]=siSlice[:,siteOutput[k]['inds']['nj'],siteOutput[k]['inds']['ni']]                                                
            print str(fitt), 'daily sic',k,(datetime.now()-startTime)       
                           

for fitt in range(len(analysisFiles)):
    f=analysisFiles[fitt]
    year=f.variables['time'][:]
    for k in siteOutput.keys():
        siteOutput[k]['nSIF'][-len(year):,fitt]=f.variables['numSIF'][:,siteOutput[k]['inds']['nj'],siteOutput[k]['inds']['ni']]
        print str(fitt), 'NumSIF', k,(datetime.now()-startTime)  
                                                                                                                  
siteOutput['obsYears']=np.arange(1979, 2014)                                                                                                                                                                                                                                           
pickle.dump(siteOutput, open("/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/siteOutput.p", "wb" ) )                                    
                         