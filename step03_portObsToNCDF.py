#put observed sea ice into netcdf!

#! /usr/bin/python

#from grass import script as gs
import os
import glob
import numpy as np
from netCDF4 import Dataset
import nio
from datetime import datetime
from datetime import timedelta
import matplotlib.pylab as plt
plt.close('all')
import cProfile, pstats, StringIO
pr = cProfile.Profile()
pr.enable()

# need to fix
# 1. time stamping 

folderpath=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/nsidcObservations/'
pathIn=u'/Volumes/Pitcairn/seaicePPF/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/'
fn_sat='/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/satelliteObs.timeseries.nc'

path=folderpath+'/supportingFiles'
os.chdir(path)

f = open('coast_25n.msk', 'r')
coastgrid = np.fromfile(f,dtype=np.uint8)
f.close()
coastgrid=coastgrid.reshape((448,304))

f = open('gsfc_25n.msk', 'r')
landmask = np.fromfile(f, dtype=np.uint8)
f.close()
landmask=landmask.reshape((448,304))

f = open('psn25lats_v2.dat', 'r')
latgrid = np.fromfile(f, dtype=np.int32)/100000.
f.close()
latgrid=latgrid.reshape((448,304))

f = open('psn25lons_v2.dat', 'r')
longrid = np.fromfile(f, dtype=np.int32)/100000.
f.close()
longrid=longrid.reshape((448,304))

#plt.figure()
#plt.imshow(longrid, vmin=-180, vmax=180)
#plt.show()
#
#plt.figure()
#plt.imshow(latgrid, vmin=0, vmax=90)
#plt.show()


yrs=range(1979, 2015)
# time is in "days since 1979-01-01 00:00:00"  
reftime=datetime(min(yrs), 1, 1, 0, 0, 0)
endtime=datetime(max(yrs), 12, 31, 0, 0, 0)

totalTime=endtime-reftime
numFiles=12784+365


#totalTime.days

timestep=np.zeros(numFiles)    
sicAll=np.zeros((numFiles, 448,304))
timeBounds=np.zeros((numFiles,2))
itters=[]
mD=[]
i=0    
startTime = datetime.now()  
dates=[]

#path=u'/Users/katherinebarnhart/Desktop/RESEARCH/seaIcePPF/b*timeseries.nc'
dirList=glob.glob(pathIn+'b*nh*.nc')
f=nio.open_file(dirList[0], 'r')

fillVal=f.variables['aice_d'].__dict__['_FillValue']
f.close()
del f

for y in yrs:
    path=folderpath+str(y)
    os.chdir(path)
    filenames = sorted(glob.glob('nt*.bin'))
    if (len(filenames)>367) or (len(filenames)==263):
        datesSort=[fn[3:11] for fn in filenames]
        udates=sorted(set(datesSort))
        newFiles=[] 
        for ud in udates:
	    tempFiles = sorted(glob.glob('nt*'+ud+'*.bin'))
	    if len(tempFiles)==1:
	        newFiles.append(tempFiles[0])
	    else:   
       	        nums=[int(fn[13:15]) for fn in tempFiles]
       	        moreFiles=sorted(glob.glob('nt*'+ud+'*'+str(max(nums))+'*.bin'))
       	        newFiles.append(moreFiles[0])
        filenames=newFiles
        
    print len(filenames)  

	    
	   # multiple versions
    for filename in filenames:
	    
    	   f = open(filename, 'r')
	   
	   data = np.fromfile(f, dtype=np.uint8)
	   f.close()
	   
	   year=str(y)
	   month_num = int(filename[7:9])
	   day=int(filename[9:11])

	   
	   
	   if (i==0):
	       print 'firstday'
	       date=datetime(y, month_num, day-1, 12)
	       dates.append(date)
	       sic=data[300:]/250.*100.
	       sic=sic.reshape((448,304))
	       
	       landMask=sic>100
	       sic[landMask]=fillVal
	       
	       sicAll[i,:,:]=sic
	       timestep[i]=0.5
	       timeBounds[i, 0]=timestep[i]-0.5
	       timeBounds[i, 1]=timestep[i]+0.5
	       itters.append(i)
	       i+=1
	       
              
	   if (i>0):
	       date=datetime(y, month_num, day, 12)
	       dates.append(date)
	       
	       localdt=dates[-1]-dates[-2]
	       missingDays=localdt.total_seconds()/(3600*24.)
	       
	       if int(missingDays)==2:
	           print i, 'filling', str(missingDays)
	           sicAll[i,:,:]=sicAll[i-1,:,:]  
	           timestep[i]=timestep[i-1]+1.
	           timeBounds[i, 0]=timestep[i]-0.5
	           timeBounds[i, 1]=timestep[i]+0.5
	           itters.append(i)
	           i+=1
	           mD.append(missingDays)
	           
    	       if missingDays>20:
	           print i, 'filling', str(missingDays)
	           missing=missingDays
    	           while missing>1:
        	           sicAll[i,:,:]=sicAll[i-1,:,:]  
        	           timestep[i]=timestep[i-1]+1.
        	           timeBounds[i, 0]=timestep[i]-0.5
        	           timeBounds[i, 1]=timestep[i]+0.5
        	           itters.append(i)
        	           i+=1
        	           mD.append(missingDays)
        	           missing=missing-1
        	           date=datetime(y, month_num, day, 12)
        	           dates.append(date)  
	                   
	   date=datetime(y, month_num, day, 12)
	   dates.append(date)        
 	   dt=date-reftime  					
	   sic=data[300:]/250.*100.
	   sic=sic.reshape((448,304))
	   landMask=sic>100
	   sic[landMask]=fillVal
	   sicAll[i,:,:]=sic
	   timestep[i]=dt.total_seconds()/(3600*24.)
	   timeBounds[i, 0]=timestep[i]-0.5
	   timeBounds[i, 1]=timestep[i]+0.5
	   print i, timestep[i] ,filename, (datetime.now()-startTime)
	   itters.append(i)
           i+=1 

 # remove extra
       
              
                            
# 2. get this into an netcdf:

dirList=glob.glob(pathIn+'*nc')
f=nio.open_file(dirList[0], 'r')

#fn_sat='/Users/katherinebarnhart/Desktop/RESEARCH/seaIcePPF/satelliteObs.timeseries.nc'
fsat=Dataset(fn_sat, 'w',format='NETCDF4')
    
# create all the dimentions, set time to unlimited 

fsat.createDimension('time', None)
fsat.createDimension('ni', 304)
fsat.createDimension('nj', 448)
fsat.createDimension('d2', 2)

fsatVars={}
for key in {'TLAT', 'TLON','time', 'time_bounds'}:
    #print 'creating ', key
    # the netCDF4 module requires that if a fill value exists, it must be declared when the variable is created.
    try:
        fsatVars[key]=fsat.createVariable(key, f.variables[key].typecode(), f.variables[key].dimensions, fill_value=f.variables[key].__dict__['_FillValue'])
    except:
        fsatVars[key]=fsat.createVariable(key, f.variables[key].typecode(), f.variables[key].dimensions)
    # sett all the attribute keys.
    atts = f.variables[key].__dict__
    for attKey in atts.keys():
        if attKey != '_FillValue':
            setattr(fsat.variables[key],attKey,atts[attKey]) 
              
setattr(fsat.variables['time'],'calendar','gregorian')
setattr(fsat.variables['time'], 'units', 'days since '+str(yrs[0])+'-01-01 00:00:00')
setattr(fsat.variables['time_bounds'],'calendar','gregorian')
setattr(fsat.variables['time_bounds'], 'units', 'days since '+str(yrs[0])+'-01-01 00:00:00')
monthAvgKey='satelliteSIC'
fsatVars[monthAvgKey]=fsat.createVariable(monthAvgKey, f.variables['aice_d'].typecode(), f.variables['aice_d'].dimensions,fill_value=f.variables['aice_d'].__dict__['_FillValue'])
#print 'creating aice_d_monthAvg'
atts = f.variables['aice_d'].__dict__
for attKey in atts.keys():
    if attKey is not '_FillValue':
        setattr(fsat.variables[monthAvgKey],attKey,atts[attKey])
setattr(fsat.variables[monthAvgKey],'long_name','Sea Ice Concentration (satellite)')   
    
fsatVars['time'][:]=timestep
fsatVars['time_bounds'][:,:]=timeBounds
fsatVars['TLAT'][:,:]=latgrid
fsatVars['TLON'][:,:]=longrid
fsatVars[monthAvgKey][:,:,:]=sicAll

fsat.close()
(datetime.now()-startTime)


#pr.disable()
#s = StringIO.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
#print s.getvalue()	
	