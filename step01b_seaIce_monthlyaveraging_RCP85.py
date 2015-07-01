# code written by K. Barnhart in 2014 and 2015 in support of
# Barnhart et al., TITLE, YEAR
# 
# this code represents the first step in analysing the CESMLE sea ice output
# it creates monthly averages of the sea ice concentration and sea ice extent.
# 
# to run it, please verify that all the modules listed below are installed on 
# your machine and change the input and output file paths.  

print 'importing modules'
import nio
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import os
import glob

# set paths for input and output folder
# input folder should contain all of the aice_d files from the CESM-LE
# for users of Yellowstone/glade, this file is located at 
# glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/
  
pathOut=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/'
pathIn=u'/Volumes/Pitcairn/seaicePPF/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/'

# find all files (nh) limits us to just the northern hemisphere
dirList=glob.glob(pathIn+'b*h*.nc')

# hard code in the names for the run parts
bgKey=u'B1850C5CN'
runPart1key=u'B20TRC5CNBDRD'
runParts23key=u'BRCP85C5CNBDRD'
runParts23keyRCP45=u'BRCP45C5CNBDRD'


owThresh=15 # sea ice concentration for "open Water"

daysPerMonth=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
cumDaysInYear=np.cumsum(daysPerMonth)

# select unique runs
uniqueRuns=[]
for fn in dirList:
    pathSplit=fn.split('/')
    fnSplit=pathSplit[-1].split('.')
    if fnSplit[2]==bgKey: # if the control run
        uniqueRuns.append(fn)
    if fnSplit[2]==runPart1key:
        uniqueRuns.append(fn)
        
# make sure monthly average do not already exist 
fileNames=[]
for fn in uniqueRuns:
    pathSplit=fn.split('/')
    fnSplit=pathSplit[-1].split('.')
    if fnSplit[2]==bgKey: # if the control run
        fnTemp=str.join(".", fnSplit[:-1])+'.timeseries.nc'
    if fnSplit[2]==runPart1key:
        fnSplit[2]=runPart1key+'-'+runParts23key   
        fnTemp=str.join(".", fnSplit[:-2])+'.timeseries.nc'
    
    fn_monthOut=pathOut+fnTemp[:-2]+'monthlyAvg.nc'
     
    if os.path.isfile(fn_monthOut)==False:
        fileNames.append(fn)
        print fn_monthOut 

# process the remaining input files. 
startTime = datetime.now()
for fn in fileNames: 
	try:
		print 'Starting averaging of', fn
		pathSplit=fn.split('/')
		fnSplit=pathSplit[-1].split('.')
	 
		# open files, get aice_d, and time and put into variables with the same name 
		if fnSplit[2]==bgKey: # if the control run
			fnTemp=str.join(".", fnSplit[:-1])+'.timeseries.nc'
			f=nio.open_file(fn, 'r')
			aice_d=f.variables['aice_d'][:,:,:]
			time=f.variables['time'][:]
		if fnSplit[2]==runPart1key:
			fnTemp=str.join(".", fnSplit[:-2])+'.timeseries.nc'
		
			# now instead you need to open THREE FILES:
			# first open the first one
			f=nio.open_file(fn, 'r')
			aice_d=f.variables['aice_d'][:,:,:]
			time=f.variables['time'][:]
	   
			# now find the other two and add them to the end 
			pathConstruct2=str.join('/',pathSplit[:-1])+'/*'+runParts23key+'*'+fnSplit[4]+'*'+fnSplit[7]+'*nc'
			findFiles=np.sort(glob.glob(pathConstruct2)) 
			for fother in findFiles:
			
				f1=nio.open_file(fother, 'r')
				time=np.append(time,f1.variables['time'][:], 0)
				aice_d=np.append(aice_d, f1.variables['aice_d'][:,:,:], 0)
				del f1
			fnSplit[2]=runPart1key+'-'+runParts23key   
			fnTemp=str.join(".", fnSplit[:-2])+'.timeseries.nc'
		fn_monthOut=pathOut+fnTemp[:-2]+'monthlyAvg.nc'  
		print fn_monthOut
		# calculate the number of years and get the x-y dimension
		numyears=time.size/365 
		ni=f.dimensions['ni']
		nj=f.dimensions['nj'] 
		area=f.variables['tarea'][:,:]
		nummonths=numyears*12
	
	
		# for each cell do the following create monthly averages of SIC and sea ice area   
	
		# initialize output
		monthlyTimestep=np.zeros(numyears*12)
		monthlyTimeBounds=np.zeros((numyears*12,2))
		monthlyAvgSIC=np.zeros((numyears*12, nj, ni))
		monthlyAvgSIA=np.zeros((numyears*12, 1))
		monthlyAvgSIA_thresh=np.zeros((numyears*12, 1))
	  
		monthlyTimestep=monthlyTimestep.astype(np.float32)
		monthlyTimeBounds=monthlyTimeBounds.astype(np.float32)
		monthlyAvgSIC=monthlyAvgSIC.astype(np.float32)
		monthlyAvgSIA=monthlyAvgSIA.astype(np.float32)
		monthlyAvgSIA_thresh=monthlyAvgSIA_thresh.astype(np.float32)
	
		# loop through each month. 
		arrayItter=0
		for i in range(numyears):
			#print 'year = ',i+int(t[0]/365)
			for j in range(12):
				indStart=int(i*365+np.remainder(cumDaysInYear[j-1], 365))
				indStop=int(i*365+cumDaysInYear[j]-1)
			
				# select the timestamps and ice concentrations for the current month
				selt=time[indStart:indStop]
				selIce=aice_d[indStart:indStop,:,:]
			
				# calculate the mean ice concentration
				meanIce=selIce.mean(axis=0)
				meanIce[selIce[0,:,:]==f.variables['aice_d'].__dict__['_FillValue']]=f.variables['aice_d'].__dict__['_FillValue']
				monthlyAvgSIC[arrayItter,:,:]=meanIce
			
				# calculate the mean ice extent
				sia=[]
				sia_thresh=[]
				conversion=1000000.*1000000.
				for jj in range(selIce.shape[0]):
				
					selice2=selIce[jj,:,:]
					selice2[selice2>1e29]=0
					temparea=area
					temparea[temparea>1e29]=0
					temparea=temparea/conversion
					temparea[selice2==0]=0
				
					sia.append((np.sum((selice2/100.)*(temparea))))
					sia_thresh.append((np.sum((selice2>owThresh)*(temparea))))

				# save values into output structures. 

				monthlyAvgSIA[arrayItter] = np.mean(sia) # 100% is represted ast 100 instead of 1
				monthlyAvgSIA_thresh[arrayItter] = np.mean(sia_thresh)
			
				monthlyTimestep[arrayItter]=round(selt.mean(axis=0))
				monthlyTimeBounds[arrayItter, 0]=selt.min(axis=0)
				monthlyTimeBounds[arrayItter, 1]=selt.max(axis=0)
				arrayItter+=1
		del time
		del aice_d  
		  
		## Create this monthly averaged file as a new netcdf
		fMonth=Dataset(fn_monthOut, 'w',format='NETCDF4')
	
		# create all the dimentions, set time to unlimited 
		for k in f.dimensions.keys():
			if f.unlimited(k)==True:
				fMonth.createDimension(k, monthlyTimestep.size)#None)
			else:
				fMonth.createDimension(k, f.dimensions[k])
				print k, f.dimensions[k]
			
		# use the netCDF4 instead of pyNIO since it seems to work much better with unlimited variables      
		fMonthVars={}
		for key in {'TLAT', 'TLON','latt_bounds','lont_bounds','time_bounds', 'time'}:
			#print 'creating ', key
			# the netCDF4 module requires that if a fill value exists, it must be declared when the variable is created.
			try:
				fMonthVars[key]=fMonth.createVariable(key, f.variables[key].typecode(), f.variables[key].dimensions, fill_value=f.variables[key].__dict__['_FillValue'])
			except:
				fMonthVars[key]=fMonth.createVariable(key, f.variables[key].typecode(), f.variables[key].dimensions,fill_value=f.variables['aice_d'].__dict__['_FillValue'])
			# sett all the attribute keys.
			atts = f.variables[key].__dict__
			for attKey in atts.keys():
				if attKey != '_FillValue':
					setattr(fMonth.variables[key],attKey,atts[attKey])         
				
		# create the montly averaged sea ice variable 
		monthAvgKey='aice_d_monthAvg'
		fMonthVars[monthAvgKey]=fMonth.createVariable(monthAvgKey, f.variables['aice_d'].typecode(), f.variables['aice_d'].dimensions,fill_value=f.variables['aice_d'].__dict__['_FillValue'])
		#print 'creating aice_d_monthAvg'
		atts = f.variables['aice_d'].__dict__
		for attKey in atts.keys():
			if attKey is not '_FillValue':
				setattr(fMonth.variables[monthAvgKey],attKey,atts[attKey])
		setattr(fMonth.variables[monthAvgKey],'long_name','mean monthly sea ice concentration  (aggregate)')   
	
	
		monthAvgKey='monthlyAvgSIA'
		fMonthVars[monthAvgKey]=fMonth.createVariable(monthAvgKey, f.variables['aice_d'].typecode(), 'time',fill_value=f.variables['aice_d'].__dict__['_FillValue'])
		#print 'creating aice_d_monthAvg'
		atts = f.variables['aice_d'].__dict__
		for attKey in atts.keys():
			if attKey is not '_FillValue':
				setattr(fMonth.variables[monthAvgKey],attKey,atts[attKey])
		setattr(fMonth.variables[monthAvgKey],'long_name','mean monthly ice area  (aggregate)') 
		setattr(fMonth.variables[monthAvgKey],'units','million square kilometers')   
	
		monthAvgKey='monthlyAvgSIA_thresh'
		fMonthVars[monthAvgKey]=fMonth.createVariable(monthAvgKey, f.variables['aice_d'].typecode(), 'time',fill_value=f.variables['aice_d'].__dict__['_FillValue'])
		#print 'creating aice_d_monthAvg'
		atts = f.variables['aice_d'].__dict__
		for attKey in atts.keys():
			if attKey is not '_FillValue':
				setattr(fMonth.variables[monthAvgKey],attKey,atts[attKey])
		setattr(fMonth.variables[monthAvgKey],'long_name','mean monthly ice area (using 15% threshold) (aggregate)')   
		setattr(fMonth.variables[monthAvgKey],'units','million square kilometers')
	
		# put data into variables, first the ones we are copying over. 
		#print 'putting data into standard variables'
		for key in {'TLAT', 'TLON','latt_bounds','lont_bounds'}:
			fMonthVars[key][:,:]=f.variables[key][:]
		# now, the ones we have created (those with time as a dimention)
		fMonthVars['time'][:]=monthlyTimestep
		fMonthVars['time_bounds'][:,:]=monthlyTimeBounds
		fMonthVars['aice_d_monthAvg'][:,:,:]=monthlyAvgSIC
	
		fMonthVars['monthlyAvgSIA'][:]=monthlyAvgSIA 
		fMonthVars['monthlyAvgSIA_thresh'][:]=monthlyAvgSIA_thresh

		# close and delete the output netCDF variable, retain f, so we can do the next part of the analysis.
		fMonth.close() 
		del fMonth

		print'finished averaging of ',fn , (datetime.now()-startTime)
	except:
		print 'oops... ', fn,  'didnt run' , (datetime.now()-startTime)