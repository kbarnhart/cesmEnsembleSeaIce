print 'importing modules'
import nio
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import os
import glob
import pickle
import matplotlib.pylab as plt


pathIn=u'/Volumes/Pitcairn/seaicePPF/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/'
pathIn45=u'/Volumes/Pitcairn/seaicePPF/p/cesm0005/CESM-CAM5-BGC-ME/ice/proc/tseries/daily/aice_d/'

pathOut=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/'
pickleFile=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/code/running.p'

# set open water thresholds and day of year for
owThresh=15 # sea ice concentration for "open Water"
midSIF_nh=250  # day of year for "middle of SIF season" 215~= August 1st
midIce_nh=70

midSIF_sh=70
midIce_sh=250

# hard code in the names for the run parts
bgKey=u'B1850C5CN'
runPart1key=u'B20TRC5CNBDRD'

runParts23keyRCP85=u'BRCP85C5CNBDRD'
runParts23keyRCP45=u'BRCP45C5CNBDRD'

# get directory info
dirList=glob.glob(pathOut+'b*timeseries.monthlyAvg.nc') # this will give us each monthlyAverage file"

# identify files and make sure they are not running already. 
running=[]
if os.path.isfile(pickleFile)==True: 
    running=pickle.load(open(pickleFile, 'rb'))
    
fileNames=[]
for fn in dirList:
    pathSplit=fn.split('/')
    fnSplit=pathSplit[-1].split('.')
    fn_analysisOut=fn[:-24]+'Analysis.nc'
    if os.path.isfile(fn_analysisOut)==False:
        if (fn not in running):
            fileNames.append(fn)
            print fn_analysisOut


daysPerMonth=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
cumDaysInYear=np.cumsum(daysPerMonth)
startTime = datetime.now()
for fn in fileNames:
    if os.path.isfile(pickleFile)==True: 
        running=pickle.load(open(pickleFile, 'rb'))
    
    if (fn in running):
        print fn+' is already running'
    else:
        running.append(fn)
        pickle.dump(running, open(pickleFile, 'wb'))    
    
        print 'Starting Analysis   of', fn
        fn_analysisOut=fn[:-24]+'Analysis.nc'
            
        pathSplit=fn.split('/')
        fnSplit=pathSplit[-1].split('.')
        
        # open files, get aice_d, and time and put into variables with the same name 
        if fnSplit[2]==bgKey: # if the control run
                        
            # construct original file. 
            fn_inFind=pathIn+'*'+fnSplit[2]+'*'+fnSplit[4]+'*'+fnSplit[7]+'*'+fnSplit[-4]+'*nc'
            fn_in=glob.glob(fn_inFind)[0]
            
            f=nio.open_file(fn_in, 'r')
            aice_d=f.variables['aice_d'][:,:,:]
            time=f.variables['time'][:]
            
            numyears=time.size/365
            
            next_yr=str(int(fnSplit[-4][-8:-4])+1).zfill(4)
            
            fn_in_next_Find=pathIn+'*'+fnSplit[2]+'*'+fnSplit[4]+'*'+fnSplit[7]+'*'+next_yr+'*nc'
            
            nextFN=glob.glob(fn_in_next_Find)
            
            if len(nextFN)>0:
                fnext=nio.open_file(nextFN[0], 'r')
                aice_d=np.append(aice_d, fnext.variables['aice_d'][:365,:,:], 0)
                time=np.append(time, fnext.variables['time'][:365],0)
                numAnalyisYears=numyears # we can only do this until the second to last year becaue
            else:
                numAnalyisYears=numyears-1 # we can only do this until the second to last year 
        
        # now instead you need to open THREE FILES:
        # also check if the corresponding RCP other file has run, so as not to re-analysse 1850/1920-2006
        # first open the first one
        if fnSplit[2][:13]==runPart1key:
			pathAnalysisSplit=fn.split('/')
			fnSplitAnalysis=pathSplit[-1].split('.')

			fn0_find=pathIn+'*'+runPart1key+'*'+fnSplit[4]+'*'+fnSplit[7]+'*nc'

			fn0=glob.glob(fn0_find)[0]
			print 'opening the first file: ', fn0

			f=nio.open_file(fn0, 'r')
			aice_d=f.variables['aice_d'][:,:,:]
			time=f.variables['time'][:]
	   		
	   		numOrigYears=time.size/365.
	   		
			## now find the other two and add them to the end 
			## large vs. medium ensemble if statement
			if fnSplit[2][18:20]=="85": # rcp 8.5
				pathConstruct2=pathIn+'*'+runParts23keyRCP85+'*'+fnSplit[4]+'*'+fnSplit[7]+'*nc'
				fnSplitAnalysis[2]=runPart1key+'-'+runParts23keyRCP45
				
			else: # rcp 4.5
				pathConstruct2= pathIn45 +'*'+runParts23keyRCP45+'*'+fnSplit[4]+'*'+fnSplit[7]+'*nc'
				fnSplitAnalysis[2]=runPart1key+'-'+runParts23keyRCP85
			findFiles=np.sort(glob.glob(pathConstruct2)) 
			for fother in findFiles:
				print 'opening ', fother
				f1=nio.open_file(fother, 'r')
				time=np.append(time,f1.variables['time'][:], 0)
				aice_d=np.append(aice_d, f1.variables['aice_d'][:,:,:], 0)
				f1.close()
				del f1
			numyears=time.size/365
			numAnalyisYears=numyears-1 # we can only do this until the second to last year

			# check for run with RCP other
			pathAnalysisSplit[-1]=str.join('.', fnSplitAnalysis)
			otherAnalysisPath=str.join('/', pathAnalysisSplit)
			otherAnalysisPath=otherAnalysisPath[:-24]+'Analysis.nc'
			if os.path.isfile(otherAnalysisPath)==True:
				otherFile=True
				
				startyear=int(numOrigYears)
			else:
				startyear=0
				

        # nh vs sh toggle        
        if fnSplit[-5][-2:] == 'nh':
            midSIF=midSIF_nh
            midIce=midIce_nh
        else:
            midSIF=midSIF_sh
            midIce=midIce_sh
            

        # calculate the number of years and get the x-y dimension
        
        ni=f.dimensions['ni']
        nj=f.dimensions['nj'] 
    
        print'analysis', (datetime.now()-startTime)
        # 2. First day, last day, and duration of open water
        # Save these things to a new netcdf 
        print 'opening the montly averages'
        fMonth=Dataset(fn, 'r',format='NETCDF4')
        monthlyAvgSIC=fMonth.variables['aice_d_monthAvg'][:,:,:]
        fMonth.close()
        
        iceMax=monthlyAvgSIC.max(axis=0)
        fillVal=f.variables['aice_d'].__dict__['_FillValue']
        landMask=f.variables['aice_d'][0,:,:]==fillVal
        
        print 'done opening files' 

        # of delayed freeze up. 
        iceVal=1000.
        firstDOY=iceVal*np.ones((numAnalyisYears, nj, ni))
        lastDOY=iceVal*np.ones((numAnalyisYears, nj, ni)) # oops. fix this.
        firstDOYcontinous=iceVal*np.ones((numAnalyisYears, nj, ni))
        lastDOYcontinous=iceVal*np.ones((numAnalyisYears, nj, ni))
        durationContinous=iceVal*np.ones((numAnalyisYears, nj, ni))
        numSIF=365.*np.ones((numAnalyisYears, nj, ni))
        springMix=iceVal*np.ones((numAnalyisYears, nj, ni))
        fallMix=iceVal*np.ones((numAnalyisYears, nj, ni))
        
        yearlyTimestep=np.zeros(numAnalyisYears)
        yearlyTimeBounds=np.zeros((numAnalyisYears,2))
        
        # if a prior analysis file exists, use it here:
        if otherFile==True:
        	print 'getting output from previously run file'
        	fOtherRCP=nio.open_file(otherAnalysisPath)
        	firstDOY[:startyear,:,:]=fOtherRCP.variables['firstDOY'][:startyear,:,:] 
        	lastDOY[:startyear,:,:]=fOtherRCP.variables['firstDOY'][:startyear,:,:] 
        	firstDOYcontinous[:startyear,:,:]=fOtherRCP.variables['firstDOYcontinous'][:startyear,:,:] 
        	lastDOYcontinous[:startyear,:,:]=fOtherRCP.variables['lastDOYcontinous'][:startyear,:,:]
        	durationContinous[:startyear,:,:]=fOtherRCP.variables['durationContinous'][:startyear,:,:]
        	numSIF[:startyear,:,:]=fOtherRCP.variables['numSIF'][:startyear,:,:]
        	springMix[:startyear,:,:]=fOtherRCP.variables['springMix'][:startyear,:,:]
        	fallMix[:startyear,:,:]=fOtherRCP.variables['fallMix'][:startyear,:,:]
        	yearlyTimestep[:startyear]=fOtherRCP.variables['time'][:startyear]
        	yearlyTimeBounds[:startyear,:]=fOtherRCP.variables['time_bounds'][:startyear,:]
        	fOtherRCP.close()
        	del fOtherRCP
        
        iceAffected=(iceMax>0.01)*(iceMax<101)
        
        notIceAffected=iceAffected==False    
        
        land=iceMax>101
        #totalAffected=iceAffected.sum(axis=None)
        
        for i in range(startyear,numAnalyisYears):
            
            # make sligtly different if Southern hemisphere.
            
            indStart=int(i*365)+midIce  #do the analysis from march to march to fix the delayed freeze up problem. 
            indStop=int((i+1)*365)+midIce
            
                     
            selt=time[indStart:indStop]
            yearConc=aice_d[indStart:indStop,:,:]
            
            yearlyTimestep[i]=round(selt.mean(axis=0))
            yearlyTimeBounds[i, 0]=round(selt.min(axis=0))
            yearlyTimeBounds[i, 1]=round(selt.max(axis=0))
            
            percDone=round((i-startyear)*100.0/(numAnalyisYears-startyear))
            print 'year = ',i,'  : ', str(percDone), '%',(datetime.now()-startTime)
        
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
                         
                                                                                    
                        selTimeDOYs=np.asarray(range(1, 366))+midIce
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
                
                firstDOY[i,notIceAffected]=fillVal
                lastDOY[i, notIceAffected]=fillVal
                firstDOYcontinous[i, notIceAffected]=fillVal
                lastDOYcontinous[i,notIceAffected]=fillVal
                springMix[i, notIceAffected]=fillVal
                fallMix[i, notIceAffected]=fillVal
                durationContinous[i, notIceAffected]=fillVal
        
                        #else: # if there is no SIF season for this cell for this year
                            # place as a nan (for now to differentiate between nan and Fill Val
                            # and to make trend fitting and mpl plotting easier. 
                     # NOV 2014 change: we've set good defaults, so keep this as is
                     #       
                #            numSIF[i, nji, nii]=np.nan
                #            firstDOY[i, nji, nii]=np.nan
      		        #lastDOY[i, nji, nii]=np.nan
      		        #firstDOYcontinous[i, nji, nii]=np.nan
      		        #lastDOYcontinous[i, nji, nii]=np.nan
      		        #springMix[i, nji, nii]=np.nan
      		        #fallMix[i, nji, nii]=np.nan
      		        #durationContinous[i, nji, nii]=np.nan
        
        # 2.5 Save this to a netCDF
        
        fnAnalysis=fn[:-24]+'Analysis.nc'
        fAn=Dataset(fnAnalysis, 'w',format='NETCDF4')
        
        # create all the dimentions, set time to unlimited 
        for k in f.dimensions.keys():
            if f.unlimited(k)==True:
                fAn.createDimension(k, None)
            else:
                fAn.createDimension(k, f.dimensions[k])
                
        # use the netCDF4 instead of pyNIO since it seems to work much better with unlimited variables      
        fAnVars={}
        for key in {'TLAT', 'TLON','latt_bounds','lont_bounds','time_bounds', 'time'}:
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
        for key in {'TLAT', 'TLON','latt_bounds','lont_bounds'}:
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
            fAnVars[key]=fAn.createVariable(key, 'f', f.variables['aice_d'].dimensions, fill_value=stdAttributes['_FillValue'])
            
            # set all the attribute keys.
            for attKey in stdAttributes.keys():
                if attKey != '_FillValue':
                    setattr(fAn.variables[key],attKey,stdAttributes[attKey])
            setattr(fAn.variables[key],'long_name',longName[i])
            setattr(fAn.variables[key],'units',units[i])
            fAnVars[key][:,:,:]=data[i] 
        
        fAn.close()    
        f.close()
        del f
        del fAn
        print'created ',fnAnalysis , (datetime.now()-startTime)
