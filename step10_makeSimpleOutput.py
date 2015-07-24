print 'importing modules'
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import glob
import nio

startTime = datetime.now()

path=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/'


bgKey=u'B1850C5CN'
runPart1key=u'B20TRC5CNBDRD'
runParts23key=u'BRCP85C5CNBDRD'
runParts23keyRCP45=u'BRCP45C5CNBDRD'

nskey=['nh','sh']
rcpName=['RCP85', 'RCP45']
rcpkey=[runPart1key+'-'+runParts23key,runPart1key+'-'+runParts23keyRCP45]
# two loops. northern/southern + rcp 8.5/4.5
for nsk in nskey:
    dirList2=glob.glob(path+'*'+bgKey+'*'+nsk+'*Analysis.nc')
    analysisFiles2=[]    
    for fn in dirList2:
        analysisFiles2.append(Dataset(fn, 'r'))
    #print fn, (datetime.now()-startTime)    

    numSIF1850=[]
    fir1850=[]
    las1850=[]
    time1850=[]
    print 'opening background: ', nsk
    for fn in dirList2:
        f=Dataset(fn, 'r')
        numSIF1850.extend(f.variables['numSIF'][:,:,:])
        fir1850.extend(f.variables['firstDOY'][:,:,:])
        las1850.extend(f.variables['lastDOY'][:,:,:])
        time1850.extend(f.variables['time'][:])
        print fn, (datetime.now()-startTime) 
        f.close()   

    nSIF1850=np.array(numSIF1850)
    first1850=np.array(fir1850)
    last1850=np.array(las1850)
    del numSIF1850
    del fir1850
    del las1850

    for ittR in range(len(rcpkey)):
        fn_monthOut=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/justNSIF_ensembleAndBG.'+nsk+'.'+rcpName[ittR]+'.nc'
        print fn_monthOut

        dirList=glob.glob(path+'*'+rcpName[ittR]+'*'+nsk+'*Analysis.nc')

        analysisFiles=[]    
        for fn in dirList:
            analysisFiles.append(Dataset(fn, 'r'))

        key='numSIF'
        f=nio.open_file(dirList[0],'r')
        landMask=f.variables[key][0,:,:]>1200
        fillVal=f.variables[key].__dict__['_FillValue']
        ni=f.dimensions['ni']
        nj=f.dimensions['nj']
        nm=len(dirList)
        numModels=len(dirList)
        lastYear=int(np.floor(f.variables['time'][:].max()/365))
        numYears=lastYear-1850

        if numYears==249:
            numYears=250
            lastYear=lastYear-1

        f.close()
        del f

        numSIF=np.nan*np.ones((numModels,numYears,nj,ni))
        first=np.nan*np.ones((numModels,numYears,nj,ni))
        last=np.nan*np.ones((numModels,numYears,nj,ni))

        i=0

        fn=dirList[0]
        f=Dataset(fn, 'r')
        ensemble_time=f.variables['time'][:]/365
        f.close()
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

        fillVal=999

        nSIF1850[nSIF1850>400]=fillVal
        nSIF[np.isnan(nSIF)]=fillVal
        nSIF[nSIF>400]=fillVal


        print  'finished getting num SIF from netcdfs'
        f=nio.open_file(dirList[0], 'r')
        numyears=nSIF.shape[1]
        ni=f.dimensions['ni']
        nj=f.dimensions['nj']

        nm=nSIF.shape[0]

        fn_monthOut=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/justNSIF_ensembleAndBG.'+nsk+'.'+rcpName[ittR]+'.nc'
        fMonth=Dataset(fn_monthOut, 'w',format='NETCDF3_64BIT')

        fMonth.createDimension('nj', nj)
        fMonth.createDimension('ni', ni)
        fMonth.createDimension('time_ensemble', numyears)
        fMonth.createDimension('time_background1850', nSIF1850.shape[0])
        fMonth.createDimension('nvertices', 4)
        fMonth.createDimension('d2', 2)

        fMonth.createDimension('nm', nm)

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
        # create the monthly averaged sea ice variable 

        cKey='nm'
        fMonthVars[cKey]=fMonth.createVariable(cKey, 'i', ('nm') ,fill_value=fillVal)
        setattr(fMonth.variables[cKey],'long_name','Model Number') 
        setattr(fMonth.variables[cKey],'units','')
        fMonthVars[cKey][:]=np.arange(nm)


        cKey='nSIF_ensemble'
        fMonthVars[cKey]=fMonth.createVariable(cKey, 'i2', ('nm', 'time_ensemble', 'nj', 'ni') ,fill_value=fillVal)
        setattr(fMonth.variables[cKey],'long_name','Number of Sea Ice Free Days (CESM-LE ensemble)') 
        setattr(fMonth.variables[cKey],'units','days')
        setattr(fMonth.variables[cKey], 'coordinates', 'nm time_ensemble TLON TLAT')
        fMonthVars[cKey][:,:,:,:]=nSIF.astype(int)

        cKey='nSIF_background1850'
        fMonthVars[cKey]=fMonth.createVariable(cKey, 'i2', ('time_background1850', 'nj', 'ni') ,fill_value=fillVal)
        setattr(fMonth.variables[cKey],'long_name','Number of Sea Ice Free Days (CESM-LE background 1850)') 
        setattr(fMonth.variables[cKey],'units','days') 
        setattr(fMonth.variables[cKey], 'coordinates', 'time_background1850 TLON TLAT')

        fMonthVars[cKey][:,:,:]=nSIF1850.astype(int)

        cKey='time_ensemble'
        fMonthVars[cKey]=fMonth.createVariable(cKey, 'i', ('time_ensemble'))
        setattr(fMonth.variables[cKey],'long_name','Time (CESM-LE ensemble)') 
        setattr(fMonth.variables[cKey],'units','year')   
        fMonthVars[cKey][:]=ensemble_time

        cKey='time_background'
        fMonthVars[cKey]=fMonth.createVariable(cKey, 'i', ('time_background1850'))
        setattr(fMonth.variables[cKey],'long_name','Time (CESM-LE background)') 
        setattr(fMonth.variables[cKey],'units','year')   
        fMonthVars[cKey][:]=time1850


        fMonth.close()
        print 'done'
