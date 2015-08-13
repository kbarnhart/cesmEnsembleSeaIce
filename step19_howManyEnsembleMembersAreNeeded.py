#! /usr/bin/python
# step 11 determine how many ensemble members you need to do the analysis correctly.

print 'importing modules'
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import glob
from scipy import stats
import math
import itertools as it
from scipy.stats import ks_2samp
import scipy  

import matplotlib.pylab as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif']=['Arial']
plt.close('all')

#import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['font.sans-serif']=['Arial']

startTime = datetime.now()

path=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/'
pathOut=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/'

path=u'/home/barnhark/seaIceEmergence/'
pathOut=u'/home/barnhark/seaIceEmergence/'

# hard code in the names for the run parts
bgKey=u'B1850C5CN'
runPart1key=u'B20TRC5CNBDRD'
runParts23keyRCP85=u'BRCP85C5CNBDRD'
runParts23keyRCP45=u'BRCP45C5CNBDRD'

nskey=['nh','sh']
rcpName=['RCP85', 'RCP45']
rcpkey=[runPart1key+'-'+runParts23keyRCP85,runPart1key+'-'+runParts23keyRCP45]
# two loops. northern/southern + rcp 8.5/4.5
for nsk in nskey:

   for ittR in range(len(rcpkey)):  
        print 'opening datasets'
        dataFN=path+'justNSIF_ensembleAndBG.'+nsk+'.'+rcpName[ittR]+'.nc'
        fnsif=Dataset(dataFN, 'r')
        
        nSIF85=np.array(fnsif.variables['nSIF_ensemble'][:,-180:,:,:],dtype=np.float16, order='C')
        nSIF85[nSIF85>1000]=np.nan
   
        CIfile=path+'justNSIF_ensembleAndBG.'+nsk+'.'+rcpName[ittR]+'.BootCI.nc'
        CI=Dataset(CIfile, 'r')
   
        CIsMeanUB=CI.variables['meanUB'][-180:,:,:]
        CIsMeanLB=CI.variables['meanLB'][-180:,:,:]
        CIsStdsUB=CI.variables['sdUB'][-180:,:,:]
        CIsStdsLB=CI.variables['sdLB'][-180:,:,:]
        
        fnAnalysis=pathOut+'numMembersNeeded.'+rcpName[ittR]+'.'+nsk+'.nc'
        print fnAnalysis
        
        numModels85=nSIF85.shape[0]
        numYears=(2100-1920)

        ni=320
        nj=104
        nm=numModels85


        import cProfile, pstats, StringIO
        pr = cProfile.Profile()
        pr.enable()

        combinations85=[]
        numComb85=[]
        num=100
        nm_it85=range(numModels85)

        print 'making combinations'

        for i in range(len(nm_it85)):
            temp=it.combinations(nm_it85, nm_it85[i]+1)

            top=math.factorial(len(nm_it85))
            b=math.factorial(len(nm_it85)-(nm_it85[i]+1))
            c=math.factorial(nm_it85[i]+1)
            numComb85.append(int(top/(b*c)))
            combinations85.append(temp)
            del temp
    
        # initalize output 
        minNum_nSIF_mean_85=np.nan*np.ones((numYears,nj,ni), dtype=np.float16)
        minNum_nSIF_std_85=np.nan*np.ones((numYears,nj,ni), dtype=np.float16)
        minNum_nSIF_both_85=np.nan*np.ones((numYears,nj,ni), dtype=np.float16)

        minNum_nSIF_85_ks_mean=np.nan*np.ones((numYears,nj,ni), dtype=np.float16)
        minNum_nSIF_85_ks_std=np.nan*np.ones((numYears,nj,ni), dtype=np.float16)
        minNum_nSIF_85_ks_both=np.nan*np.ones((numYears,nj,ni), dtype=np.float16)


        ## then do analysis
        # we should be able to do this outside of the loop written above. just use np.mean on the correct axis.

        mean30=np.nanmean(nSIF85, axis=0)
        std30=np.nanstd(nSIF85, axis=0)

        mean_bool=np.nan*np.ones((numModels85, numYears, nj, ni), dtype=np.int16)
        std_bool=np.nan*np.ones((numModels85, numYears, nj, ni), dtype=np.int16)

        mean_p_ks=np.ones((numModels85, numYears, nj, ni), dtype=np.float16) # set to ones, so that the first comparison had bad pval
        stds_p_ks=np.zeros((numModels85, numYears, nj, ni), dtype=np.float16)

        mean_stat_ks=np.ones((numModels85, numYears, nj, ni), dtype=np.float16) # set to one, so that the first comparison had bad statistic
        stds_stat_ks=np.ones((numModels85, numYears, nj, ni), dtype=np.float16)

        print 'analyzing random samples'        
        for nmm in range(numModels85):
            startTime = datetime.now()
            combos=combinations85[nmm]

            numitter=min(num, numComb85[nmm])
            #print nmm, numitter, numComb85[nmm], startTime

            tempMeans_bool=np.zeros((numYears, nj, ni))
            tempStds_bool=np.zeros((numYears, nj, ni))

            means=np.nan*np.ones((num, numYears, nj, ni), dtype=np.float16)
            stds=np.nan*np.ones((num, numYears, nj, ni), dtype=np.float16)

            # make random inds generator, get numitter values out of scaled uniform distribution
            if numComb85[nmm]>numitter:
                stepInds=np.sort(np.floor(numComb85[nmm]*np.random.rand(numitter)))
                steps=np.hstack((stepInds[0], np.diff(stepInds)))
            else:
                steps=np.ones((numitter))

            for nc in range(numitter):
    
                for s in range(int(steps[nc])): # this is the random part
                    inds=combos.next()

                tempMeans=np.nanmean(nSIF85[inds, :, :, :], axis=0)
                tempStds=np.nanstd(nSIF85[inds, :, :, :], axis=0)
        
                # bootstrap confidence interval
    
                tempMeans_bool+=(tempMeans>(CIsMeanLB))*(tempMeans<(CIsMeanUB)) # number within the range
                tempStds_bool+=(tempStds>(CIsStdsLB))*(tempStds<(CIsStdsUB))

                means[nc, :,:,:]=tempMeans
                stds[nc, :,:,:]=tempStds

            mean_bool[nmm, :,:,:]=tempMeans_bool/float(numitter)
            std_bool[nmm, :,:,:]=tempStds_bool/float(numitter)
            print nmm+1, ' models ', numitter, ' itterations,', 'time needed:', (datetime.now()-startTime)
            if nmm>0:
                print "making KS comparisons ", (datetime.now()-startTime)
                for yr in range(numYears):
                     #print yr, numYears, (datetime.now()-startTime)
                     for nii in range(ni):
                         for njj in range(nj):
                             if mean30[0,njj,nii]<363.:
                                 try:
                                     thisMean=means[:,yr, njj, nii]
                                     thisStd=stds[:,yr, njj, nii]

                                     lastMean=meanslastModel[:,yr, njj, nii]
                                     lastStd=stdslastModel[:,yr, njj, nii]

                                     ksStatMean, pvalMean = stats.ks_2samp(thisMean[~np.isnan(thisMean)], lastMean[~np.isnan(lastMean)])
                                     ksStatStd, pvalStd = stats.ks_2samp(thisStd[~np.isnan(thisStd)], lastStd[~np.isnan(lastStd)])

                                     mean_p_ks[nmm, yr, njj, nii]=pvalMean
                                     stds_p_ks[nmm, yr, njj, nii]=pvalStd
                        
                                     mean_stat_ks[nmm, yr, njj, nii]=ksStatMean
                                     stds_stat_ks[nmm, yr, njj, nii]=ksStatStd
                         
                                 except:
                                     dum=1
            meanslastModel=np.copy(means)
            stdslastModel=np.copy(stds)
            print " done with KS comparisons ", (datetime.now()-startTime)
     

        # choose the last time that the variable is outside of the +/- 2.5% range                    
        # for each year, for each cell, choose the number of models that you need to get within the desired range. 

print 'selecting the number of models needed'
for yr in range(numYears):
    print yr, (datetime.now()-startTime)
    for nii in range(ni):
        for njj in range(nj):
            if mean30[0,njj,nii]<363.:
                # use first method. 95% of values are within the bootstrapped confidence intervals of the 30 year value
    
                sel_mean_bool=mean_bool[:,yr, njj, nii]
                sel_std_bool=std_bool[:,yr, njj, nii]
        
                lastInd_mean=np.where(sel_mean_bool>0.95)[0]
                lastInd_std=np.where(sel_std_bool>0.95)[0]
        
                if len(lastInd_mean)>0:
                    minNum_nSIF_mean_85[yr, njj, nii]=min(30, lastInd_mean[0]+1)
                else:
                    minNum_nSIF_mean_85[yr, njj, nii]=1 # this happens in the sea ice covered area. 
        
                if len(lastInd_std)>0:
                    minNum_nSIF_std_85[yr, njj, nii]=min(30, lastInd_std[0]+1)
                else:
                    minNum_nSIF_std_85[yr, njj, nii]=1 
            
                minNum_nSIF_both_85[yr, njj, nii]=np.max([minNum_nSIF_std_85[yr, njj, nii],minNum_nSIF_mean_85[yr, njj, nii]])  
             
                # use second method, KS test.
                #If the K-S statistic is small or the p-value is high, then we cannot 
                # reject the hypothesis that the distributions of the two samples are the same.
        
                sel_mean_ks=mean_p_ks[:,yr, njj, nii]
                sel_std_ks=stds_p_ks[:,yr, njj, nii]
        
                lastInd_mean_ks=np.where(sel_mean_ks[:-1]>0.05)[0]
                lastInd_std_ks=np.where(sel_std_ks[:-1]>0.05)[0]
        
        
                # find the last time that the p value is at above 0.05
                if len(lastInd_mean_ks)==0: 
                    minNum_nSIF_85_ks_mean[yr, njj, nii]=1
                else:              
                    if len(lastInd_mean_ks)<len(sel_mean_ks[:-1]):
                        minNum_nSIF_85_ks_mean[yr, njj, nii]=min(30, lastInd_mean_ks[-1]+2)
                    else:
                        minNum_nSIF_85_ks_mean[yr, njj, nii]=32 # when convergence is not reached.
    
                if sum(sel_mean_ks)==len(sel_mean_ks):
                    minNum_nSIF_85_ks_mean[yr, njj, nii]=1
    
                if len(lastInd_std_ks)==0: 
                    minNum_nSIF_85_ks_std[yr, njj, nii]=1
                else:   
                    if len(lastInd_std_ks)<len(sel_std_ks[:-1]):
                        minNum_nSIF_85_ks_std[yr, njj, nii]=min(30, lastInd_std_ks[-1]+2)
                    else:
                        minNum_nSIF_85_ks_std[yr, njj, nii]=32
                     
                if sum(sel_std_ks)==len(sel_std_ks):
                    minNum_nSIF_85_ks_std[yr, njj, nii]=1    
                minNum_nSIF_85_ks_both[yr, njj, nii]=np.max([minNum_nSIF_85_ks_std[yr, njj, nii],minNum_nSIF_85_ks_mean[yr, njj, nii]])  
            
                                  
        # write output to netcdf
        print 'writng output to netcdf'

        fAn=Dataset(fnAnalysis, 'w',format='NETCDF4')

        # create all the dimentions, set time to unlimited 
        for k in fnsif.dimensions.keys():
            kdim=fnsif.dimensions[k]
            if kdim.isunlimited()==True:
                fAn.createDimension(k, None)
                print k
            else:
                if k=='time_ensemble':
                    dimSize=180

                else:
                    dimSize=len(kdim)

                fAn.createDimension(k, dimSize)  
            del kdim

        fAn.createDimension('nitt', num)

        # use the netCDF4 instead of pyNIO since it seems to work much better with unlimited variables      
        fAnVars={}
        for key in {'TLAT', 'TLON','latt_bounds','lont_bounds', 'time_ensemble'}:
            print 'creating ', key
            kvar=fnsif.variables[key]
            # the netCDF4 module requires that if a fill value exists, it must be declared when the variable is created.
            try:
                fAnVars[key]=fAn.createVariable(key, kvar.dtype.name[0], kvar.dimensions, fill_value=kvar.missing_value)
            except:
                fAnVars[key]=fAn.createVariable(key, kvar.dtype.name[0], kvar.dimensions)

            # sett all the attribute keys.
            atts = kvar.__dict__
            for attKey in atts.keys():
                if attKey != '_FillValue':
                    setattr(fAn.variables[key],attKey,atts[attKey])  
            del kvar

        # put data into variables, first the ones we are copying over. 
        print 'putting data into standard variables'
        for key in {'TLAT', 'TLON','latt_bounds','lont_bounds'}:
            fAnVars[key][:,:]=fnsif.variables[key][:]

        # output by year
        key='time_ensemble'
        var=fnsif.variables[key][:][-180:]
        fAnVars[key][:]=var

        # key='time_bounds'
        # var=fnsif.variables[key][:][-180:,:]
        # fAnVars[key][:,:]=var


        dataKey=['ks_statitic_mean','ks_statitic_std']
        data=[mean_p_ks,stds_p_ks]
        units=['ks statistic','ks statistic']
        longName=['KS statistic mean',
        'KS statistic std']
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
            fAnVars[key]=fAn.createVariable(key, 'f', ('nm','time_ensemble','nj', 'ni'), fill_value=stdAttributes['_FillValue'])

            # set all the attribute keys.
            for attKey in stdAttributes.keys():
                if attKey != '_FillValue':
                    setattr(fAn.variables[key],attKey,stdAttributes[attKey])
            setattr(fAn.variables[key],'long_name',longName[i])
            setattr(fAn.variables[key],'units',units[i])
            fAnVars[key][:,:,:,:]=data[i] 

        dataKey=['numMembers_mean','numMembers_std','numMembers_both','numMembers_ksTest_mean','numMembers_ksTest_std','numMembers_ksTest_both']
        data=[minNum_nSIF_mean_85,minNum_nSIF_std_85,minNum_nSIF_both_85,minNum_nSIF_85_ks_mean,minNum_nSIF_85_ks_std,minNum_nSIF_85_ks_both]
        units=['number of models','number of models','number of models','number of models','number of models','number of models']
        longName=['Ensemble Members Required for Analysis (simple, mean only)',
        'Ensemble Members Required for Analysis (simple, std only)',
        'Ensemble Members Required for Analysis (simple, both mean and std)',
        'Ensemble Members Required for Analysis (ks-test, mean only)',
        'Ensemble Members Required for Analysis (ks-test, std only)',
        'Ensemble Members Required for Analysis (ks-test, both mean and std)']
        stdAttributes={'_FillValue': np.array([  1.00000002e+30], dtype=float),
        'cell_measures': 'area: tarea',
        'cell_methods': 'time: mean',
        'comment': 'none',
        'coordinates': 'TLON TLAT time',
        'missing_value': np.array([  1.00000002e+30], dtype=float),
        'time_rep': 'averaged'}

        outDims=fnsif.variables['nSIF_ensemble'].dimensions[1:]

        for i in range(len(dataKey)):
            key=dataKey[i]
            print 'creating ', key
            # the netCDF4 module requires that if a fill value exists, it must be declared when the variable is created.
            fAnVars[key]=fAn.createVariable(key, 'f', outDims, fill_value=stdAttributes['_FillValue'])

            # set all the attribute keys.
            for attKey in stdAttributes.keys():
                if attKey != '_FillValue':
                    setattr(fAn.variables[key],attKey,stdAttributes[attKey])
            setattr(fAn.variables[key],'long_name',longName[i])
            setattr(fAn.variables[key],'units',units[i])
            fAnVars[key][:,:,:]=data[i] 

        # output by max per year
        dataKey=['numMembers_mean_allYears','numMembers_std_allYears','numMembers_both_allYears',
        'numMembers_ksTest_mean_allYears','numMembers_ksTest_std_allYears','numMembers_ksTest_both_allYears']
        data=[np.max(minNum_nSIF_mean_85, axis=0),
        np.max(minNum_nSIF_std_85, axis=0),
        np.max(minNum_nSIF_both_85, axis=0),
        np.max(minNum_nSIF_85_ks_mean, axis=0),
        np.max(minNum_nSIF_85_ks_std, axis=0),
        np.max(minNum_nSIF_85_ks_both, axis=0)]

        units=['number of models','number of models','number of models','number of models','number of models','number of models']
        longName=['Ensemble embers Required for Analysis for All Years (simple, mean only)',
        'Ensemble Members Required for Analysis for All Years(simple, std only)',
        'Ensemble Members Required for Analysis for All Years (simple, both mean and std)',
        'Ensemble Members Required for Analysis for All Years (ks test, mean only)',
        'Ensemble Members Required for Analysis for All Years (ks test, std only)',
        'Ensemble Members Required for Analysis for All Years (ks test, both mean and std)']
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
            fAnVars[key]=fAn.createVariable(key, 'f', ('nj', 'ni'), fill_value=stdAttributes['_FillValue'])

            # set all the attribute keys.
            for attKey in stdAttributes.keys():
                if attKey != '_FillValue':
                    setattr(fAn.variables[key],attKey,stdAttributes[attKey])
            setattr(fAn.variables[key],'long_name',longName[i])
            setattr(fAn.variables[key],'units',units[i])
            fAnVars[key][:,:]=data[i] 

  

# pull values for examples:

yr=2000
ny=yr-1920
nj=67
ni=206

combinations85=[]
numComb85=[]
num=100
nm_it85=range(numModels85)

print 'making combinations'
for i in range(len(nm_it85)):
    temp=it.combinations(nm_it85, nm_it85[i]+1)

    top=math.factorial(len(nm_it85))
    b=math.factorial(len(nm_it85)-(nm_it85[i]+1))
    c=math.factorial(nm_it85[i]+1)
    numComb85.append(int(top/(b*c)))
    combinations85.append(temp)
del temp



## then do analysis

print 'analyzing random samples' 
means=np.nan*np.ones((numModels85, num), dtype=np.float16)
stds=np.nan*np.ones((numModels85, num), dtype=np.float16)       
for nmm in range(numModels85):
    startTime = datetime.now()
    combos=combinations85[nmm]

    numitter=min(num, numComb85[nmm])
    print nmm, numitter, numComb85[nmm], startTime

    tempMeans_bool=np.zeros((numYears, nj, ni))
    tempStds_bool=np.zeros((numYears, nj, ni))

    # make random inds generator, get numitter values out of scaled uniform distribution
    if numComb85[nmm]>numitter:
        stepInds=np.sort(np.floor(numComb85[nmm]*np.random.rand(numitter)))
        steps=np.hstack((stepInds[0], np.diff(stepInds)))
    else:
        steps=np.ones((numitter))

    for nc in range(numitter):

        for s in range(int(steps[nc])): # this is the random part
            inds=combos.next()

        tempMeans=np.nanmean(nSIF85[inds, ny, nj, ni], axis=0)
        tempStds=np.nanstd(nSIF85[inds, ny, nj, ni], axis=0)

        # bootstrap confidence interval

        means[nmm, nc]=tempMeans
        stds[nmm, nc]=tempStds

CIsMeanUB=CI.variables['meanUB'][ny+70,nj,ni]
CIsMeanLB=CI.variables['meanLB'][ny+70,nj,ni]
CIsStdsUB=CI.variables['sdUB'][ny+70,nj,ni]
CIsStdsLB=CI.variables['sdLB'][ny+70,nj,ni]

numMean=minNum_nSIF_mean_85[ny, nj, ni]
numStd=minNum_nSIF_std_85[ny, nj, ni]
 

plt.figure()    
plt.hlines(CIsMeanUB, 0, 31)   
plt.hlines(CIsMeanLB, 0, 31) 
plt.hlines(means[-1][0], 0, 31)   
plt.boxplot(means.T)
plt.title('Example mean distributions ('+str(int(numMean))+' members needed), YEAR=' + str(yr)+ ' NI='+str(ni)+ 'NJ='+str(nj))
plt.ylabel('Number of open water days')
plt.xlabel('Number of subsampled ensemble members')
plt.savefig('SI_FigXx_numNeeded_Mean.'+rcpName[ittR]+'.'+nsk+'.pdf', format='pdf')
#plt.show()

plt.figure()    
plt.hlines(CIsStdsUB, 0, 31)   
plt.hlines(CIsStdsLB, 0, 31) 
plt.hlines(stds[-1][0], 0, 31)   
plt.boxplot(stds.T)
plt.title('Example standard deviation distributions('+str(int(numStd))+' members needed), YEAR=' + str(yr)+ ' NI='+str(ni)+ 'NJ='+str(nj))
plt.ylabel('Number of open water days')
plt.xlabel('Number of subsampled ensemble members')
plt.savefig('SI_FigXx_numNeeded_STD.'+rcpName[ittR]+'.'+nsk+'.pdf', format='pdf')

#plt.show()  


        fAn.close()    
        fnsif.close()
        CI.close()
        del fnsif
        del fAn
        del CI
# pr.disable()
# s = StringIO.StringIO()
# sortby = 'cumulative'
# ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
# ps.print_stats()
# print s.getvalue()