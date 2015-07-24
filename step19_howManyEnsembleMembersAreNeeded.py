# step 11 determine how many ensemble members you need to do the analysis correctly.

print 'importing modules'
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import glob
import nio
from scipy import stats
import math
import itertools as it
from scipy.stats import ks_2samp
import scikits.bootstrap as bootstrap
import scipy  

#import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['font.sans-serif']=['Arial']

startTime = datetime.now()

path=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/'

pathIn=u'/Volumes/Pitcairn/seaicePPF/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/'
pathIn45=u'/Volumes/Pitcairn/seaicePPF/p/cesm0005/CESM-CAM5-BGC-ME/ice/proc/tseries/daily/aice_d/'
pathOut=u'/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/'


# hard code in the names for the run parts
bgKey=u'B1850C5CN'
runPart1key=u'B20TRC5CNBDRD'
runParts23keyRCP85=u'BRCP85C5CNBDRD'
runParts23keyRCP45=u'BRCP45C5CNBDRD'

#os.chdir(path)
#dirList45=glob.glob(pathOut+

dirList85=glob.glob(path+'*BRCP85C5CNBDRD*nh*.Analysis.nc')
dirListBG=glob.glob(path+'b.e11.B1850C5CN.*.nh.Analysis.nc')#path=u'/Users/katherinebarnhart/Desktop/RESEARCH/seaIcePPF/'

analysisFiles85=[]    
for fn in dirList85:
    analysisFiles85.append(Dataset(fn, 'r'))
    print fn, (datetime.now()-startTime)
    
analysisFilesBG=[]    
for fn in dirListBG:
    analysisFilesBG.append(Dataset(fn, 'r'))
    print fn, (datetime.now()-startTime)    
    
numModels85=len(dirList85)
numYears=(2100-1920)

ni=320
nj=104
nm=numModels85

numSIF85=np.nan*np.ones((numModels85,numYears,nj,ni))
first85=np.nan*np.ones((numModels85,numYears,nj,ni))
last85=np.nan*np.ones((numModels85,numYears,nj,ni))

i=0
for fn in dirList85:
    f=Dataset(fn, 'r')
    
    ns=f.variables['numSIF'][-180:,:,:]
    numSIF85[i,:,:,:]=ns
    
    fs=f.variables['firstDOY'][-180:,:,:]
    first85[i,:,:,:]=fs
    
    ls=f.variables['lastDOY'][-180:,:,:]
    last85[i,:,:,:]=ls
    
    print fn, (datetime.now()-startTime)
    i+=1
    f.close()

nSIF85=np.array(numSIF85,dtype=np.float16, order='C')
first85=np.array(first85)
last85=np.array(last85)

nSIF85[nSIF85>1000]=np.nan

del first85
del last85

import cProfile, pstats, StringIO
pr = cProfile.Profile()
pr.enable()

combinations85=[]
numComb85=[]
num=100
nm_it85=range(numModels85)

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

numCI=1000
CIsMeanUB=np.nan*np.zeros((numYears, nj, ni), dtype=np.float16)
CIsMeanLB=np.nan*np.zeros((numYears, nj, ni), dtype=np.float16)
CIsStdsUB=np.nan*np.zeros((numYears, nj, ni), dtype=np.float16)
CIsStdsLB=np.nan*np.zeros((numYears, nj, ni), dtype=np.float16)

for yr in range(numYears):
    print yr, (datetime.now()-startTime)
    for nii in range(ni):
        for njj in range(nj): 
            vals=nSIF85[:, yr, njj, nii]
            if np.std(vals)>0:
                if mean30[0,njj,nii]<363.:
                    CIsMean = bootstrap.ci(vals, scipy.mean, n_samples=numCI, alpha=0.05) 
                    CIsStds = bootstrap.ci(vals, scipy.std, n_samples=numCI, alpha=0.05) 
                    CIsMeanUB[yr, njj, nii]=CIsMean.max()
                    CIsMeanLB[yr, njj, nii]=CIsMean.min()
                    CIsStdsUB[yr, njj, nii]=CIsStds.max()
                    CIsStdsLB[yr, njj, nii]=CIsStds.min()
               
                
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
    
        tempMeans_bool+=(tempMeans>(CIsMean.min()))*(tempMeans<(CIsMean.max())) # number within the range
        tempStds_bool+=(tempStds>(CIsStds.min()))*(tempStds<(CIsStds.max()))

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
            	# use first method. 95% of values are within +/- 2.5% of 30 year value
            
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
                
                if len(lastInd_mean_ks)<len(sel_mean_ks[:-1]):
                    minNum_nSIF_85_ks_mean[yr, njj, nii]=min(30, lastInd_mean_ks[-1]+2)
                else:
                    minNum_nSIF_85_ks_mean[yr, njj, nii]=32 # when convergence is not reached.
                
                if sum(sel_mean_ks)==len(sel_mean_ks):
                    minNum_nSIF_85_ks_mean[yr, njj, nii]=1
                
                
                
                if len(lastInd_std_ks)<len(sel_mean_ks[:-1]):
                    minNum_nSIF_85_ks_std[yr, njj, nii]=min(30, lastInd_std_ks[-1]+2)
                else:
                    minNum_nSIF_85_ks_std[yr, njj, nii]=32 
                if sum(sel_std_ks)==len(sel_std_ks):
                    minNum_nSIF_85_ks_std[yr, njj, nii]=1    
                minNum_nSIF_85_ks_both[yr, njj, nii]=np.max([minNum_nSIF_85_ks_std[yr, njj, nii],minNum_nSIF_85_ks_mean[yr, njj, nii]])  
                                     
# write output to netcdf
f=nio.open_file(dirList85[0], 'r')
fnAnalysis=pathOut+'numMembersNeeded.nc'
fAn=Dataset(fnAnalysis, 'w',format='NETCDF4')

# create all the dimentions, set time to unlimited 
for k in f.dimensions.keys():
    if f.unlimited(k)==True:
        fAn.createDimension(k, None)
        print k
    else:
        fAn.createDimension(k, f.dimensions[k])  

fAn.createDimension('nitt', num)
fAn.createDimension('nm', numModels85)
    
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

# output by year
key='time'
var=f.variables[key][:][-180:]
fAnVars[key][:]=var
key='time_bounds'
var=f.variables[key][:][-180:,:]
fAnVars[key][:,:]=var


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
    fAnVars[key]=fAn.createVariable(key, 'f', ('nm','time','nj', 'ni'), fill_value=stdAttributes['_FillValue'])
    
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

for i in range(len(dataKey)):
    key=dataKey[i]
    print 'creating ', key
    # the netCDF4 module requires that if a fill value exists, it must be declared when the variable is created.
    fAnVars[key]=fAn.createVariable(key, 'f', f.variables['numSIF'].dimensions, fill_value=stdAttributes['_FillValue'])
    
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

fAn.close()    
f.close()
del f
del fAn


#import matplotlib.pylab as plt
#    
#plt.figure()    
#plt.hlines(means[-1][0]*1.025, 0, 31)   
#plt.hlines(means[-1][0]*0.975, 0, 31) 
#plt.hlines(means[-1][0], 0, 31)   
#plt.boxplot(means)
#plt.title('Mean Distribution')
#plt.ylabel('Number of OW days')
#plt.xlabel('Number of ensemble members')
#plt.show()
#
#plt.figure()    
#plt.hlines(stds[-1][0]*1.025, 0, 31)   
#plt.hlines(stds[-1][0]*0.975, 0, 31) 
#plt.hlines(stds[-1][0], 0, 31)   
#plt.boxplot(stds)
#plt.title('Standard Deviation Distribution')
#plt.ylabel('Number of OW days')
#plt.xlabel('Number of ensemble members')
#plt.show()  

pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print s.getvalue()