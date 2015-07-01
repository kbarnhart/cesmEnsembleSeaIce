#OWextentThroughTime_plots
print 'importing modules'
import numpy as np
import matplotlib.pyplot as plt

import datetime
from matplotlib.mlab import csv2rec

import matplotlib.patches as mpatches

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif']=['Arial']
plt.close('all')
#from matplotlib import rc
#rc('text', usetex=True)
plt.ion()

import sys
sys.path.append('/Users/katherinebarnhart/python/myPyModules/')
import statsForSeaIceBayes as sis

import cProfile, pstats, StringIO
pr = cProfile.Profile()

import pickle
plt.close('all')
output=pickle.load(open("/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/OWExtent.p", "rb" ))                             

BGiceArea=output['BGiceAreaThresh']
BGtime=output['BGtime']
time=output['time']
iceArea=output['iceAreaThresh']

areaConversion=(1000.*1000.)*1000000. #meters square to millions of km square

# Observations
a=csv2rec("/Volumes/Pitcairn/seaicePPF/northernHemisphere/code/NH_seaice_extent_final.csv")

obsYearAll=[]
obsDay=[]
obsMonth=[]
obsArea=[]
obsTime=[]
obsDate=[]
for i in range(1, len(a)):
    obsYearAll.append(float(a[i][0]))
    obsMonth.append(float(a[i][1]))
    obsDay.append(float(a[i][2]))
    obsArea.append(float(a[i][3])*1000.*1000.*1000000.)
    obsTime.append(datetime.datetime(int(a[i][0]),int(a[i][1]),int(a[i][2])))

# get yearly min and max
obsYear=np.unique(obsYearAll)[1:] # remove the few points in 1978
obsMo=[1,2,3,4,5,6,7,8,9,10,11,12]
yearlyArea=[]
obsMin=[]
obsMax=[]
yearlyByMonth=[]
avgByMonth=[]
for i in range(len(obsYear)):
    yr=obsYear[i]
    inds=np.where(obsYearAll==yr)[0]
    data=np.asarray(obsArea)[inds]
    yearlyArea.append(data)
    obsMin.append(data.min())
    obsMax.append(data.max())
    
    months=np.asarray(obsMonth)[inds]
    monthData=[]
    monthAvg=[]
    for j in range(len(obsMo)):
        m=obsMo[j]
        minds=np.where(months==m)[0]
        mdata=data[minds]
        monthData.append(mdata)
        monthAvg.append(mdata.mean())
    yearlyByMonth.append(monthData)
    avgByMonth.append(monthAvg)
avgByMonth=np.asarray(avgByMonth)   
##
obsMin=np.asarray(obsMin)
obsMax=np.asarray(obsMax)

## Model Output
numModel=iceArea.shape[0]
numYear=time.size/365
numBGyear=BGtime.size/365
year=np.unique(np.floor(time/365))[:-1]
BGyear=np.unique(np.floor(BGtime/365))[:-1]

# reshape toget minimums and maximums per model, per time

BGiceReshape=np.reshape(BGiceArea, (numBGyear,365))
BGmin=BGiceReshape.min(axis=-1)
BGmax=BGiceReshape.max(axis=-1)

iceReshape=np.reshape(iceArea, (numModel,numYear,365))
iceMin=iceReshape.min(axis=-1)
iceMax=iceReshape.max(axis=-1)

# now for monthly means

daysPerMonth=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
dayOfMonth=[0,31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
cumDaysInYear=np.cumsum(daysPerMonth)
monthStartStop=np.cumsum(dayOfMonth)
iceMonthMean=np.ones((30, 251, 12))
BGmonthMean=np.ones((BGiceArea.size/365, 12))

for i in range(len(daysPerMonth)):
    sInd=monthStartStop[i]
    eInd=monthStartStop[i+1]
    
    BGmonthMean[:,i]=BGiceReshape[:,sInd:eInd].mean(axis=-1)
    iceMonthMean[:,:,i]=iceReshape[:,:,sInd:eInd].mean(axis=-1)
    

############################################
############################################
############################################
###########################################
startTime=datetime.datetime.now()
pr.enable()

# Bayesian and non bayesian mean comparison

# bayesian inputs
dmesh = 0.1# 5000*1000*1000 # 10k sq km
mesh_max=365 #1.5*BGdata.max()
mesh_min=0

startYear_ind=72
endYear_ind=250

M=1000
alpha=5

num_samples=30#iceMin.shape[0]
numYear=251
time=np.arange(1850, 2101)

monthText=['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
names=[]#=['Minimum Extent', 'Maximum Extent']
bgs=[]#[BGmin, BGmax]
d=[]#=[iceMin, iceMax]
ob=[]#=[obsMin, obsMax]


shiftYear_newDLM=[2005, 2002	,1999	,1989,	1988,	1998	,1987	,1981	,1983	,1982	,1998	,2000]
emergeYear_newDLM=[2025,	2028,	2030,	2034,	2032,	2027,	2017,	2017,	2016,	2015,	2017,	2017]
###############################

for i in range(12):
    bgs.append(BGmonthMean[:,i])
    d.append(iceMonthMean[:,:,i])
    ob.append(avgByMonth[:,i])
    names.append(monthText[i]+ ' Average Extent')
    
for di in range(len(names)):
    BGdata=bgs[di]
    data=d[di]
    obs=ob[di]
    
    if sum(np.isnan(BGdata))>0:
        ninds=np.where(np.isnan(BGdata))[0]
        for ni in ninds:
            BGdata[ni]=np.nanmean(BGdata)

    
    dmesh = 5000*1000*1000 # 10k sq km
    mesh_max=1.5*np.nanmax(BGdata)
    mesh_min=0
    

      
    # auto generated things:
    numX=int(np.floor((mesh_max-mesh_min)/dmesh))
    mesh=np.linspace(mesh_min,mesh_max,num=numX)
    
    #pi_temp, expectedValue_temp, sample_mean_temp, sample_var_temp, posterior_var_temp, bg_mean_temp, ub_temp, lb_temp, e_year_temp, e_ind_temp =sis.bayes(BGdata, data, mesh, M, alpha, time, startYear_ind, endYear_ind)
                

    fig5=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig5.patch.set_alpha(0.0)
    ax5 = fig5.add_subplot(111) 
    
    plt.fill_between([1920,2110], BGdata.min()/areaConversion,BGdata.max()/areaConversion,color='k',alpha=.2) 
    l_bgmean = plt.hlines(BGdata.mean()/areaConversion, 1920, 2100, colors='k' ,lw=0.5, label='Background Mean', linestyle='dotted')
    
    plt.fill_between(year, data.min(axis=0)/areaConversion,data.max(axis=0)/areaConversion,color='CornflowerBlue',alpha=.6, label='Model Range') 
    
    l_modmean, = plt.plot(year, data.mean(axis=0)/areaConversion, color='k', label='Model Mean Value') 
    #l_bev, = plt.plot(time,expectedValue_temp/areaConversion, 'b', label='Bayesian Expected Value')
       
    plt.vlines(shiftYear_newDLM[di], -0.1,20, color='k', alpha=0.5, lw=2, linestyle='dashed')
    plt.vlines(emergeYear_newDLM[di],  -0.1,20, color='k', alpha=0.5, lw=2, linestyle='dashed'  )

    l_obs, = plt.plot(obsYear, obs/areaConversion,  'o', markeredgecolor='#4C2E0F', markerfacecolor='#FF9933', markeredgewidth=0.1, markersize=3, label='Observations')


    #plt.fill_between(year, lb_temp/areaConversion,ub_temp/areaConversion,color='Yellow',alpha=.6, label='Model Range') 

    
    bluePatch= mpatches.Patch(color='CornflowerBlue',alpha=.6, label='Model Ensemble Range') 
    greyPatch= mpatches.Patch(color='k',alpha=.2, label='Background Range') 
    #yellowPatch= mpatches.Patch(color='Yellow',alpha=.6, label='Bayesian Confidence Interval') 
    
    plt.legend(handles=[l_bgmean, greyPatch, l_modmean, bluePatch,l_obs], fontsize=6, fancybox=True, numpoints=1, loc=3)
    
    ax5.set_xlim([1920, 2100])
    ax5.set_ylim([-0.1,18])
    
    ax5.set_ylabel('Sea Ice Extent [millions of km$^{2}$]', fontsize=8)
    ax5.set_xlabel('Year', fontsize=8)
    plt.text(2030, 15 , 'Emergence Year: '+str(emergeYear_newDLM[di]), fontsize=6)
    plt.text(2030, 10 , 'Shift Year: '+str(shiftYear_newDLM[di]), fontsize=6)
    
    ax5.set_title(names[di]+'\n CESM-LE (1920-2100)', fontsize=10)  
    plt.tight_layout()
    plt.savefig('/Users/katherinebarnhart/Desktop/MANUSCRIPTS/2014PPFSeaIce/SI_figs/exitYearFig_pub_v2_'+ names[di]+'.pdf', format='pdf')


    #plt.show()    