print 'importing modules'
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import glob
import pickle

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif']=['Arial']
plt.close('all')

import sys
sys.path.append(u'/Users/katherinebarnhart/git/cesmEnsembleSeaIce')


import haversine as h


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



startTime = datetime.now()

# drew point is at (70.8, -157)
# this corresponds to 
f=Dataset(dirList[0], 'r')


lats=f.variables['TLAT'][:,:]
lons=f.variables['TLON'][:,:]

aice=f.variables['aice_d_monthAvg'][10,:,:]

latSelect=lats
latSelect[aice>100]=np.nan

lonSelect=lons
lonSelect[aice>100]=np.nan


siteOutput={}
siteOutput['SeaOfOkhotsk']={'loc':(57, 147.571006)}
siteOutput['DrewPoint']={'loc':(71, -153.867593)}
siteOutput['CentralHudsonBay']={'loc':(60.680056, -86.563720)}
siteOutput['NorwegianSea']={'loc':(73.865368, -4.054318)}
siteOutput['Bering Sea']={'loc':(58.068153, -176.423135)}
siteOutput['LaptevSea']={'loc':(77.931939, 114.168586)}
siteOutput['EasternSiberianSea']={'loc':(71, 175.296515)}
siteOutput['BeringStrait']={'loc':(67.156341, -168.994759)}
siteOutput['BaffinBay']={'loc':(73.157264, -66.192510)}
siteOutput['BarentsSea']={'loc':(74.713197, 33.804760)}


# add open water points
siteOutput['OceanA']={'loc':(84.891602, -61.583258)}
siteOutput['OceanB']={'loc':(81.001454, 143.201901)}
siteOutput['OceanC']={'loc':(76.128028, -149.298103)}
siteOutput['OceanD']={'loc':(84.118592, 56.190179)}
siteOutput['OceanE']={'loc':(73.410644, -177.071537)}
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
    

numDays=(2100-1849)*365

numSites=len(siteOutput)
thresh=15
keys=siteOutput.keys()
keys=[k for k in keys if k!='obsYears']


# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
map = Basemap(projection='npstere',boundinglat=45,lat_0=90,lon_0=-100,resolution='l')
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))
# contour data over the map.

#x, y = map(lons.flatten(), lats.flatten()) 
#NIf=NI.flatten()
#NJf=NJ.flatten()
#
#for i in range(len(x)):
#    if np.isnan(x[i])==False: 
#        map.plot(x[i],y[i], 'k.')
#        if np.remainder(i, 44)==0:
#            plt.text(x[i],y[i],str(NJf[i])+','+str(NIf[i]),fontsize=9)
#

for key in keys:
    x, y = map(siteOutput[key]['cellLL'][1], siteOutput[key]['cellLL'][0])    
    map.plot(x,y, 'ro')
    text=key
    plt.text(x,y,text,fontsize=9)
plt.title('Assumption Check Locations')
plt.savefig(u'/Users/katherinebarnhart/git/cesmEnsembleSeaIce/checkLocations.pdf', format='pdf')
plt.show()
import csv
outFile=u'/Users/katherinebarnhart/git/cesmEnsembleSeaIce/assumptionCheckPoints.csv'
csvfile=open(outFile, 'wb')  
wr = csv.writer(csvfile, delimiter=',')
wr.writerow(['sitename','ni','nj'])

for key in keys:
    wr.writerow([key, siteOutput[key]['inds']['ni']+1, siteOutput[key]['inds']['nj']+1])

csvfile.close()

# print location list
