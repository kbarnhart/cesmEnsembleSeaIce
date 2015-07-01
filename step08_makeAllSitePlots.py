#from netCDF4 import Dataset
#import glob
#from datetime import datetime
#import numpy as np
#path=u'/Volumes/Svalbard/seaIcePPF/'
##path=u'/Users/katherinebarnhart/Desktop/RESEARCH/seaIcePPF/'
##os.chdir(path)
#dirList=glob.glob(path+'b.e11.f09_g16.001.*.timeseries.Analysis.nc')
#dirListOrig=glob.glob(path+'b.e11.f09_g16.001.*.timeseries.nc')
#
#fn1850=u'/Volumes/Svalbard/seaIcePPF/b.e11.B1850C5CN.f09_g16.005.cice.h1.aice_d_nh.timeseries.nc'
#
#startTime = datetime.now()
#
## drew point is at (70.8, -157)
## this corresponds to 
#f=Dataset(dirList[0], 'r')
#lats=f.variables['TLAT'][:,:]
#lons=f.variables['TLON'][:,:]
#
#
#nj=104
#ni=320
#
#NI, NJ= np.meshgrid(range(ni), range(nj))
#siteOutput['']=()
#siteOutput['']=()
#siteOutput['']=()
#siteOutput['']=()
#siteOutput['']=()
#siteOutput['']=()
#siteOutput['']=()
 
from scipy.io.wavfile import write as wavwrite
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif']=['Arial']
plt.close('all')
#plt.ion()
import pickle
import numpy as np
from datetime import datetime

siteOutput=pickle.load(open("/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/siteOutput.p", "rb" ) )                                    
startTime = datetime.now()

numDays=(2100-1849)*365
numModels=siteOutput['Drew Point']['nSIF'].shape[1]

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
    text=key+'('+str(siteOutput[key]['inds']['nj'])+','+str(siteOutput[key]['inds']['ni'])+')'
    plt.text(x,y,text,fontsize=9)
plt.title('Site Locations')
plt.show()
plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/siteLocations.pdf', format='pdf')

print (datetime.now()-startTime)
C=[[255,255,255],[255,255,255],[254,255,255],[254,255,255],[254,255,255],[253,255,255],[253,255,255],[253,255,255],[252,255,255],[252,255,255],[252,255,255],[251,255,255],[251,255,255],[251,255,255],[250,255,255],[250,255,255],[250,255,255],[249,255,255],[249,255,255],[249,255,255],[248,255,255],[248,255,255],[248,255,255],[247,255,255],[247,255,255],[247,255,255],[246,255,255],[246,255,255],[245,255,255],[245,255,255],[245,255,255],[244,255,255],[244,255,255],[244,255,255],[243,255,255],[243,255,255],[243,255,255],[242,255,255],[242,255,255],[242,255,255],[241,255,255],[241,255,255],[241,255,255],[240,255,255],[240,255,255],[240,255,255],[239,255,255],[239,255,255],[239,255,255],[238,255,255],[238,255,255],[238,255,255],[237,255,255],[237,255,255],[237,255,255],[236,255,255],[236,255,255],[236,255,255],[235,255,255],[235,255,255],[235,255,255],[234,255,255],[234,255,255],[234,255,255],[233,255,255],[232,254,255],[232,254,255],[231,254,255],[230,253,255],[230,253,255],[229,252,255],[228,252,255],[227,251,255],[227,251,255],[226,251,255],[225,250,255],[225,250,255],[224,249,255],[223,249,255],[222,249,255],[222,248,255],[221,248,255],[220,247,255],[220,247,255],[219,246,255],[218,246,255],[217,246,255],[217,245,255],[216,245,255],[215,244,255],[215,244,255],[214,243,255],[213,243,255],[212,243,255],[212,242,255],[211,242,255],[210,241,255],[210,241,255],[209,240,255],[208,240,255],[207,240,255],[207,239,255],[206,239,255],[205,238,255],[205,238,255],[204,238,255],[203,237,255],[202,237,255],[202,236,255],[201,236,255],[200,235,255],[200,235,255],[199,235,255],[198,234,255],[197,234,255],[197,233,255],[196,233,255],[195,232,255],[195,232,255],[194,232,255],[193,231,255],[192,231,255],[192,230,255],[191,230,255],[190,229,255],[190,229,255],[189,229,255],[188,228,255],[187,227,255],[184,226,255],[182,224,254],[180,223,254],[177,221,254],[175,220,254],[172,218,254],[170,217,253],[168,215,253],[165,214,253],[163,212,253],[161,211,252],[158,209,252],[156,208,252],[153,206,252],[151,205,252],[149,203,251],[146,202,251],[144,200,251],[142,199,251],[139,197,251],[137,196,250],[135,194,250],[132,193,250],[130,191,250],[128,190,249],[125,188,249],[123,187,249],[120,185,249],[118,184,249],[116,182,248],[113,181,248],[111,180,248],[109,178,248],[106,177,247],[104,175,247],[101,174,247],[99,172,247],[97,171,247],[94,169,246],[92,168,246],[90,166,246],[87,165,246],[85,163,246],[83,162,245],[80,160,245],[78,159,245],[76,157,245],[73,156,244],[71,154,244],[68,153,244],[66,151,244],[64,150,244],[61,148,243],[59,147,243],[57,145,243],[54,144,243],[52,142,242],[49,141,242],[47,139,242],[45,138,242],[42,136,242],[40,135,241],[38,133,241],[37,131,238],[36,129,234],[35,127,231],[35,125,227],[34,123,223],[34,121,219],[33,119,216],[33,117,212],[32,115,208],[31,113,204],[31,110,200],[30,108,197],[30,106,193],[29,104,189],[29,102,185],[28,100,182],[27,98,178],[27,96,174],[26,94,170],[26,92,166],[25,90,163],[24,88,159],[24,85,155],[23,83,151],[23,81,147],[22,79,144],[22,77,140],[21,75,136],[20,73,132],[20,71,129],[19,69,125],[19,67,121],[18,65,117],[17,63,113],[17,60,110],[16,58,106],[16,56,102],[15,54,98],[15,52,95],[14,50,91],[13,48,87],[13,46,83],[12,44,79],[12,42,76],[11,40,72],[10,38,68],[10,35,64],[9,33,61],[9,31,57],[8,29,53],[8,27,49],[7,25,45],[6,23,42],[6,21,38],[5,19,34],[5,17,30],[4,15,26],[3,13,23],[3,10,19],[2,8,15],[2,6,11],[1,4,8],[1,2,4],[0,0,0]]
colors=np.asarray(C)
colors=np.flipud(colors)/255.
cmap = mpl.colors.ListedColormap(colors)
tix=[70, 100, 150, 200, 250]
ticlabels=[str(ti+1850) for ti in tix]

boxPos=np.arange(1850, 2101)
boxtix=[1850, 1900, 1950, 2000, 2050, 2100]
bticlabels=[str(ti) for ti in boxtix]

cmap2=plt.get_cmap('Greys_r')                                                                                                                                                                                                                                                                                           
for key in keys: 
    print 'making plts', key
    plt.close('all')
    sic=siteOutput[key]['sic']                        
    sic1850=siteOutput[key]['sic1850']
    sic[sic>101]=np.nan    
    sic1850[sic1850>101]=np.nan
    
    nSIF=siteOutput[key]['nSIF']
    nSIF1850=siteOutput[key]['nSIF1850']
    nSIFobs=siteOutput[key]['nsifobs']
    
    obsYears=siteOutput['obsYears']
    
    nyears=numDays/365
    years1850=len(sic1850)/365    
    sicRS=np.reshape(sic[:,0], (nyears,365))
    allsicRS=np.reshape(sic, (nyears,365, numModels))
    meanAllsic=np.nanmean(allsicRS, axis=-1)
    sic1850RS=np.reshape(sic1850, (years1850,365))    
    
   
    nanMask=np.isnan(allsicRS)
    isIce=allsicRS>thresh
    isIce[nanMask]=np.nan
    probIce=np.nanmean(isIce, axis=-1)
    
    p25=[]
    p75=[]
    for yit in range(len(boxPos)):
        p25.append(np.percentile(nSIF[yit,:], 25))
        p75.append(np.percentile(nSIF[yit,:], 75))
    
    p25_1850= np.percentile(nSIF1850, 25)
    p75_1850= np.percentile(nSIF1850, 75)     
         
    data1=(sic[:,0]>15)*1.
    data2=(sic1850>15)*1.
    wavwrite('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/OWnoise1850-2100.'+key+'.wav', 365*8, data1)
    wavwrite('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/OWnoise1850BG'+key+'.wav', 365*8, data2)

    fig8, ax8 = plt.subplots(1, 2, sharey=True, sharex=False, num=None, figsize=(7,3.5), dpi=300, facecolor='w', edgecolor='w')
    fig8.sca(ax8[0])
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig8.sca(ax8[1])
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig8.patch.set_alpha(0.0)
    
    #ax5.set_aspect('equal', 'datalim')
    fig8.sca(ax8[1])
    plt.imshow(meanAllsic.T, vmin=0, vmax=100, cmap=cmap)
    #plt.colorbar()
    ax8[1].set_ylabel('Day of Year', fontsize=8)
    ax8[1].set_xlabel('Year', fontsize=8)
    ax8[1].set_title('Ensemble Mean Sea Ice Concentration at\n '+key+' (1920-2100)', fontsize=10)
    ax8[1].set_xticks(tix)
    ax8[1].set_xticklabels(ticlabels)
    ax8[1].set_xlim([70,250])
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=6)
    cb.set_label('Sea Ice Concentration', size=6)
    
    fig8.sca(ax8[0])
    ax8[0].set_title('1850 constant climate at \n'+key+' ('+str(years1850)+' model years)', fontsize=10)    
    plt.imshow(sic1850RS.T, vmin=0, vmax=100, cmap=cmap)
    #plt.colorbar()
    ax8[0].set_ylabel('Day of Year', fontsize=8)
    ax8[0].set_xlabel('Year', fontsize=8)
    ax8[0].set_xlim([0,1797])    
    #plt.tight_layout()
    plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/allSIC.'+key+'.pdf', format='pdf')    
            
    fig3=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig3.patch.set_alpha(0.0)
    ax3 = fig3.add_subplot(111)
    #ax3.set_aspect('equal', 'datalim')
    plt.imshow(sic1850RS.T, vmin=0, vmax=100, cmap=cmap)
    #plt.colorbar()
    ax3.set_ylabel('Day of Year', fontsize=8)
    ax3.set_xlabel('Year', fontsize=8)
    ax3.set_title('1850 constant climate at \n'+key+' ('+str(years1850)+' model years)', fontsize=10)
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=6)
    cb.set_label('Sea Ice Concentration', size=6)
    plt.tight_layout()
    plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/constant1850SIC.'+key+'.pdf', format='pdf')
    
    fig4=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig4.patch.set_alpha(0.0)
    ax4 = fig4.add_subplot(111)
    #ax4.set_aspect('equal', 'datalim')
    plt.imshow(sicRS.T, vmin=0, vmax=100, cmap=cmap)
    #plt.colorbar()
    ax4.set_ylabel('Day of Year', fontsize=8)
    ax4.set_xlabel('Year', fontsize=8)
    ax4.set_title('Ensemble Member #1 Sea Ice Concentration at\n '+key+' (1920-2100)', fontsize=10)
    ax4.set_xlim([70,250])
    ax4.set_xticks(tix)
    ax4.set_xticklabels(ticlabels)
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=6)
    cb.set_label('Sea Ice Concentration', size=6)
    plt.tight_layout()
    plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/run001SIC.'+key+'.pdf', format='pdf')
    #plt.show()
    
    fig5=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig5.patch.set_alpha(0.0)
    ax5 = fig5.add_subplot(111)
    #ax5.set_aspect('equal', 'datalim')
    plt.imshow(meanAllsic.T, vmin=0, vmax=100, cmap=cmap)
    #plt.colorbar()
    ax5.set_ylabel('Day of Year', fontsize=8)
    ax5.set_xlabel('Year', fontsize=8)
    ax5.set_title('Ensemble Mean Sea Ice Concentration at\n '+key+' (1920-2100)', fontsize=10)
    ax5.set_xlim([70,250])
    ax5.set_xticks(tix)
    ax5.set_xticklabels(ticlabels)
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=6)
    cb.set_label('Sea Ice Concentration', size=6)
    plt.tight_layout()
    plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/ensembleMeanSIC.'+key+'.pdf', format='pdf')
    #plt.show()
    
    fig6=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig6.patch.set_alpha(0.0)
    ax6 = fig6.add_subplot(111)
    #ax6.set_aspect('equal', 'datalim')
    plt.imshow(probIce.T, vmin=0, vmax=1, cmap=cmap2)
    #plt.colorbar()
    ax6.set_ylabel('Day of Year', fontsize=8)
    ax6.set_xlabel('Year', fontsize=8)
    ax6.set_title('Ensemble Probability SIC>15%\n at '+key+' (1920-2100)', fontsize=10)
    ax6.set_xticks(tix)
    ax6.set_xticklabels(ticlabels)
    ax6.set_xlim([70,250])
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=6)
    cb.set_label('Probability', size=6)
    plt.tight_layout()
    plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/ensembleProbSI.'+key+'.pdf', format='pdf')
     
    flierprops = dict(marker='+', markerfacecolor='blue', markeredgecolor='blue', markersize=6,
                    linestyle='none')
    medianprops = dict(linestyle='-', linewidth=1, color='red')
    whiskerprops=dict(linestyle='-', linewidth=0.5, color='black')
    meanpointprops = dict(marker='o', markeredgecolor='black',
                        markerfacecolor='grey')
    
    fig2=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig2.patch.set_alpha(0.0)
    ax2 = fig2.add_subplot(111) 
    #ax2.set_aspect('equal', 'datalim')
    plt.fill_between([1900,2110], nSIF1850.min(),nSIF1850.max(),color='k',alpha=.2) 
    ax2.hlines(nSIF1850.min(), 1900, 2105, colors='k', lw=0.5)
    ax2.hlines(nSIF1850.mean(), 1900, 2105, colors='k' ,lw=0.5)
    ax2.hlines(nSIF1850.max(), 1900, 2105, colors='k', lw=0.5)
    
    ax2.boxplot(nSIF.T, positions=boxPos, whis=[10, 90], showbox=False, showcaps=False,flierprops=flierprops, medianprops=medianprops,whiskerprops=whiskerprops)
    ax2.plot(obsYears, nSIFobs, 'yo', markersize=6)
    ax2.set_xticks(boxtix)
    ax2.set_xticklabels(bticlabels)
    ax2.set_xlim([1915, 2105])
    ax2.set_ylim([-10, 365])
    #ax2.set_xlim([70,250])
    ax2.set_ylabel('Number of Sea Ice Free Days', fontsize=8)
    ax2.set_xlabel('Year', fontsize=8)
    ax2.set_title('Number of Open Water Days\n '+key+' (1850-2100)', fontsize=10)    
    plt.tight_layout()
    plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/numSIF.'+key+'.pdf', format='pdf')
    
    fig7=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig7.patch.set_alpha(0.0)
    ax7 = fig7.add_subplot(111) 
    #ax2.set_aspect('equal', 'datalim')
    plt.fill_between([1900,2110], nSIF1850.min(),nSIF1850.max(),color='k',alpha=.2) 

    ax7.hlines(nSIF1850.min(), 1900, 2105, colors='k', lw=0.5)
    ax7.hlines(nSIF1850.mean(), 1900, 2105, colors='k' ,lw=0.5)
    ax7.hlines(nSIF1850.max(), 1900, 2105, colors='k', lw=0.5)
    
    ax7.boxplot(nSIF.T,positions=boxPos, whis=[10, 90], showbox=False, showcaps=False,flierprops=flierprops, medianprops=medianprops,whiskerprops=whiskerprops)
    #ax7.boxplot(nSIF1850, positions=[1850], whis=[10, 90], showbox=False, showcaps=False,flierprops=flierprops, medianprops=medianprops,whiskerprops=whiskerprops)#, showmeans=True)
    ax7.plot(obsYears, nSIFobs, 'yo', markersize=6)
    ax7.set_xticks(boxtix)
    ax7.set_xticklabels(bticlabels)
    ax7.set_xlim([1976, 2015])
    ax7.set_ylim([-10, 365])
    ax7.set_ylabel('Number of Sea Ice Free Days', fontsize=8)
    ax7.set_xlabel('Year', fontsize=8)
    ax7.set_title('Number of Open Water Days\n '+key+' (1979-2013)', fontsize=10)    
    plt.tight_layout()
    plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/numSIF.obs.'+key+'.pdf', format='pdf')
    
    
    
    fig6=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig6.patch.set_alpha(0.0)
    ax6 = fig6.add_subplot(111)
    #ax6.set_aspect('equal', 'datalim')
    plt.imshow(probIce.T, vmin=0, vmax=1, cmap=cmap2)
    #plt.colorbar()
    ax6.set_ylabel('Day of Year', fontsize=8)
    ax6.set_xlabel('Year', fontsize=8)
    ax6.set_title('Ensemble Probability SIC>15%\n at '+key+' (1920-2100)', fontsize=10)
    ax6.set_xticks(tix)
    ax6.set_xticklabels(ticlabels)
    ax6.set_xlim([70,250])
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=6)
    cb.set_label('Probability', size=6)
    plt.tight_layout()
    plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/ensembleProbSI.'+key+'.pdf', format='pdf')
     
    flierprops = dict(marker='+', markerfacecolor='blue', markeredgecolor='blue', markersize=6,
                    linestyle='none')
    medianprops = dict(linestyle='-', linewidth=1, color='red')
    whiskerprops=dict(linestyle='-', linewidth=0.5, color='black')
    meanpointprops = dict(marker='o', markeredgecolor='black',
                        markerfacecolor='grey')
    
    fig2=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig2.patch.set_alpha(0.0)
    ax2 = fig2.add_subplot(111) 
    ax2.hlines(nSIF1850.mean(), 1900, 2105, colors='k' ,lw=0.5)
    ax2.hlines(nSIF1850.std(), 1900, 2105, colors='k' ,lw=0.5, linestyle='dotted')
    
    ax2.plot(boxPos,nSIF.mean(axis=-1), 'r')
    ax2.plot(boxPos,nSIF.std(axis=-1), 'b')
    #ax2.boxplot(nSIF1850, positions=[1850], whis=[10, 90], showbox=False, showcaps=False,flierprops=flierprops, medianprops=medianprops,whiskerprops=whiskerprops)#, showmeans=True)
    ax2.plot(obsYears, nSIFobs, 'yo', markersize=6)
    ax2.set_xticks(boxtix)
    ax2.set_xticklabels(bticlabels)
    ax2.set_xlim([1915, 2105])
    ax2.set_ylim([0, 365])
    #ax2.set_xlim([70,250])
    ax2.set_ylabel('Number of Sea Ice Free Days', fontsize=8)
    ax2.set_xlabel('Year', fontsize=8)
    ax2.set_title('Number of Open Water Days\n '+key+' (1850-2100)', fontsize=10)    
    plt.tight_layout()
    plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/numSIF.meanStd.'+key+'.pdf', format='pdf')


        
    fig2=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=6)
    fig2.patch.set_alpha(0.0)
    ax2 = fig2.add_subplot(111) 
    plt.fill_between([1900,2110], nSIF1850.min(),nSIF1850.max(),color='k',alpha=.2) 
    plt.fill_between(boxPos, p25_1850,p75_1850,color='w',alpha=.4) 
    
    ax2.hlines(p75_1850, 1900, 2105, colors='k', lw=0.5, linestyle='dotted',alpha=.4)
    ax2.hlines(nSIF1850.mean(), 1900, 2105, colors='k' ,lw=1,alpha=.4)
    ax2.hlines(p25_1850, 1900, 2105, colors='k', lw=0.5, linestyle='dotted',alpha=.4)
    
    plt.fill_between(boxPos, nSIF.min(axis=-1),nSIF.max(axis=-1),color='CornflowerBlue',alpha=.6) 
    plt.fill_between(boxPos, p25,p75,color='w',alpha=.5,edgecolor='w') 
    
    #plt.fill_between(boxPos, p25,p25,color='w') 
    #plt.fill_between(boxPos, p75,p75,color='w') 
    
    #ax2.plot(boxPos,p25, 'k',lw=0.5)
    #ax2.plot(boxPos,p75, 'k',lw=0.5)
    ax2.plot(boxPos,nSIF.mean(axis=-1), 'k')
    
    ax2.plot(obsYears, nSIFobs, 'o', markeredgecolor='#4C2E0F', markerfacecolor='#FF9933', markeredgewidth=0.1, markersize=3)
    ax2.set_xticks(boxtix)
    ax2.set_xticklabels(bticlabels)
    ax2.set_xlim([1920, 2100])
    ax2.set_ylim([0, 365])
    #ax2.set_xlim([70,250])
    ax2.set_ylabel('Number of Sea Ice Free Days', fontsize=8)
    ax2.set_xlabel('Year', fontsize=8)
    ax2.set_title('Number of Open Water Days\n '+key+' (1850-2100)', fontsize=10)    
    plt.tight_layout()
    plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/numSIF.percentiles.'+key+'.pdf', format='pdf')

# show how distributions change through time at each site (by decade?). 

# first a pdf
    try:
        bins=np.arange(0, 367)
        nPlot=19
        legend=[]
        fig1=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
        plt.tick_params(axis='both', which='major', labelsize=6)
        plt.tick_params(axis='both', which='minor', labelsize=6)
        fig1.patch.set_alpha(0.0)
        ax1 = fig1.add_subplot(111)
        #ax1.set_aspect('equal', 'datalim')
        ax1.set_color_cycle(plt.cm.jet(np.linspace(0,1,nPlot+1)))
        legend.append(1850)
        n, bins, patches = plt.hist(nSIF1850, bins, normed=1, histtype='step', linewidth=1, edgecolor='k' )
        for hItt in range(nPlot): # one per decade
            selnsif=nSIF[70+hItt*10,:]
            n, bins, patches = ax1.hist(selnsif, bins, normed=1, histtype='step',linewidth=0.5,)
            legend.append(1920+hItt*10)
        
        ax1.set_xlabel('Number of open water days \n', fontsize=8)
        ax1.set_ylabel('Probability', fontsize=8)
        plt.legend(legend,fontsize=6)
        ax1.set_title('PDFs of open water at \n'+key+' (1920-2100)', fontsize=10)    
        ax1.set_ylim([0,1])
        ax1.set_xlim([0,365])
        plt.tight_layout()
       # plt.show()
        plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/pdfThroughTime.'+key+'.pdf', format='pdf')
    except:
        dum=1
        
    try:    
        fig7=plt.figure(num=None, figsize=(3.5,3.5), dpi=300, facecolor='w', edgecolor='w')
        plt.tick_params(axis='both', which='major', labelsize=6)
        plt.tick_params(axis='both', which='minor', labelsize=6)
        fig7.patch.set_alpha(0.0)
        ax7 = fig7.add_subplot(111)
        #ax7.set_aspect('equal', 'datalim')
        ax7.set_color_cycle(plt.cm.jet(np.linspace(0,1,nPlot+1)))
        legend=[]
        legend.append(1850)
        
        n, bins, patches = plt.hist(nSIF1850, bins, normed=1, histtype='step', linewidth=1, edgecolor='k' ,cumulative=True)
        for hItt in range(nPlot): # one per decade
            selnsif=nSIF[70+hItt*10,:]
            n, bins, patches = ax7.hist(selnsif, bins, normed=1, histtype='step',linewidth=0.5,cumulative=True)
            legend.append(1920+hItt*10)
        
        
        ax7.set_xlabel('Number of Sea Ice Free Days', fontsize=8)
        plt.legend(legend,fontsize=6)
        ax7.set_ylabel('Probability', fontsize=8)
        ax7.set_title('CDFs of open water at \n'+key+' (1920-2100)', fontsize=10)  
        ax7.set_ylim([0,1])
        ax7.set_xlim([0,365])  
        plt.tight_layout()
        plt.savefig('/Volumes/Pitcairn/seaicePPF/northernHemisphere/figures/cdfThroughTime.'+key+'.pdf', format='pdf')
    except:
        dum=1
    
        
    
    
    
    

## stack all sea ice concentrations to get probability of ice at that time. 
## will need to cut to the 1920-2100 section for this. 
#   
## show how pdfs of first/last, num, change through time.     
#     
#       
        
            
#fig1,ax1 = plt.subplts(1,2)
#ax1[0].set_color_cycle(plt.cm.spectral(np.linspace(0,1,1*len(dirList))))
#ax1[0].set_xlabel('Year')
#ax1[0].set_ylabel('Number of OW days')
#
#ax1[1].set_color_cycle(plt.cm.spectral(np.linspace(0,1,2*len(dirList))))
##ax1[1].set_xlabel('Year')
#ax1[1].set_ylabel('Day of Year')                   
#
#first=[]
#last=[]
#

#    
#    
#    
#    f.close()
#    del f
#    
#plt.show()
#

#
#fig3,ax3 = plt.subplts(1,1)
#ax3.boxplt(zip(*last))
#ax3.boxplt(zip(*first))
#plt.show()

