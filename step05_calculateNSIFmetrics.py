# calculate mean, median, std of all 

print 'importing modules'
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import glob
import nio
from scipy import stats

#import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['font.sans-serif']=['Arial']

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
	print 'opening background: ', nsk
	for fn in dirList2:
		f=Dataset(fn, 'r')
		numSIF1850.extend(f.variables['numSIF'][:,:,:])
		fir1850.extend(f.variables['firstDOY'][:,:,:])
		las1850.extend(f.variables['lastDOY'][:,:,:])
		print fn, (datetime.now()-startTime) 
		f.close()   

	nSIF1850=np.array(numSIF1850)
	first1850=np.array(fir1850)
	last1850=np.array(las1850)
	del numSIF1850
	del fir1850
	del las1850
		
	for ittR in range(len(rcpkey)):
		fnAnalysis='/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/allModels.Summary.'+nsk+ '.'+ rcpName[ittR] +'.nc'
		print fnAnalysis
		dirList=glob.glob(path+'*'+rcpName[ittR]+'*'+nsk+'*Analysis.nc')
	
		analysisFiles=[]    
		for fn in dirList:
			analysisFiles.append(Dataset(fn, 'r'))
			#print fn, (datetime.now()-startTime)
		

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
		del f
	

		numSIF=np.nan*np.ones((numModels,numYears,nj,ni))
		first=np.nan*np.ones((numModels,numYears,nj,ni))
		last=np.nan*np.ones((numModels,numYears,nj,ni))

		i=0
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


		origData={'nSIF1850':nSIF1850,
		'nSIF':nSIF,
		'first1850':first1850,
		'last1850':last1850,
		'first':first,
		'last':last}

		dataset={}
		for k in origData.keys():
			d=origData[k]
			# do appropriate masking:
			d[d>900]=np.nan # set days with all ice to nan for appropriate

			# statistics
			# if not all models have values, we don't want to make a statistic, therefore don't use nan
			# changed mind. do want statistic
			tempMean=np.nanmean(d, axis=0)
			tempMedian=np.median(d, axis=0)
			tempStd=np.nanstd(d, axis=0) 

			if len(tempMean.shape)==2:
				tempMean[landMask]=fillVal   
				tempMedian[landMask]=fillVal
				tempStd[landMask]=fillVal 
			else:
				tempMean[:,landMask]=fillVal   
				tempMedian[:,landMask]=fillVal
				tempStd[:,landMask]=fillVal 
	
			dataset[k+'Mean']={'data':tempMean}
			dataset[k+'Median']={'data':tempMedian}
			dataset[k+'Std']={'data':tempStd}

			if k[0]=='n':
				dataset[k+'Std']['units']='number of days'
				dataset[k+'Mean']['units']='number of days'    
				dataset[k+'Median']['units']='number of days'
	
			else:
				dataset[k+'Std']['units']='day of year'
				dataset[k+'Mean']['units']='day of year'    
				dataset[k+'Median']['units']='day of year'
			
			dataset[k+'Std']['longName']='Standard deviation of the number of '+k
			dataset[k+'Mean']['longName']='Mean number of '+k    
			dataset[k+'Median']['longName']='Median number of '+k

		print k, (datetime.now()-startTime)

		# calculate comparisons to 1850 (once that modle has been analyzed)
		for k in ['nSIF', 'first', 'last']:
			for kadd in ['Std', 'Median', 'Mean']:
				key=k+kadd
				key2=k+'1850'+kadd
	
				data=dataset[key]['data']
				data2=dataset[key2]['data']
	
				compdata=np.ones(data.shape)
				for i in range(data.shape[0]):
					compdata[i,:,:]=(data[i,:,:]-data2)/data2
				compdata[:,landMask]=fillVal
				isnan=np.isnan(compdata)
				compdata[isnan]=fillVal  
	
				isinf=np.isinf(compdata)
				compdata[isinf]=fillVal
	
				dataset[key+'Comp_prop']={'data':compdata,
				'units':'Difference with 1850 normalized by 1850 value',
				'longName': 'Comparison of '+dataset[key]['longName']+'with 1850'}
		
				compdata2=np.ones(data.shape)
				for i in range(data.shape[0]):
					compdata2[i,:,:]=data[i,:,:]-data2
				compdata2[:,landMask]=fillVal
				dataset[key+'Comp_abs']={'data':compdata2,
				'units':'Difference with 1850',
				'longName': 'Comparison of '+dataset[key]['longName']+'with 1850'}

		######### calculate slopes                   
						 
		slope=np.nan*np.ones((nm, nj, ni))
		slope_masked=np.nan*np.ones((nm, nj, ni))
		R=np.nan*np.ones((nm, nj, ni))
		pval=np.nan*np.ones((nm, nj, ni)) 
		icemask=np.nanmean(nSIF1850, axis=0)<350

		final_year=np.nan*np.ones((nm, nj, ni))

		time=np.arange(1850,lastYear+2)
		startInd=129 # (1979-2014)
		stopInd=165
		for nii in range(ni):
			print 'calculating trendlines', str(nii)+'/'+str(ni), (datetime.now()-startTime)
			for njj in range(nj):
				if  (icemask[njj, nii]==True): 
					for nmm in range(nm):             
								
						s, icpts, r, p, ster = stats.linregress(time[startInd:stopInd],nSIF[nmm,startInd:stopInd, njj, nii])               
						slope[nmm, njj, nii]=s
						R[nmm, njj, nii]=r
						pval[nmm, njj, nii]=p 
						if p<0.05:
							slope_masked[nmm, njj, nii]=s
						#### choose last year that there is are more than 15 days of ice per year.
						fyind=np.where(nSIF[nmm,:, njj, nii]>182)[0]
						if len(fyind)==0:
							final_year[nmm, njj, nii]=np.nan
						else:
							final_year[nmm, njj, nii]=time[fyind[0]]

				
		final_year_mean=np.nanmean(final_year, axis=0)

		isnan=np.isnan(final_year_mean)
		final_year_mean[isnan]=fillVal  
		isinf=np.isinf(final_year_mean)
		final_year_mean[isinf]=fillVal   

		isnan=np.isnan(final_year)
		final_year[isnan]=fillVal  
		isinf=np.isinf(final_year)
		final_year[isinf]=fillVal  

		# get rid of infs before averaging.  
		isinf=np.isinf(slope_masked)
		slope_masked[isinf]=np.nan 
	
		isinf=np.isinf(slope)
		slope[isinf]=np.nan                                                                                                                                                                           
																																																																																																																																	   
		meanSlope=np.nanmean(slope, axis=0)
		stdSlope=np.nanstd(slope, axis=0)  
																		 
		meanSlope_masked=np.nanmean(slope_masked, axis=0)
		stdSlope_masked=np.nanstd(slope_masked, axis=0)                                                    
	
		isnan=np.isnan(slope_masked)
		slope_masked[isnan]=fillVal  
		isinf=np.isinf(slope_masked)
		slope_masked[isinf]=fillVal                                                                                                                                                      
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																		   
		isnan=np.isnan(slope)
		slope[isnan]=fillVal  
		isinf=np.isinf(slope)
		slope[isinf]=fillVal                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
																																																																																																																																																																																																																																																								   
		isnan=np.isnan(meanSlope_masked)
		meanSlope_masked[isnan]=fillVal  
		isinf=np.isinf(meanSlope_masked)
		meanSlope_masked[isinf]=fillVal  

		isnan=np.isnan(stdSlope_masked)
		stdSlope_masked[isnan]=fillVal  
		isinf=np.isinf(stdSlope_masked)
		stdSlope_masked[isinf]=fillVal                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																			
		# 2.5 Save this to a netCDF

		f=nio.open_file(dirList[0], 'r')

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
			fAnVars['time'][:]=f.variables['time'][-numYears:]
			fAnVars['time_bounds'][:,:]=f.variables['time_bounds'][-numYears:,:]

		dataKey=dataset.keys()

		stdAttributes={'_FillValue': np.array([  1.00000002e+30], dtype=float),
		'cell_measures': 'area: tarea',
		'cell_methods': 'time: mean',
		'comment': 'none',
		'coordinates': 'TLON TLAT time',
		'missing_value': np.array([  1.00000002e+30], dtype=float),
		'time_rep': 'averaged'}

		for i in range(len(dataKey)):
			key=dataKey[i]
			d=dataset[key]
			data=d['data']
			print 'creating ', key
			# the netCDF4 module requires that if a fill value exists, it must be declared when the variable is created.
			
			isnan=np.isnan(data)
			data[isnan]=fillVal  
			isinf=np.isinf(data)
			data[isinf]=fillVal   
			
			if len(data.shape)==3:
				fAnVars[key]=fAn.createVariable(key, 'f',  ('time', 'nj', 'ni'), fill_value=stdAttributes['_FillValue'])
				fAnVars[key][:,:,:]=data
				('time', 'nj', 'ni')
			elif len(data.shape)==2:
				fAnVars[key]=fAn.createVariable(key, 'f',  ('nj', 'ni'), fill_value=stdAttributes['_FillValue'])
				fAnVars[key][:,:]=data
		
			# set all the attribute keys.
			for attKey in stdAttributes.keys():
				if attKey != '_FillValue':
					setattr(fAn.variables[key],attKey,stdAttributes[attKey])
			setattr(fAn.variables[key],'long_name',d['longName'])
			setattr(fAn.variables[key],'units',d['units'])

		fAn.createDimension('nm', nm)

		cKey='meanSlope'
		fAn.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
		setattr(fAn.variables[cKey],'long_name','Mean of Slope (1979-2014)') 
		setattr(fAn.variables[cKey],'units','year')   
		setattr(fAn.variables[cKey], 'coordinates', 'TLON TLAT')

		cKey='stdSlope'
		fAn.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
		setattr(fAn.variables[cKey],'long_name','Standard Deviation of Slope (1979-2014)') 
		setattr(fAn.variables[cKey],'units','year')   
		setattr(fAn.variables[cKey], 'coordinates', 'TLON TLAT')

																						 
		cKey='meanSlope_masked'
		fAn.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
		setattr(fAn.variables[cKey],'long_name','Mean of Slope (masked by Pval) (1979-2014)') 
		setattr(fAn.variables[cKey],'units','year')   
		setattr(fAn.variables[cKey], 'coordinates', 'TLON TLAT')

		cKey='stdSlope_masked'
		fAn.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
		setattr(fAn.variables[cKey],'long_name','Standard Deviation of Slope (masked by Pval) (1979-2014)') 
		setattr(fAn.variables[cKey],'units','year')   
		setattr(fAn.variables[cKey], 'coordinates', 'TLON TLAT')
																																																											  
																																																																																											 
		cKey='pval'
		fAn.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
		setattr(fAn.variables[cKey],'long_name','Trend pValue') 
		setattr(fAn.variables[cKey],'units','')   
		setattr(fAn.variables[cKey], 'coordinates', 'nm TLON TLAT')

		cKey='R'
		fAn.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
		setattr(fAn.variables[cKey],'long_name','R2') 
		setattr(fAn.variables[cKey],'units','')
		setattr(fAn.variables[cKey], 'coordinates', 'nm TLON TLAT')

		cKey='slope'
		fAn.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
		setattr(fAn.variables[cKey],'long_name','Trend Slope') 
		setattr(fAn.variables[cKey],'units','Days per Year')
		setattr(fAn.variables[cKey], 'coordinates', 'nm TLON TLAT')

		cKey='slope_masked'
		fAn.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
		setattr(fAn.variables[cKey],'long_name','Trend Slope Masked by Pval <0.05 (1979-2014)') 
		setattr(fAn.variables[cKey],'units','Days per Year')
		setattr(fAn.variables[cKey], 'coordinates', 'nm TLON TLAT')

		cKey='final_year_mean'
		fAn.createVariable(cKey, 'f', ('nj','ni'),fill_value=fillVal)
		setattr(fAn.variables[cKey],'long_name','Mean of Final Year with nSIF<182') 
		setattr(fAn.variables[cKey],'units','year')   
		setattr(fAn.variables[cKey], 'coordinates', 'TLON TLAT')

		cKey='final_year'
		fAn.createVariable(cKey, 'f', ('nm', 'nj','ni'),fill_value=fillVal)
		setattr(fAn.variables[cKey],'long_name','Final Year with nSIF<182') 
		setattr(fAn.variables[cKey],'units','year')
		setattr(fAn.variables[cKey], 'coordinates', 'nm TLON TLAT')

		fAn.variables['final_year'][:,:,:]=final_year
		fAn.variables['final_year_mean'][:,:]=final_year_mean

		fAn.variables['slope_masked'][:,:,:]=slope_masked
		fAn.variables['slope'][:,:,:]=slope
		fAn.variables['meanSlope'][:,:]=meanSlope
		fAn.variables['stdSlope'][:,:]=stdSlope

		fAn.variables['meanSlope_masked'][:,:]=meanSlope_masked
		fAn.variables['stdSlope_masked'][:,:]=stdSlope_masked

		fAn.variables['pval'][:,:,:]=pval
		fAn.variables['R'][:,:,:]=R      
				
		fAn.close()    
		f.close()
		del f
		del fAn
		print'finished analysis of ',fn , (datetime.now()-startTime)