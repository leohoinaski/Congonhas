#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 13:39:11 2024



@author: leohoinaski
"""



import os
import netCDF4 as nc
import temporalStatistics as tst
import numpy as np
import pandas as pd
import figureMaker as figMaker
import matplotlib as mpl
import ismember

NO2 = {
  "Pollutant": "$NO_{2}$",
  "Unit": '$\u03BCg.m^{-3}$',
  "conv": 1880,
  "tag":'NO2',
  #"Criteria": 260, # 260, 240, 220, 200
}

CO = {
  "Pollutant": "CO",
  "Unit": 'ppb',
  "conv": 1000, # Conversão de ppm para ppb
  "tag":'CO',
}

O3 = {
  "Pollutant": "$O_{3}$",
  "Unit": 'ppm',
  "conv": 1,
  "tag":'O3'
}

SO2 = {
  "Pollutant": "$SO_{2}$",
  "Unit": '$\u03BCg.m^{-3}$',
  "conv": 2620,
  "tag":'SO2'
}

PM10 = {
  "Pollutant": "$PM_{10}$",
  "Unit": '$\u03BCg.m^{-3}$',
  "conv": 1,
  "tag":'PM10',
}

PM25 = {
  "Pollutant": "$PM_{2.5}$",
  "Unit": '$\u03BCg.m^{-3}$',
  "conv": 1,
  "tag":'PM25',
}


pollutants = [CO]

pol = 'PM10'

GDNAM = 'Con_3km'

IBGE_CODE = 3118007

cwdPath = os.path.abspath(os.getcwd())

outPath = os.path.dirname(os.path.abspath(os.getcwd()))+'/outputs'
os.makedirs(outPath,exist_ok=True)

rootPath = os.path.dirname(os.path.dirname(cwdPath))

dataPath = rootPath + '/data/'+GDNAM

borderShapePath = os.path.dirname(rootPath)+'/shapefiles/BR_Municipios_2020.shp'

def listDatasets(dataPath,GDNAM,pol):
    
    prefixed = [filename for filename in os.listdir(dataPath) if 
                filename.startswith("BRAIN_BASEMIS_"+GDNAM)]
    matching = [s for s in prefixed if pol in s]
    return matching
        



matching = listDatasets(dataPath,GDNAM,pol)

print('Pollutant = ' + pol)

for count, file in enumerate(matching):
    
    print('File = ' + file)
    source = file.split('_')[5]
    
    ds = nc.MFDataset(dataPath+'/'+file)
    xlon,ylat = ds['LON'][:], ds['LAT'][:]
    
   
    datesTime =  tst.datePrepBRAIN(ds)
    
    if count==0:
        datePfct = pd.date_range(start=str(datesTime.year[0])+'-01-01 00:00:00',
                                 end=str(datesTime.year[datesTime.shape[0]-1])+'-12-31 23:00:00', 
                                 freq='h')
        dates = pd.DatetimeIndex(datePfct)
        datePfct=pd.DataFrame()
        datePfct['year'] = dates.year
        datePfct['month'] = dates.month
        datePfct['day'] = dates.day
        datePfct['hour'] = dates.hour
        datePfct['day_of_week'] = dates.day_name()
        datePfct['datetime']=dates
        
    
    lia,loc = ismember.ismember(datePfct.datetime.values,datesTime.datetime.values)
    
    data = ds[pol][loc,:,:]
    
    
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["azure","lightgray","crimson","darkred"])
    figMaker.cityEmissTimeSeries(data,xlon,ylat,datesTime,borderShapePath,
                                 outPath,pol,str(IBGE_CODE),source,cmap)
    
    
    # Análise por ano
    yearlySum,yearly = tst.yearlySum(datesTime.iloc[loc,:],data)
    yearly[pol]=np.nansum(yearlySum)
    yearly.to_csv(outPath+'/'+pol+'_YEAR_'+source+'.csv', index=False)
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["azure","lightgray","crimson","darkred"])
    figMaker.spatialFigure(yearlySum[0,0,:,:],xlon,ylat,pol+' YEAR '+source,
                  cmap,borderShapePath,outPath,pol,IBGE_CODE,source)
    
    
    # Análise por mês
    monthlySum,monthly = tst.monthlySum(datesTime.iloc[loc,:],data)
    monthly[pol]=np.nansum(np.nansum(np.nansum(monthlySum,axis=1),axis=1),axis=1)
    monthly.reset_index().to_csv(outPath+'/'+pol+'_MONTH_'+source+'.csv', index=False)
    monthlySum[monthlySum==0]=np.nan
    maxmonth = monthlySum[:,0,:,:].argmax(axis=0).astype(float)
    maxmonth[np.isnan(monthlySum).all(axis=0)[0,:,:]] = np.nan
    maxmonth[(monthlySum==monthlySum).all(axis=0)[0,:,:]] = 99
    figMaker.maxPixelFigure(maxmonth,xlon,ylat,pol+' MaxMonth '+source,
                   cmap,borderShapePath,outPath,pol,IBGE_CODE,source,
                   'Month')
    
    # Análise por hora
    hourlySum,hourly = tst.hourlySum(datesTime.iloc[loc,:],data)
    hourly[pol]=np.nansum(np.nansum(np.nansum(hourlySum,axis=1),axis=1),axis=1)
    hourly.reset_index().to_csv(outPath+'/'+pol+'_HOUR_'+source+'.csv', index=False)
    maxhour= hourlySum[:,0,:,:].argmax(axis=0).astype(float)
    hourlySum[hourlySum==0]=np.nan
    maxhour[np.isnan(hourlySum).all(axis=0)[0,:,:]] = np.nan
    maxhour[(hourlySum==hourlySum).all(axis=0)[0,:,:]] = np.nan
    figMaker.maxPixelFigure(maxhour,xlon,ylat,pol+' MaxHour '+source,
                            cmap,borderShapePath,outPath,pol,IBGE_CODE,source,
                            'Hour')
    
    # Análise por dia da semana
    dayOfWeekSum,dayOfWeek = tst.dayOfWeekSum(datesTime.iloc[loc,:],data)
    dayOfWeek[pol]=np.nansum(np.nansum(np.nansum(dayOfWeekSum,axis=1),axis=1),axis=1)
    dayOfWeek.reset_index().to_csv(outPath+'/'+pol+'_DAYofWEEK_'+source+'.csv', index=False)
    maxdayOfWeek= dayOfWeekSum[:,0,:,:].argmax(axis=0).astype(float)
    dayOfWeekSum[dayOfWeekSum==0]=np.nan
    maxdayOfWeek[np.isnan(dayOfWeekSum).all(axis=0)[0,:,:]] = np.nan
    maxdayOfWeek[(dayOfWeekSum==dayOfWeekSum).all(axis=0)[0,:,:]] = np.nan
    figMaker.maxPixelFigure(maxdayOfWeek,xlon,ylat,pol+' MaxDayOfWeek '+source,
                            cmap,borderShapePath,outPath,pol,IBGE_CODE,source,
                            'DayOfWeek')

    # Acumulando em uma matriz com todas as fontes
    if count == 0:
        
        
        
        dataTotal = np.empty([datePfct.shape[0], 
                                    xlon.shape[0],xlon.shape[1]])
        
        dataTotal[lia,:,:] = ds[pol][loc,0,:,:]
        
        dataBySourceYear = np.empty([len(matching), 
                                    dataTotal.shape[1],dataTotal.shape[2]])
        dataBySourceYear[count,:,:] = yearlySum[0,0,:,:]
        
        dataBySourceMonth = np.empty([len(matching),12,
                                    dataTotal.shape[1],dataTotal.shape[2]])
        dataBySourceMonth[count,:,:,:] = monthlySum[:,0,:,:]
        
        dataBySourceHour = np.empty([len(matching), 24,
                                    dataTotal.shape[1],dataTotal.shape[2]])
        dataBySourceHour[count,:,:,:] = hourlySum[:,0,:,:]
        
        dataBySourcedayOfWeek= np.empty([len(matching), 7,
                                    dataTotal.shape[1],dataTotal.shape[2]])
        dataBySourcedayOfWeek[count,:,:,:] = dayOfWeekSum[:,0,:,:]
        
        
    else:
        dataTotalNew = np.empty([datePfct.shape[0], 
                                    xlon.shape[0],xlon.shape[1]])
        
        lia,loc = ismember.ismember(datePfct.datetime.values,datesTime.datetime.values)
        
        dataTotalNew[lia,:,:] = ds[pol][loc,0,:,:]
        
        dataTotal = dataTotal + dataTotalNew
        
        dataBySourceYear[count,:,:] = yearlySum[0,0,:,:]
        dataBySourceMonth[count,:,:,:] = monthlySum[:,0,:,:]
        dataBySourceHour[count,:,:,:] = hourlySum[:,0,:,:]
        dataBySourcedayOfWeek[count,:,:,:] = dayOfWeekSum[:,0,:,:]

# Análise da soma das emissões
source = 'TOTAL'
yearlySum,yearly = tst.yearlySum(datePfct,dataTotal)
yearly[pol]=np.nansum(yearlySum)
yearly.to_csv(outPath+'/'+pol+'_YEAR_'+source+'.csv', index=False)
cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["azure","lightgray","crimson","darkred"])
figMaker.spatialFigure(yearlySum[0,0,:,:],xlon,ylat,pol+' YEAR '+source,
              cmap,borderShapePath,outPath,pol,IBGE_CODE,source)


monthlySum,monthly = tst.monthlySum(datePfct,dataTotal)
monthly[pol]=np.nansum(np.nansum(np.nansum(monthlySum,axis=1),axis=1),axis=1)
monthly.to_csv(outPath+'/'+pol+'_MONTH_'+source+'.csv', index=False)
maxmonth = monthlySum[:,0,:,:].argmax(axis=0)
figMaker.maxPixelFigure(maxmonth,xlon,ylat,pol+' MaxMonth '+source,
               cmap,borderShapePath,outPath,pol,IBGE_CODE,source,
               'Month')

hourlySum,hourly = tst.hourlySum(datePfct,dataTotal)
hourly[pol]=np.nansum(np.nansum(np.nansum(hourlySum,axis=1),axis=1),axis=1)
hourly.to_csv(outPath+'/'+pol+'_HOUR_'+source+'.csv', index=False)
maxhour= hourlySum[:,0,:,:].argmax(axis=0)
figMaker.maxPixelFigure(maxhour,xlon,ylat,pol+' MaxHour '+source,
                        cmap,borderShapePath,outPath,pol,IBGE_CODE,source,
                        'Hour')

dayOfWeekSum,dayOfWeek = tst.dayOfWeekSum(datePfct,dataTotal)
dayOfWeek[pol]=np.nansum(np.nansum(np.nansum(dayOfWeekSum,axis=1),axis=1),axis=1)
dayOfWeek.to_csv(outPath+'/'+pol+'_DAYofWEEK_'+source+'.csv', index=False)
maxdayOfWeek= dayOfWeekSum[:,0,:,:].argmax(axis=0)
figMaker.maxPixelFigure(maxdayOfWeek,xlon,ylat,pol+' MarDayOfWeek '+source,
                        cmap,borderShapePath,outPath,pol,IBGE_CODE,source,
                        'DAYofWEEK')
