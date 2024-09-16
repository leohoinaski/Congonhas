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
        

def averageTables(dataPath,pol,outPath):
    
    matching = listDatasets(dataPath,GDNAM,pol)
    
    print('Pollutant = ' + pol)
    
    for count, file in enumerate(matching):
        
        print('File = ' + file)
        source = file.split('_')[5]
        
        ds = nc.MFDataset(dataPath+'/'+file)
        xlon,ylat = ds['LON'][:], ds['LAT'][:]
        
       
        datesTime =  tst.datePrepBRAIN(ds)
        
        yearlySum,yearly = tst.yearlySum(datesTime,ds[pol])
        yearly[pol]=np.nansum(yearlySum)
        yearly.to_csv(outPath+'/'+pol+'_YEAR_'+source+'.csv', index=False)
        cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["azure","lightgray","crimson","darkred"])
        figMaker.spatialFigure(yearlySum[0,0,:,:],xlon,ylat,pol+' YEAR '+source,
                      cmap,borderShapePath,outPath,pol,IBGE_CODE,source)
        
        
        monthlySum,monthly = tst.monthlySum(datesTime,ds[pol])
        monthly[pol]=np.nansum(np.nansum(np.nansum(monthlySum,axis=1),axis=1),axis=1)
        monthly.to_csv(outPath+'/'+pol+'_MONTH_'+source+'.csv', index=False)
        maxmonth = monthlySum[:,0,:,:].argmax(axis=0)
        figMaker.maxPixelFigure(maxmonth,xlon,ylat,pol+' MaxMonth '+source,
                       cmap,borderShapePath,outPath,pol,IBGE_CODE,source,
                       'Month')

        hourlySum,hourly = tst.hourlySum(datesTime,ds[pol])
        hourly[pol]=np.nansum(np.nansum(np.nansum(hourlySum,axis=1),axis=1),axis=1)
        hourly.to_csv(outPath+'/'+pol+'_HOUR_'+source+'.csv', index=False)
        maxhour= hourlySum[:,0,:,:].argmax(axis=0)
        figMaker.maxPixelFigure(maxhour,xlon,ylat,pol+' MaxHour '+source,
                                cmap,borderShapePath,outPath,pol,IBGE_CODE,source,
                                'Hour')

        dayOfWeekSum,dayOfWeek = tst.dayOfWeekSum(datesTime,ds[pol])
        dayOfWeek[pol]=np.nansum(np.nansum(np.nansum(dayOfWeekSum,axis=1),axis=1),axis=1)
        dayOfWeek.to_csv(outPath+'/'+pol+'_DAYofWEEK_'+source+'.csv', index=False)
        maxdayOfWeek= dayOfWeekSum[:,0,:,:].argmax(axis=0)
        figMaker.maxPixelFigure(maxdayOfWeek,xlon,ylat,pol+' MarDayOfWeek '+source,
                                cmap,borderShapePath,outPath,pol,IBGE_CODE,source,
                                'DAYofWEEK')
        
        if count == 1:
            dataTotal = ds[pol]
        else:
            dataTotal = dataTotal + ds[pol]

