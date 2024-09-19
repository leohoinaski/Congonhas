#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 09:59:08 2024

@author: leohoinaski
"""
import os 
import numpy as np
import netCDF4 as nc

def listDatasets(dataPath,GDNAM,pol):
    
    prefixed = [filename for filename in os.listdir(dataPath) if 
                filename.startswith("BRAIN_BASEMIS_"+GDNAM)]
    matching = [s for s in prefixed if pol in s]
    return matching

def all_equal(sequence):
    return len(set(sequence)) == 1
        
def checkMatEquals(monthlySum):
    if len(monthlySum.shape)==4:
        test = np.empty((monthlySum.shape[2],monthlySum.shape[3])).astype(bool)
        for ii in range(0,monthlySum.shape[2]):
            for jj in range(0,monthlySum.shape[3]):
                test[ii,jj] = all_equal(monthlySum[:,0,ii,jj])
                
    if len(monthlySum.shape)==3:
        test = np.empty((monthlySum.shape[1],monthlySum.shape[2])).astype(bool)
        for ii in range(0,monthlySum.shape[1]):
            for jj in range(0,monthlySum.shape[2]):
                test[ii,jj] = all_equal(monthlySum[:,ii,jj])
    return test

def agrmaxArray(data):
    if len(data.shape)==4:
        test = np.empty((data.shape[2],data.shape[3]))
        for ii in range(0,data.shape[2]):
            for jj in range(0,data.shape[3]):
                test[ii,jj] = np.nanargmax(data[:,0,ii,jj])
                
    if len(data.shape)==3:
        test = np.empty((data.shape[1],data.shape[2]))
        for ii in range(0,data.shape[1]):
            for jj in range(0,data.shape[2]):
                test[ii,jj] = np.nanargmax(data[:,ii,jj])
    return test

def createNETCDF(folderOut,name,data,ds,pollutant,xlon,ylat,datesTime):
    # print('===================STARTING netCDFcreator_v1.py=======================')
    # datesTime['TFLAG']=0
    # for ii in range(0,data.shape[0]):
    #     datesTime['TFLAG'][ii] = np.int32(str(datesTime.year[ii])+\
    #         str(datesTime.month[ii]).zfill(2)+\
    #             str(datesTime.day[ii]).zfill(2)+\
    #                 str(datesTime.hour[ii]).zfill(2))
          
    f2 = nc.Dataset(folderOut+'/'+name,'w') #'w' stands for write 
    for gatr in ds.ncattrs() :
        print(gatr)
        try:
            setattr(f2, gatr, ds.__getattribute__(gatr))
        except:
            print('bad var')
    f2.NVARS= data.shape[1]
    f2.HISTORY =''
    setattr(f2, 'VAR-LIST', pollutant)
    f2.NVARS= 1
    f2.NCOLS = data.shape[2]
    f2.NROWS = data.shape[1]
    #f2.NVARS = data.shape[1]
    f2.SDATE = str(datesTime.datetime[0])
    f2.FILEDESC = 'Concentration of ' +pollutant +' created by Leonardo Hoinaski - '
    f2.HISTORY = ''
    # # Specifying dimensions
    #tempgrp = f.createGroup('vehicularEmissions_data')
    f2.createDimension('TSTEP', None )
    f2.createDimension('DATE-TIME', 2)
    f2.createDimension('LAY', 1)
    f2.createDimension('VAR', 1)
    f2.createDimension('ROW', data.shape[1])
    f2.createDimension('COL', data.shape[2])
    # Building variables
    TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP'))
    # Passing data into variables
    TFLAG[:] = datesTime['datetime']
    LON = f2.createVariable('LON', 'f4', ( 'ROW','COL'))
    LAT = f2.createVariable('LAT', 'f4', ( 'ROW','COL'))
    LAT[:,:] =  ylat
    LON[:,:] = xlon
    LON.units = 'degrees '
    LAT.units = 'degrees '
    globals()[pollutant] = f2.createVariable(pollutant, np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    globals()[pollutant][:,0,:,:] = data[:,:,:]
    globals()[pollutant].units = ds[pollutant].units
    print('finishing')
    f2.close()
    return f2