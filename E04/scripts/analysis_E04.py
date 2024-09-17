# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:06:38 2024

@author: Camilo Bastos Ribeiro
"""

import xarray as xr
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import temporalStatistics as tst
import os

def emis_freq(dir_folder, var, op, freq='monthly'):
    
    """
    Function to aggregate emissions based on the specified mathematical operation and 
    temporal frequency.
    
    Parameters:
        dir_folder (str): Path to the folder containing the NetCDF files.
        var (str): Variable to be processed.
        op (function): Mathematical operation to be applied (e.g., np.sum, np.mean, np.median, etc.).
        freq (str): Aggregation frequency: 'monthly', 'weekly', 'hourly', or 'yearly'.
        
    Returns:
        dict or DataFrame: A dictionary where the keys are the months, days of the week, or hours of the day,
                           and the values are DataFrames with the operation applied to the emissions. 
                           If freq is 'yearly', it returns a single DataFrame.
    
    Dependencies:
    -------------
    - os: Interacting with the operating system (listing files in dir).
    - pandas (pd): Data manipulation and analysis.
    - numpy (np): Numerical operations.
    - xarray (xr): Working with NetCDF files.
    """
    
    final_df = {} if freq != 'yearly' else pd.DataFrame()
    
    for file in os.listdir(dir_folder):
        
        if file.endswith('.nc'):
            sector_name = file.split('_')[0]
            
            dir_data = os.path.join(dir_folder, file)
            data = xr.open_dataset(dir_data)
            
            tflag = data['TFLAG'].values
            time = pd.to_datetime(tflag, format='%Y%m%d%H')
            
            if freq == 'monthly':
                time_group = time.month
                unique_times = np.unique(time_group)
            
            elif freq == 'weekly':
                time_group = time.weekday
                unique_times = np.arange(7)
            
            elif freq == 'hourly':
                time_group = time.hour
                unique_times = np.arange(24)
            
            elif freq == 'yearly':
                time_group = None
            
            #Op for defined freq
            if freq != 'yearly':
                for t in unique_times:
                    time_indices = np.where(time_group == t)[0]
                    
                    #Op in the time inverval
                    pol_2d = pd.DataFrame(
                        op(np.array(data[var][time_indices, 0, :, :]), axis=0).flatten()
                    ).rename(columns={0: sector_name})
                    
                    #Add values for all sectors in the df
                    if t not in final_df:
                        final_df[t] = pol_2d
                    else:
                        final_df[t] = pd.concat([final_df[t], pol_2d], axis=1)
            
            else:
                # Op with emissions for entire year
                time_indices = range(len(time))
                pol_2d = pd.DataFrame(
                    op(np.array(data[var][time_indices, 0, :, :]), axis=0).flatten()
                ).rename(columns={0: sector_name})
                max_value = pol_2d.max()
                print(f'max value {var} = {max_value}')
                
                # Add the values to final df (yealy)
                if final_df.empty:
                    final_df = pol_2d
                else:
                    final_df = pd.concat([final_df, pol_2d], axis=1)
    
    # Get the major emitter
    if freq != 'yearly':
        for t, df in final_df.items():
            df['major_emitter'] = df.idxmax(axis=1)
    else:
        # Get the major emitter for 'yearly' freq
        final_df['major_emitter'] = final_df.idxmax(axis=1)
    
    return final_df

def latlon_2d(dir_data):
   
    """
    Function to convert latlon from NetCDF into 2D arrays.

    Parameters:
    ----------
    dir_data : str
        Path to the NetCDF file containing the data with the coordinates.

    Returns:
    -------
    lon2d : numpy.ndarray
        1D array containing the extracted and transformed longitude coordinates.
    
    lat2d : numpy.ndarray
        1D array containing the extracted and transformed latitude coordinates.

    Dependencies:
    -------------
    - xarray (xr)
    - tst (must contain the functions ioapiCoords and eqmerc2latlon)
    """

    #read data
    data = xr.open_dataset(dir_data)
    
    #processing coordinates
    xv, yv, lon, lat = tst.ioapiCoords(data)
    xlon, ylat = tst.eqmerc2latlon(data, xv, yv)
    lon2d = xlon.flatten()
    lat2d = ylat.flatten()
    
    return lon2d, lat2d
