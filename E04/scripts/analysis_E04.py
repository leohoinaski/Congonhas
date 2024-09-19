# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:06:38 2024

@author: Camilo Bastos Ribeiro
"""

import xarray as xr
import numpy as np
import pandas as pd
import temporalStatistics as tst
import os

def aggEmis(dir_folder, var, op, freq='monthly'):
    
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
            
            print(f"Processing sector: {sector_name}")
            
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
                    max_value = pol_2d.max()
                    print(f'max value of {var} in the tspep {t} = {max_value}')
                    
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
                print(f'max value of {var} = {max_value}')
                
                # Add the values to final df (yearly)
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
        
    # Add the total emissions column after identifying the major emitter
    if freq != 'yearly':
        for t, df in final_df.items():
            df['total_emissions'] = df.drop(columns=['major_emitter']).sum(axis=1)
    else:
        final_df['total_emissions'] = final_df.drop(columns=['major_emitter']).sum(axis=1)
    
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


def highTime(dfs, time_labels):
    
    """
    Identifies the time period with the highest emission for each pixel (row) and each sector (column).
    
    Parameters:
    - dfs: list of DataFrames, each corresponding to a different time period (e.g., months, days, hours),
           where each column represents an emission sector.
    - time_labels: list of strings representing the labels for each time period, in the same order as the DataFrames.
    
    Returns:
    - A DataFrame where each cell contains the label of the time period with the highest emission for the respective pixel and sector.
    """
   
    if len(dfs) != len(time_labels):
        raise ValueError("The number of DataFrames must match the number of time labels.")

    # Stack DataFrames along a new axis (adds a 'time' dimension)
    stacked_data = np.stack([df.to_numpy() for df in dfs], axis=2)
    
    # Find the index of the time period with the highest value for each pixel (row) and sector (column)
    max_indices = np.argmax(stacked_data, axis=2)

    # Create a DataFrame where each cell contains the label of the time period with the highest emission
    high_time_df = pd.DataFrame(
        data=np.vectorize(lambda idx: time_labels[idx])(max_indices),
        columns=dfs[0].columns,
        index=dfs[0].index
    )

    return high_time_df


def timeSeriesEmis(dir_folder, var, op, freq='monthly'):
   
    """
    Processes NetCDF files to compute the total emissions in the doamin 
    over a specified time frequency (monthly or daily).

    Parameters
    ----------
    dir_folder (str): Path to the folder containing the NetCDF files.
    var (str): Variable to be processed.
    op (function): Mathematical operation to be applied (e.g., np.sum, np.mean, np.median, etc.).
    freq (str): time frequency: 'monthly'or 'daily'.
        
    Returns
    -------
    final_df : pandas.DataFrame
        A DataFrame where each column represents the emissions for a particular sector, 
        and each row corresponds to a time period (months for 'monthly', days for 'daily'). 
        An additional 'total' column provides the summed emissions across all sectors for 
        each time period. 
    """
    
    final_df = pd.DataFrame()
    
    for file in os.listdir(dir_folder):
        
        if file.endswith('.nc'):
            sector_name = file.split('_')[0]
            
            print(f"Processing sector: {sector_name}")
            
            dir_data = os.path.join(dir_folder, file)
            data = xr.open_dataset(dir_data)
            
            tflag = data['TFLAG'].values
            time = pd.to_datetime(tflag, format='%Y%m%d%H')
            
            if freq == 'monthly':
                time_group = time.month
                unique_times = np.unique(time_group)
                time_index = range(1, 13)
           
            elif freq == 'daily':
                time_group = time.dayofyear
                unique_times = np.unique(time_group)
                time_index = range(1, 366)
            
            sector_sum = pd.DataFrame(0, index=time_index, columns=[sector_name])

            for t in unique_times:
                time_indices = np.where(time_group == t)[0]
                
                # Sum emission in the domain
                summed_emissions = op(np.array(data[var][time_indices, 0, :, :]), axis=0).sum()
                
                # Att the summed value to the sector
                sector_sum.loc[t, sector_name] = summed_emissions
            
            # Concat the sector in the final_df
            if final_df.empty:
                final_df = sector_sum
            else:
                final_df = pd.concat([final_df, sector_sum], axis=1)
    
    #Get total emission in the domain
    final_df['total'] = final_df.sum(axis=1)
    
    return final_df