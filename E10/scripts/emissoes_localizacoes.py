# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 14:48:50 2024

@author: joseh
"""

import geopandas as gpd
import contextily as cx
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.geometry import Point
import pandas as pd
import ast

#%%


#%%

# Carregar o shapefile
shp = gpd.read_file("C:/Users/joseh/Downloads/Fontes_Difusas_pol_v00.shp")

shp.set_crs('EPSG:31982', inplace=True)
shp = shp.to_crs('EPSG:4674')

shp['geometry'] = shp['geometry'].translate(xoff=12.0)

shp_cong = gpd.read_file('C:\BolsaCongonhas\Git\Congonhas_LCQAr\shp\Shapefile_Congonhas.shp')

area_con = pd.read_csv(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\Con_Area_Geometry.csv')
volume_con = pd.read_csv(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\Con_Volume_Geometry.csv')

#area_con['geometry'] = area_con['geometry'].apply(ast.literal_eval)

#geo_area = gpd.GeoDataFrame(area_con, geometry = area_con.geometry.apply(Polygon))
#geo_volume = gpd.GeoDataFrame(volume_con, geometry = volume_con.geometry.apply(Polygon))

#%%

# Criar o plot
fig, ax = plt.subplots()

shp.plot(ax=ax, edgecolor='black', facecolor='none', alpha=1)
area_con.plot(ax=ax, edgecolor='blue', facecolor='none', alpha=1)
#volume_con.plot(ax=ax, edgecolor='red', facecolor='none', alpha=1)

shp_cong.boundary.plot(ax=ax, alpha=0.2)

# Adicionar o fundo de mapa com Contextily
cx.add_basemap(ax, crs=shp.crs, source=cx.providers.CartoDB.Positron, alpha=1)

plt.show()

#%%

vale_temp = pd.read_excel(r'C:\Users\joseh\Downloads\tx-recortelavraVALE-PM10PM25-05-11-2024.xls', sheet_name=2)

fig, ax = plt.subplots()

plt.plot(vale_temp.index,vale_temp['Emission Rate'])

plt.show()

#%%

csn_temp = pd.read_excel(r'C:\Users\joseh\Downloads\tx-Patio1CSN-PM10PM25-05-11-2024.xlsx', sheet_name=2)

fig, ax = plt.subplots()

plt.plot(csn_temp.index,csn_temp['Emission Rate'])

plt.show()

#%%

vale_temp2 = pd.read_excel(r'C:\Users\joseh\Downloads\tx-recortelavraVALE-PM10PM25-05-11-2024.xls', sheet_name=2)
vale_temp10 = pd.read_excel(r'C:\Users\joseh\Downloads\tx-recortelavraVALE-PM10PM25-05-11-2024.xls', sheet_name=1)

csn_temp2 = pd.read_excel(r'C:\Users\joseh\Downloads\tx-Patio1CSN-PM10PM25-05-11-2024.xlsx', sheet_name=2)
csn_temp10 = pd.read_excel(r'C:\Users\joseh\Downloads\tx-Patio1CSN-PM10PM25-05-11-2024.xlsx', sheet_name=1)

import scipy.stats

plt.scatter(csn_temp['Emission Rate'],vale_temp['Emission Rate'])

cor_csn = scipy.stats.spearmanr(csn_temp2['Emission Rate'],csn_temp10['Emission Rate'], nan_policy='omit')
cor_vale = scipy.stats.spearmanr(vale_temp2['Emission Rate'],vale_temp10['Emission Rate'], nan_policy='omit')





