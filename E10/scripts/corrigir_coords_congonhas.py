# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 15:54:27 2024

@author: joseh
"""

#%%

import numpy as np
import pandas as pd
import xarray as xr
import pyproj
from pyproj import Proj, transform
import shapely
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import Polygon
import geopandas as gpd
import ast


#%% Puxar planilha

pd.set_option('display.float_format', '{:.15f}'.format)  # Isso garante que até 15 casas decimais sejam exibidas

br_coord = pd.read_csv(r'C:/BolsaCongonhas/Git/Congonhas/data/Planilhas/BR_Coords.csv')

#%% Função para converter UTM para latitude e longitude
def utm_to_latlong(row, x_col, y_col):
    X = row[x_col]
    Y = row[y_col]
    longitude, latitude = transform(utm_proj, wgs84_proj, X, Y)
    return pd.Series({'longitude': longitude, 'latitude': latitude})

#%% Mudar coordenadas

# Definir o sistema de coordenadas UTM (exemplo: zona 23S para o hemisfério sul)
utm_proj = Proj(proj='utm', zone=23, south=True, ellps='WGS84')

# Definir o sistema de coordenadas geográficas (WGS84)
wgs84_proj = Proj(proj='latlong', datum='WGS84')

# Aplicar a conversão para cada par de colunas Xn e Yn
for i in range(1, 58):
    br_coord[[f'X{i}', f'Y{i}']] = br_coord.apply(utm_to_latlong, axis=1, x_col=f'X{i}', y_col=f'Y{i}')
    
#%% Salvar em csv e excel

br_coord.to_csv(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\coordenadas_corrigidas.csv',index=False)
#br_coord.to_excel(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\coordenadas_corrigidas.xlsx',index=False)

#%% Fazer geometria

def juntar_lat_lon(row):
    # Junta as colunas de latitude e longitude em uma lista de pares
    return [[row[f'X{i}'], row[f'Y{i}']] 
            for i in range(1, 58)
            if not pd.isna(row[f'X{i}']) and not pd.isna(row[f'Y{i}'])]

# Cria uma nova coluna com as listas de pares de coordenadas
br_coord['geometry'] = br_coord.apply(juntar_lat_lon, axis=1)

# Mantém apenas as colunas 'Tipo', 'ID' e a nova coluna 'Coordenadas'
br_coord_geometrias = br_coord[['Type', 'ID', 'geometry']]

#%% Encontrar centroide

from shapely.geometry import Polygon, Point, LineString

# Função para calcular o centróide
def calcular_centroide(coords):
    # Verifica se temos pelo menos 3 coordenadas para formar um polígono, caso contrário usa LineString
    if len(coords) >= 3:
        poligono = Polygon(coords)
        return poligono.centroid.coords[0]
    elif len(coords) >= 2:
        linha = LineString(coords)
        return linha.centroid.coords[0]
    else:
        return coords[0] if coords else (None, None)  # Retorna o próprio ponto se só houver uma coordenada

# Aplica a função para cada conjunto de coordenadas
br_coord_geometrias['Centroide'] = br_coord_geometrias['geometry'].apply(calcular_centroide)
br_coord_geometrias[['X', 'Y']] = pd.DataFrame(br_coord_geometrias['Centroide'].tolist(), index=br_coord_geometrias.index)


#%%

br_coord_geometrias.to_csv(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\BR_Geometry.csv',index=False)
#br_coord_geometrias.to_excel(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\BR_Geometry.xlsx',index=False)

#%%

con_area = pd.read_excel(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\Con_AREAS_GEOMETRYS.xlsx')

#%%
    
def juntar_lat_lon(row):
    # Junta as colunas de latitude e longitude em uma lista de pares
    return [[row[f'X{i}'], row[f'Y{i}']] 
            for i in range(1, 5)
            if not pd.isna(row[f'X{i}']) and not pd.isna(row[f'Y{i}'])]

con_area[['X2', 'Y2']] = con_area[['X1', 'Y1']].values + con_area[['Length_X', 'Rotation_Angle']].values
con_area[['X3', 'Y3']] = con_area[['X1', 'Y1']].values + con_area[['Length_X', 'Length_Y']].values
con_area[['X4', 'Y4']] = con_area[['X1', 'Y1']].values + con_area[['Rotation_Angle', 'Length_Y']].values

# Definir o sistema de coordenadas UTM (exemplo: zona 23S para o hemisfério sul)
utm_proj = Proj(proj='utm', zone=23, south=True, ellps='WGS84')

# Definir o sistema de coordenadas geográficas (WGS84)
wgs84_proj = Proj(proj='latlong', datum='WGS84')

# Aplicar a conversão para cada par de colunas Xn e Yn
for i in range(1, 5):
    con_area[[f'X{i}', f'Y{i}']] = con_area.apply(utm_to_latlong, axis=1, x_col=f'X{i}', y_col=f'Y{i}')

# Cria uma nova coluna com as listas de pares de coordenadas
con_area['geometry'] = con_area.apply(juntar_lat_lon, axis=1)

# Mantém apenas as colunas 'Tipo', 'ID' e a nova coluna 'Coordenadas'
con_area_geometrias = con_area[['ID', 'geometry']]

#%% Salvar Areas

con_area_geometrias.to_csv(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\Con_Area_Geometry.csv',index=False)

#%% Criar geometria volume

con_volume = pd.read_csv(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\Con_VOLUME_Geometry.csv')

def geometria_volume(lon,lat,syinit):
    
    series = gpd.GeoSeries(
        Point(lon, lat)
        )
    
    geo_series = series.buffer(syinit/2)

    return geo_series

lista_geometrias = []

for i in range(len(con_volume)):
     
    lon = con_volume['X1'][i]
    lat = con_volume['Y1'][i]
    syinit = con_volume['SigmaY'][i]
    
    geometry = geometria_volume(lon,lat,syinit)
    geometry = geometry.set_crs(epsg = 32723)
    geometry = geometry.to_crs(epsg=4326)
    
    lista_geometrias.append(geometry)
    
#%%

geometria = pd.concat(lista_geometrias).reset_index(drop=True)

#%%

con_volume = pd.read_csv(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\CongonhasSources\Con_VOLUME.csv', encoding = 'latin')

geo_volume = gpd.GeoDataFrame(con_volume, geometry = geometria)

#%% Transformar g/s de Line-Volume pra g/sm²

con_linevolume = pd.read_excel(r"C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\Congonhas_Line_g_s.xlsx")

# Função para garantir que o valor é uma lista e não uma string
def convert_to_list(geom):
    if isinstance(geom, str):  # Verifica se o valor está em formato de string
        return ast.literal_eval(geom)  # Converte a string em lista
    return geom  # Caso já seja lista, retorna como está

# Aplica a função para garantir que todas as geometrias estão em formato de lista
con_linevolume['geometry'] = con_linevolume['geometry'].apply(convert_to_list)

# Converter as listas de coordenadas em geometria de linhas (LineString)
con_linevolume['geometry'] = con_linevolume['geometry'].apply(lambda x: Polygon(x) if len(x) >= 3 else None)

# Converter o DataFrame em GeoDataFrame
geo_linevolume = gpd.GeoDataFrame(con_linevolume, geometry='geometry', crs="EPSG:4326")

# Transformar para UTM Zone 23
geo_linevolume = geo_linevolume.to_crs("EPSG:32623")  # Use EPSG:32723 para o hemisfério sul

# Calcular o comprimento das linhas
geo_linevolume['length_m'] = geo_linevolume['geometry'].area

# Selecionar as colunas que começam com 'Ln'
ln_columns = [col for col in geo_linevolume.columns if col.startswith('Ln')]

# Multiplicar essas colunas pelas colunas 'X' e 'Y'
for col in ln_columns:
    geo_linevolume[col] = geo_linevolume[col] / ( geo_linevolume['length_m'] )

geo_linevolume.to_csv(r'C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\Con_Area_area.csv',index=False)

#%% Transformar em shapefiles com a tabela de atributos da planilha

# Função para garantir que o valor é uma lista e não uma string
def convert_to_list(geom):
    if isinstance(geom, str):  # Verifica se o valor está em formato de string
        return ast.literal_eval(geom)  # Converte a string em lista
    return geom  # Caso já seja lista, retorna como está

con_point = pd.read_csv(r"C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\EDAs_Congonhas\Con_POINTS.csv", encoding = 'latin')
con_line = pd.read_csv(r"C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\EDAs_Congonhas\Con_LINE.csv", encoding = 'latin')
con_area = pd.read_csv(r"C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\EDAs_Congonhas\Con_AREA.csv", encoding = 'latin')
con_volume = pd.read_csv(r"C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\EDAs_Congonhas\Con_VOLUME.csv", encoding = 'latin')

# Aplica a função para garantir que todas as geometrias estão em formato de lista
con_line['geometria'] = con_line['geometry'].apply(convert_to_list)
con_area['geometria'] = con_area['geometry'].apply(convert_to_list)

# Converter as listas de coordenadas em geometrias points, lines e polygons
con_point['geometria'] = con_point.apply(lambda row: Point(row['Longitude'], row['Latitude']), axis=1)
con_line['geometria'] = con_line['geometria'].apply(LineString)
con_area['geometria'] = con_area['geometria'].apply(lambda x: Polygon(x) if len(x) >= 3 else None)
con_volume['geometria'] = con_volume.apply(lambda row: Point(row['Longitude'], row['Latitude']), axis=1)

# Suponha que seu DataFrame seja df e a coluna de geometria seja 'geometry'
con_point = gpd.GeoDataFrame(con_point, geometry=con_point['geometria'])
con_line = gpd.GeoDataFrame(con_line, geometry=con_line['geometria'])
con_area = gpd.GeoDataFrame(con_area, geometry=con_area['geometria'])
con_volume = gpd.GeoDataFrame(con_volume, geometry=con_volume['geometria'])

# Defina o CRS (Sistema de Referência de Coordenadas)
con_point.set_crs(epsg=4326, inplace=True)
con_line.set_crs(epsg=4326, inplace=True)
con_area.set_crs(epsg=4326, inplace=True)
con_volume.set_crs(epsg=4326, inplace=True)

# Salvar como shapefile
con_point.to_file("C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\shapes\con_point.shp", driver="ESRI Shapefile")
con_line.to_file("C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\shapes\con_line.shp", driver="ESRI Shapefile")
con_area.to_file("C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\shapes\con_area.shp", driver="ESRI Shapefile")
con_volume.to_file("C:\BolsaCongonhas\Git\Congonhas\data\Planilhas\shapes\con_volume.shp", driver="ESRI Shapefile")

