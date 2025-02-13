# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 07:39:03 2025

@author: joseh
"""

#%% Bibliotecas utilizadas

import os
import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import temporalStatistics as ts

#%% Abrir a pasta emission_data

lista_arquivos = [arquivo for arquivo in os.listdir("C:\BolsaCongonhas\Git\Congonhas\E04\emission_data") if arquivo.endswith((".nc", ".ncf"))]

#%% Função para somar as 24 horas de emissão

def soma24horas(poluente, arquivo_nc):
    '''
    Essa função faz a soma das 24 horas de emissão em um determinado dia de um determinado poluente 

    Parâmetros:
    - ds (Dataset): Arquivo NetCDF.
    - poluente (str): Nome da variável dentro do NetCDF (padrão: "PM10").

    Retorna:
    - ds_novo (xarray.Dataset): Dataset atualizado com a soma das 24h adicionada.
    '''

    ds = xr.open_dataset("C:\\BolsaCongonhas\\Git\\Congonhas\\E04\\emission_data\\" + arquivo_nc) 

    xv,yv,lon,lat = ts.ioapiCoords(ds)
    xlon, ylat = ts.eqmerc2latlon(ds, xv, yv)
    
    ylat[:,0] # pega os valores de todas as linhas da primeira coluna
    xlon[0,:] # pega os valores de todas as colunas da primeira linha
    
    if poluente in ds.variables:
        
        # Calcula a soma das 24 primeiras horas na altura 1
        soma_24h = ds[poluente].isel(LAY=0, TSTEP=slice(0, 24)).sum(dim="TSTEP", keepdims=True)
        
        # Adiciona os valores de latitude ao ROW e longitude ao COL
        soma_24h = soma_24h.assign_coords({"ROW": ylat[:,0], "COL": xlon[0,:]})
        
        # Pega qual o tipo de emissor
        soma_24h.attrs["Tipo"] = ds.FILEDESC

        return soma_24h
    
#%% Função para identificar maior emissor

def indentifierMajorSource(lista_dataset):
    '''

    Parâmetros:
    - lista_dataset (list): Lista com todos os datasets.
 
    Retorna:
    - indices (core.dataarray.DataArray): Dataarray com as células preenchidas pelo valor de índice.
    '''
    
    # Filtra os DataArrays que não são None
    lista_dataset = [da for da in lista_dataset if da is not None]
    
    # Empilha os DataArrays ao longo de uma nova dimensão
    da = xr.concat(lista_dataset, dim="source", join="override")
    
    # Define o código para cada tipo de fonte
    codigo_para_tipo = {i: da.attrs.get("Tipo", f"Código {i}") for i, da in enumerate(lista_dataset)}
    
    # Obtém o índice do DataArray que tem o valor máximo para cada célula
    indices = da.argmax(dim="source")
    
    # Adiciona atributos para facilitar a identificação
    indices.attrs["description"] = "Índice do DataArray com maior valor em cada célula"
    
    # Adiciona ao encoding a legenda de códigos para os tipos de fonte
    indices.encoding = codigo_para_tipo

    return indices
   

#%% 

lista_ds = [] # Criação de uma lista dos datasets
pol = 'CO' # Nome do poluente

# Criar uma estrutura de repetição para passar por todos os arquivos netcdf
for num in range(0, len(lista_arquivos)):
    # Aplica a função soma24horas para um determinado poluente de um arquivo netcdf
    ds_24 = soma24horas(pol, lista_arquivos[num])
    
    # Adiciona à lista a soma de 24 horas dos poluentes
    lista_ds.append(ds_24)
    
# Aplica a função para ter o array com o índice dos maiores valores
ds_indices = indentifierMajorSource(lista_ds) 

#%% Printar mapa

fig, ax = plt.subplots(figsize=(8, 6))

# Adiciona o shapefile
shp = gpd.read_file('C:\BolsaCongonhas\Git\Congonhas_LCQAr\shp\Shapefile_Congonhas.shp')

# Reprojetar para EPSG:4326, o sistema de coordenadas geográficas (lat/lon)
shp = shp.to_crs(epsg=4326)
    
img = ax.imshow(ds_indices.data.squeeze(), extent=[
    ds_indices['COL'].min().values, ds_indices['COL'].max().values, 
    ds_indices['ROW'].min().values, ds_indices['ROW'].max().values
])

shp.plot(ax=ax, edgecolor='black', facecolor='none', alpha=0.5)

# Criar uma lista de rótulos para a legenda
labels = [f"{codigo}: {tipo}" for codigo, tipo in ds_indices.encoding.items()]

# Adicionando a legenda ao gráfico
ax.legend(labels, title="Fontes", loc="upper right", fontsize=10)
