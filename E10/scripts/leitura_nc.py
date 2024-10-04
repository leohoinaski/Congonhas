# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 10:03:41 2024

@author: José Henrique Hess

Student in Sanitary and Environmental Engineering
Department of Sanitary and Environmental Engineering of UFSC
"""

#%%

from netCDF4 import Dataset
import os
import netCDF4 as nc
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import temporalStatistics as ts
import contextily as cx

#%% COMO ESTÁ A QUALIDADE DO AR NA REGIÃO DE CONGONHAS?

'''
Função para desenvolver figura da variabilidade espacial.
Função para desenvolver figura com locais de não atendimento dos padrões
Função para desenvolver figura com número de violações no domínio
Função para criar tabela com estatística de cada poluente para cada bairro da cidade. Média anual, máxima, mínima, número de violações.
'''

def variabilidadeEspacial(ds, num, cor, shp):
    
    if num < 3:
        i = 0
        j = num
    else:
        i = 1
        j = num - 3
    
    poluente = ds.attrs['VAR-LIST']
    
    mean_data = ds[poluente].mean(dim='TSTEP')

    mean_data_2d = mean_data.isel(LAY=0)

    mean_values = mean_data_2d.values.astype(float)

    # Exibir a imagem
    img = ax[i,j].imshow(mean_values, cmap=cor, extent=[
        ds['LON'].min().values, ds['LON'].max().values, 
        ds['LAT'].min().values, ds['LAT'].max().values
    ])
    
    # Adicionar Shapefile de Congonhas
    shp.plot(ax=ax[i,j], edgecolor='black', facecolor='none', alpha=0.5)
    
    # Adicionar o fundo de mapa com Contextily
    # cx.add_basemap(ax[i,j], crs=shp.crs.to_string(), source=cx.providers.CartoDB.Positron, alpha=0)
    
    # Remover os valores dos eixos x e y
    if i == 0:
        ax[i,j].set_xticks([])
    if j != 0:
        ax[i,j].set_yticks([])
        
    # Adicionar o colorbar acima do gráfico
    cbar = plt.colorbar(img, ax=ax[i,j], orientation='vertical', pad=0.05, fraction=0.045)
    cbar.ax.yaxis.set_ticks_position('right')
    cbar.ax.yaxis.set_label_position('right')
    
    # Adicionar um título
    ax[i,j].set_title(poluente + ' (' + ds[poluente].attrs['units'] +') ' )

def naoAtendimentoPadrao(ds):
    a = ds
    
    return a

def numeroViolacoes(ds):
    a = ds
    
    return a

def tabelaEstatisticaPoluente(ds):
    a = ds
    
    return a


#%%

# Define o caminho da pasta onde tem as emissoes em netCDF
caminho = 'C:\BolsaCongonhas\Git\Congonhas\data\Con_3km'

# Cria uma lista dos nomes de cada arquivo que há na pasta de caminho
lista_poluentes = os.listdir(caminho)

lista_ds = []

# Carregar o arquivo NetCDF
for num in range(0, len(lista_poluentes)):
    ds = xr.open_dataset(caminho + "\\" + lista_poluentes[num])
    
    lista_ds.append(ds)

#%% Figura da Variabilidade Espacial

fig, ax = plt.subplots(2,3,figsize=(8, 6))

lista_cores = ['jet','jet','jet','jet','jet','jet']

# Adicionar o shapefile
shp = gpd.read_file('C:\BolsaCongonhas\Git\Congonhas_LCQAr\shp\Shapefile_Congonhas.shp')

for num in range(0, len(lista_ds)):
    
    color = lista_cores[num]
    
    variabilidadeEspacial(lista_ds[num], num, color, shp)
    
# Adicionar labels global
fig.suptitle('Média das concentrações para cada poluente de 01-jan-2023 a 04-mar-2023', fontsize=12)
fig.supxlabel('Longitude', fontsize=12)
fig.supylabel('Latitude', fontsize=12)
    
# Salvar Plot
fig.savefig('C:\\BolsaCongonhas\\Git\\Congonhas\\E10\\figures\\variabilidade_espacial.png', dpi=300, bbox_inches='tight')

#%%

#%%

