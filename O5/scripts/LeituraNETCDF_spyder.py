# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:03:45 2024

@author: José Henrique Hess

Student in Sanitary and Environmental Engineering
Department of Sanitary and Environmental Engineering of UFSC

"""
#%%

from netCDF4 import Dataset
import os
import netCDF4 as nc
import temporalStatistics as ts
import xarray as xr
import numpy as np

'''

Dimensões do dataset:

TSTEP: 25 - Espaço de tempo de 0 até 24
VAR: 66 - Número de variáveis, pelo que entendi de layers com diferentes concentrações de cada poluente
DATE-TIME: 2 - Não sei
LAY: 1 - Camada de altura
ROW: 148 - Latitude
COL: 148 - Longitude

Dimensões de um poluente

TSTEP, LAY, ROW, COL
25, 1, 148, 148

'''

#%% 

def listarPoluente(lista_arquivos: list, caminho: str, poluente: str) -> list:
    """
    Este código abre uma pasta com datasets.
    Os datasets necessariamente tem que ter a mesma fonte de emissão.
    Deles é concatenado todas as datas para um poluente selecionado.

    Parameters
    ----------
    lista_arquivos : list
        DESCRIPTION.
    caminho : str
        DESCRIPTION.
    poluente : str
        DESCRIPTION.

    Returns
    -------
    lista_tstep : list
        DESCRIPTION.

    """

    lista_tstep = []
    
    for i in range (0, len(lista_arquivos)):
        
        if lista_arquivos[i].endswith(".nc"):
            
            # Define o caminho completo até o arquivo
            caminho_completo = caminho + '\\' + lista_arquivos[i]
            
            # Carregar o arquivo NetCDF
            ds = xr.open_dataset(caminho_completo)
            
            for j in range(0,24):
                
                ds_pol_tempo = ds[poluente][j,:,:,:]
            
                lista_tstep.append(ds_pol_tempo)
    
    #
    xv,yv,lon,lat = ts.ioapiCoords(ds)
    xlon, ylat = ts.eqmerc2latlon(ds, xv, yv)
    
    for i in range(0,len(lista_tstep)):
       
        # Substituir 'col' por 'longitude' e 'row' por 'latitude'
        lista_tstep[i] = lista_tstep[i].assign_coords({"ROW": ylat[:,0], "COL": xlon[0,:]})
    
    return lista_tstep


def pegarPrimeiroDataset(lista_arquivos: list, caminho: str) -> xr.Dataset():
    """
    

    Parameters
    ----------
    lista_arquivos : list
        DESCRIPTION.
    caminho : str
        DESCRIPTION.

    Returns
    -------
    ds : xr.Dataset
        DESCRIPTION.

    """
    
    j = 0
        
    while j != 9999: 
        
        if lista_arquivos[j].endswith(".nc"):
            
            # Define o caminho completo até o arquivo
            caminho_completo = caminho + '\\' + lista_arquivos[j]
            
            # Carregar o arquivo NetCDF
            ds = xr.open_dataset(caminho_completo)
             
            j = 9998
            
        j = j + 1
        
    return ds
        
    
    
def adequarMetadadosParaVerdi(lista_tstep: list, ds: xr.Dataset(), poluente: str) -> xr.Dataset():
    """
    

    Parameters
    ----------
    lista_tstep : list
        DESCRIPTION.
    ds : xr.Dataset()
        DESCRIPTION.
    poluente : str
        DESCRIPTION.

    Returns
    -------
    combined : xr.Dataset()
        DESCRIPTION.

    """

    # Combine todos os DataArrays ao longo de uma nova dimensão
    combined = xr.concat(lista_tstep, dim='TSTEP')

    # Verifique se combined é um Dataset
    if isinstance(combined, xr.DataArray):
        combined = combined.to_dataset(name=poluente)

    # Copie os atributos globais
    combined.attrs = ds.attrs

    # Copie os atributos das variáveis
    for var_name in ds.data_vars:
        if var_name in combined.data_vars:
            combined[var_name].attrs = ds[var_name].attrs

    combined.attrs['VAR-LIST'] = poluente

    combined = combined.transpose('TSTEP', 'LAY', 'ROW', 'COL')
        
    return combined
  

#%%

# Define o caminho da pasta onde tem as emissoes em netCDF
caminho = 'C:\BolsaCongonhas\Emissoes'

# Cria uma lista dos nomes de cada arquivo que há na pasta de caminho
fontes = os.listdir(caminho)
    
for fonte in fontes:
    
    if fonte == 'veicular':
        
        # Define o caminho completo até determinada fonte
        caminho_fonte = caminho + '\\' + fonte
        
        # Define uma lista de todos os arquivos no caminho    
        caminho_completo = os.listdir(caminho_fonte)
        
        # Carregar o arquivo NetCDF
        ds = xr.open_dataset(caminho_fonte + "\\" + caminho_completo[0])
        
        # Listar as variáveis e retirar o TFLAG
        poluentes = list(ds.variables)
        poluentes.remove('TFLAG')
        
        for poluente in poluentes:
            
            print(poluente)
            
            # Chama a função que cria a lista de determinado poluente e determinada fonte nas datas disponíveis
            lista_tstep = listarPoluente(caminho_completo, caminho_fonte, poluente)
                 
            # Chama a função que pega o primeiro dataset da lista de arquivos
            ds = pegarPrimeiroDataset(caminho_completo, caminho_fonte)
            
            # Chama a função que adequa os metadados para o Verdi
            combined = adequarMetadadosParaVerdi(lista_tstep, ds, poluente)
            
            # Salva os metadadados em determinada pasta
            combined.to_netcdf("C:\BolsaCongonhas\Git\Congonhas_LCQAr\emissoes" + "\\" + fonte + "\\" + fonte + "_" + poluente + ".nc" )

#%%

def somarTodosPoluentes(ds):

    ds_sum = ds.sum(dim='TSTEP', keep_attrs=True)
    first_tstep = ds['TSTEP'].isel(TSTEP=0)
    ds_sum = ds_sum.expand_dims('TSTEP')
    ds_sum['TSTEP'] = [first_tstep.values]
    
    return ds_sum

def maiorFontePoluidora(ds):
    
    a = ds
    
    return a

def somarNOx(dsNO, dsNO2):
    
    ds1_aligned, ds2_aligned = xr.align(dsNO, dsNO2, join="inner")

    # Agora, faça a soma
    dsNOx = ds1_aligned + ds2_aligned
    
    return dsNOx

'''
for pol in poluentes:
    print(pol)
    # Selecting variable
    if fileType == 'BRAVESdatabase2CMAQ' and pol==PM10:
        data = ds['PMC'][:]
        polu = {'tag': 'PM10',
            'Unit': ds['PMC'].units}

    #REVISAR ISTO
        #data=ATOTI*1+ATOTJ*1+ATOTK*0.5

    elif fileType == 'GLOB_GEOSchem_MG_3km' and pol==PM10:

        sources='FINN'

        intermediateFilePrefix = 'FINN2D'
        pspec =['POC','PEC','PSO4','PNO3','PMOTHR']
        for kk,pref in enumerate(prefixed):
            dsi = nc.MFDataset(pref)
            data = np.zeros((dsi[pspec[0]][:].shape[0], 1,dsi[pspec[0]][:].shape[2], dsi[pspec[0]][:].shape[3]))
            print(data.shape)
            for ps in pspec:
                data[:,0,:,:] = data[:,0,:,:]+ np.nansum(dsi[ps][:],axis=1)
            netCDFEmiswriter(dsi,data,sources,polu,'FINN2D_',folderOut)
        polu = {'tag': pol['tag'],
            'Unit': dsi[ps].units}
        prefixed2 = prefixed = sorted([filename for filename in os.listdir(path[count]) if filename.startswith(fileType)])
        #prefixed2 =  sorted([filename for filename in os.listdir('/home/artaxo/CMAQ_REPO/PREP/emis/finn2cmaq-master/hourly/2021/09') if filename.startswith('FINN2D_'+sources+>
        print(prefixed2)
        ds = nc.MFDataset(prefixed2)
        data = ds[ps][:]
        polu = {'tag': pol['tag'],
            'Unit': ds[ps].units}
        intermediateFileRemover(folderOut,intermediateFilePrefix)

    else:
        data = ds[:]
'''

#%%

poluentes_veicular = ['CO', 'PMFINE', 'NO2','NO','NOx', 'SO2', 'PM10']

fonte = 'veicular'

for poluente in poluentes_veicular:
        
    if poluente == 'NOx':
        
        dsNO = xr.open_dataset("C:\BolsaCongonhas\Git\Congonhas_LCQAr\emissoes" + "\\" + fonte + "\\" + fonte + "_NO.nc")
        dsNO2 = xr.open_dataset("C:\BolsaCongonhas\Git\Congonhas_LCQAr\emissoes" + "\\" + fonte + "\\" + fonte + "_NO2.nc")

        ds = somarNOx(dsNO, dsNO2)
    
    else:
        # Carregar o arquivo NetCDF
        ds = xr.open_dataset("C:\BolsaCongonhas\Git\Congonhas_LCQAr\emissoes" + "\\" + fonte + "\\" + fonte + "_" + poluente + ".nc")

    ds = somarTodosPoluentes(ds)
    
    ds.to_netcdf("C:\BolsaCongonhas" + "\\" + fonte + "_" + poluente + ".nc" )

#%%

'''
for i in range (0, len(arquivos)):
    
    if arquivos[i].endswith(".nc"):
        
        ds = listarPoluente(arquivos, caminho, 'PMOTHR')
            
        
            
#%% Uso das funçoes de temporalStatistics para corrigir as coordenadas do dataset

xv,yv,lon,lat = ts.ioapiCoords(ds)
xlon, ylat = ts.eqmerc2latlon(ds, xv, yv)

#%%

# ylat[:,0] pega os valores de todas as linhas da primeira coluna
# xlon[0,:] pega os valores de todas as colunas da primeira linha

for i in range(0,len(lista_tstep)):
   
    # Substituir 'col' por 'longitude' e 'row' por 'latitude'
    lista_tstep[i] = lista_tstep[i].assign_coords({"ROW": ylat[:,0], "COL": xlon[0,:]})
    
#%%

j = 0
    
while j != 9999: 
    
    if arquivos[j].endswith(".nc"):
        
        ds = listarPoluente(j, arquivos, caminho, 'PMOTHR')
         
        j = 9998
        
    j = j + 1

# Combine todos os DataArrays ao longo de uma nova dimensão
combined = xr.concat(lista_tstep, dim='TSTEP')

# Verifique se combined é um Dataset
if isinstance(combined, xr.DataArray):
    combined = combined.to_dataset(name='PMOTHR')

# Copie os atributos globais
combined.attrs = ds.attrs

# Copie os atributos das variáveis
for var_name in ds.data_vars:
    if var_name in combined.data_vars:
        combined[var_name].attrs = ds[var_name].attrs

combined.attrs['VAR-LIST'] = 'PMOTHR'

#%%

combined = combined.transpose('TSTEP', 'LAY', 'ROW', 'COL')
'''
# Exporte o DataArray combinado para um arquivo NetCDF
ds.to_netcdf("combined_data.nc")

#%%
    
# Define o caminho completo até o primeiro arquivo
caminho_completo = 'C:\\BolsaCongonhas\\Git\\Congonhas_LCQAr\\combined_data.nc'

# Carregar o arquivo NetCDF
ds_todos = xr.open_dataset(caminho_completo)

#%%

wbd_tstep = []

for i in range (0, len(ds_todos.TSTEP)):
    
    wbd = ds_todos.PMOTHR[i,:,:,:]
    wbd_tstep.append(wbd)

#%%

import matplotlib.pyplot as plt
import geopandas as gpd

caminho_shapefile = 'C:\BolsaCongonhas\Git\Congonhas_LCQAr\shp\Shapefile_Congonhas.shp'
shapefile = gpd.read_file(caminho_shapefile)

# Supondo que seu DataArray esteja armazenado na variável `data_array`
# Crie uma figura e um eixo
plt.figure(figsize=(8, 6))

# Use extent para definir os valores de longitude e latitude nos eixos
plt.imshow(ds_todos[16], origin='lower', cmap='YlOrRd',
           extent=[xlon.min(), xlon.max(), ylat.min(), ylat.max()])

# Adicione o shapefile ao gráfico
shapefile.plot(ax=plt.gca(), edgecolor='blue', facecolor='none', linewidth=1)

# Adicione uma barra de cores para referência
plt.colorbar(label='Valor da variável')

# Adicione rótulos aos eixos
plt.xlabel('Longitude (COL)')
plt.ylabel('Latitude (ROW)')

# Adicione um título ao gráfico (opcional)
plt.title('Imagem WindBlowDust')

# Exibe o gráfico
plt.show()

#%%

import matplotlib.pyplot as plt
import geopandas as gpd

caminho_shapefile = 'C:\BolsaCongonhas\Git\Congonhas_LCQAr\shp\Shapefile_Congonhas.shp'
shapefile = gpd.read_file(caminho_shapefile)

f = 0

for i in range(0,len(wb_tstep)):
    # Supondo que seu DataArray esteja armazenado na variável `data_array`
    # Crie uma figura e um eixo
    plt.figure(figsize=(8, 6))
    
    # Use extent para definir os valores de longitude e latitude nos eixos
    plt.imshow(wb_tstep[i], origin='lower', cmap='YlOrRd',
               extent=[xlon.min(), xlon.max(), ylat.min(), ylat.max()])
    
    # Adicione o shapefile ao gráfico
    shapefile.plot(ax=plt.gca(), edgecolor='blue', facecolor='none', linewidth=1)
    
    # Adicione uma barra de cores para referência
    plt.colorbar(label='Valor da variável')
    
    # Adicione rótulos aos eixos
    plt.xlabel('Longitude (COL)')
    plt.ylabel('Latitude (ROW)')
    
    # Adicione um título ao gráfico (opcional)
    plt.title('Imagem WindBlowDust hora ' + str(i) + ' dia ' + str(f))
    
    # Exibe o gráfico
    plt.show()
    
    print(i)
    
    if i % 24:
        f = f + 1
        print(f)


