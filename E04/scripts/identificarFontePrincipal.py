# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 07:39:03 2025

@author: 
    Camilo Bastos Ribeiro
    José Henrique Hess
"""

#%% Bibliotecas utilizadas

import os
import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import temporalStatistics as tst
import pandas as pd
import numpy as np
from shapely.geometry import Point
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import shutil
import netCDF4 as nc
import BRAINutils as bu

#%% Ignorando avisos de atualização

import warnings
warnings.filterwarnings('ignore')

#%% Função para agregar emissões

def aggEmis(dir_folder, var, op, freq):
    
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
    
    # Cria final_df, sendo um dicionário caso freq for anual e dataframe caso não
    final_df = {} if freq != 'yearly' else pd.DataFrame()
    
    # Definir os caminhos específicos das pastas desejadas
    pasta_veiculo = os.path.join(dir_folder, "BRAVES")
    pasta_queimada = os.path.join(dir_folder, "FINN", "hourly")
    pasta_industrial = os.path.join(dir_folder, "IND2CMAQ")
    pasta_biogenica = os.path.join(dir_folder, "MEGAN")
    pasta_smoke = os.path.join(dir_folder, "smoke")
    pasta_unpaved = os.path.join(dir_folder, "unpaved_emission")
    pasta_wbd = os.path.join(dir_folder, "windBlowDustBR")
    
    pastas = [pasta_veiculo, pasta_queimada, pasta_industrial, pasta_biogenica, 
              pasta_smoke, pasta_unpaved,pasta_wbd]
    
    # Cria uma iteração das pastas dentro da lista pastas
    for pasta in pastas:
    
        # Cria uma iteração dos arquivos dentro da pasta dir_folder
        for file in os.listdir(pasta):
            
            # Condição if para caso arquivo seja netcdf
            if file.endswith(('.nc','.ncf')) and file.startswith(
                    ('BRAVESdatabase2CMAQ','IND2CMAQ','MEGANv31',
                     'agts','modified','windBlowDust_PM10', 'GLOB_GEOS')):
                
                # Para facilitar o nome, do arquivo, pega a primeira parte antes do _
                sector_name = file.split('_')[0] 
                
                if sector_name == 'agts':
                    sector_name = file.split('.')[1] 
                    print(sector_name)
                
                # Cria um dicionário para mapear o nome dos arquivos
                mapeamento = {
                    'AgrWstBrn': 'AgrWstBrn',
                    'DomAvi': 'DomAvi',
                    'DomShip': 'DomShip',
                    'IntAvi': 'IntAvi',
                    'IntShip': 'IntShip',
                    'Lstock': 'Lstock',
                    'Resi': 'Resi',
                    'Solvents': 'Solvents',
                    'Waste': 'Waste',
                    'GLOB': 'Queimadas',
                    'BRAVESdatabase2CMAQ': 'Vehicular',
                    'IND2CMAQ': 'Industrial',
                    'MEGANv31.Con': 'Biogenichals',
                    'modified': 'Unpaved',
                    'windBlowDust': 'Difuse'
                }
                
                # Altera a variável sector_name utilizando o dicionário mapeamento
                sector_name = mapeamento.get(sector_name, sector_name)
                  
                print(f"{sector_name}")
                print(file.split('_')[-1])
                print('\n')
                    
                dir_data = os.path.join(pasta, file) # Cria o diretório correto para o arquivo
                data = xr.open_dataset(dir_data) # Abre o netcdf
                
                # Cria uma lista dos poluentes que não tem PM10 como variável
                lista_PM10 = ['AgrWstBrn','DomAvi','DomShip','IntAvi','IntShip',
                              'Lstock','Resi','Solvents','Waste',
                              'Vehicular','Industrial','Unpaved','Queimadas']
                
                # Condição if para verificar se a variável analisada é PM10 e se o arquivo está na lista_PM10
                if var == 'PM10' and sector_name in lista_PM10:
                    if sector_name == 'Industrial' or sector_name == 'Unpaved' or sector_name == 'Queimadas':
                        data = addPM10(sector_name,data,pasta) # Utiliza a função addPM10 para adicionar PM10 nas fontes industriais e rodovias
                    elif sector_name == 'Vehicular':
                        data["PM10"] = data['PMC'] # Cria o PM10 nos veículos a partir do PMC
                    else:
                        data["PM10"] = data['PM10_INV'] # Cria o PM10 a partir do PM10_INV, INV = INVENTORY
                
                # Condição if para verificar se var está presente nas variáveis do netcdf
                if var in data.variables:
                    
                    # Pega os valores de tflag do netcdf
                    tflag = data['TFLAG'].values
                    
                    # Cria uma fatia do tflag com os primeiros valores de cada linha e todas as horas
                    tflag = tflag[:, 0, :]
                    
                    ano = tflag[:, 0] // 1000   # Pega os 4 primeiros dígitos como ano
                    dia_juliano = tflag[:, 0] % 1000  # Últimos 3 dígitos são o dia do ano
                    horas = tflag[:, 1] // 10000 # Pega as duas primeiras posições
                    
                    # Cria um dataframe com as colunas year, day_of_year e hour, para adicionar as variáveis geradas
                    df = pd.DataFrame({
                        'year': ano,
                        'day_of_year': dia_juliano,
                        'hour': horas
                    })
    
                    # Cria uma string no formato 'YYYYDDD' e converter para datetime
                    time = pd.to_datetime(df['year'].astype(str) + df['day_of_year'].astype(str).str.zfill(3), format='%Y%j') 
                    
                    # Adiciona as respectivas horas tirado 3 horas devido ao GMT em Congonhas
                    time = time + pd.to_timedelta(df['hour'], unit='h') - pd.to_timedelta(3, unit='h')
                    
                    if freq == 'monthly':
                        time_group = time.dt.month 
                        unique_times = np.unique(time_group)
                    
                    elif freq == 'weekly':
                        time_group = time.dt.weekday
                        unique_times = np.arange(7)
                    
                    elif freq == 'hourly':
                        time_group = time.dt.hour
                        unique_times = np.arange(24)
                    
                    elif freq == 'yearly':
                        time_group = None
                    
                    #Op for defined freq
                    if freq != 'yearly':
                        for t in unique_times:
                            time_indices = np.where(time_group == t)[0]
                            
                            if sector_name == 'Difuse':
                                #Op in the time inverval
                                pol_2d = pd.DataFrame(
                                    op(np.array(data[var][time_indices, :, :]), axis=0).flatten()
                                ).rename(columns={0: sector_name})
                            else:
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
                        if sector_name == 'Difuse':
                            pol_2d = pd.DataFrame(
                                op(np.array(data[var][time_indices, :, :]), axis=0).flatten()
                            ).rename(columns={0: sector_name})
                        else:
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
    
    if type(final_df) == dict:
        for t in unique_times:
            
            
            # Renomear colunas repetidas para evitar sobrescrita
            final_df[t].columns = pd.Index([f"{col}_{i}" for i, col in enumerate(final_df[t].columns)])
            
            # Criar um dicionário para agrupar colunas pelo prefixo (antes do "_")
            groups = {}
            for col in final_df[t].columns:
                base_name = col.split("_")[0]  # Pega o nome original da coluna
                groups.setdefault(base_name, []).append(col)
            
            # Calcular a média para cada grupo de colunas
            final_df[t] = pd.DataFrame({col: op(final_df[t][cols], axis=1) for col, cols in groups.items()})
        
    else:
        # Renomear colunas repetidas para evitar sobrescrita
        final_df.columns = pd.Index([f"{col}_{i}" for i, col in enumerate(final_df.columns)])
        
        # Criar um dicionário para agrupar colunas pelo prefixo (antes do "_")
        groups = {}
        for col in final_df.columns:
            base_name = col.split("_")[0]  # Pega o nome original da coluna
            groups.setdefault(base_name, []).append(col)
        
        # Calcular a média para cada grupo de colunas
        final_df = pd.DataFrame({col: op(final_df[cols], axis=1) for col, cols in groups.items()})

    # Get the major emitter
    if freq != 'yearly':
        for t, df in final_df.items():
            df['major_emitter'] = df.idxmax(axis=1)
            df.loc[df.iloc[:, :-1].max(axis=1) == 0, 'major_emitter'] = np.nan
    else:
        # Get the major emitter for 'yearly' freq
        final_df['major_emitter'] = final_df.idxmax(axis=1)
        final_df.loc[final_df.iloc[:, :-1].max(axis=1) == 0, 'major_emitter'] = np.nan
    
    # Add the total emissions column after identifying the major emitter
    if freq != 'yearly':
        for t, df in final_df.items():
            df['total_emissions'] = df.drop(columns=['major_emitter']).sum(axis=1)
    else:
        final_df['total_emissions'] = final_df.drop(columns=['major_emitter']).sum(axis=1)
    
    return final_df
    
#%% Coletar lat lon

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
    
#%% Função para identificar maior emissor

def highEmitter(dfs, lat, lon, shp, freq, var):
    
    """
    

    Parameters
    ----------
    dfs : Dictionary ou DataFrame 
            que possui os valores de emissão para cada poluente em cada pixel, 
            qual o maior emissor e qual a soma de emissão
    lat : Array of float64 (e.g. (21904,)) 
            com os valores de latitude para todos os pixels
    lon : Array of float64 (e.g. (21904,)) 
            com os valores de longitude para todos os pixels
        DESCRIPTION.
    shp : DataFrame
            com coluna geometry que possui um POLYGON
    freq : string
            que indica a frequência temporal (e.g. hourly, weekly, monthly, yearly)
    var : string
            que indica o poluente que está sendo avaliado

    Returns
    -------
    None.

    """
    
    # Condição if para verificar se dfs é um dicionário
    if type(dfs) == dict:
        lista_dfs = list(dfs.values()) # Transforma os valores do dicionário em uma lista
    else:
        lista_dfs = [] # Cria uma lista chama lista_dfs
        lista_dfs.append(dfs) # Adiciona o DataFrame dfs à lista lista_dfs 
    
    lista_gdfs = [] # Cria uma lista chamada lista_gdfs que receberá geoDataFrames
    
    color_map = {
        'AgrWstBrn': '#FF5733',
        'DomAvi': '#33FF57',
        'DomShip': '#3357FF',
        'IntAvi': '#FFFF33',
        'IntShip': '#FF33FF',
        'Lstock': '#33FFFF',
        'Resi': '#FF8C33',
        'Solvents': '#8C33FF',
        'Waste': '#33FF8C',
        'Queimadas': '#FF3333',
        'Vehicular': '#3333FF',
        'Industrial': '#FF33A1',
        'Biogênicas': '#A133FF',
        'Unpaved': '#33A1FF',
        'Difuse': '#A1FF33'
    }
    
    # Iteração for para i como o índice de lista_dfs e df como o DataFrame selecionado
    for i, df in enumerate(lista_dfs):
        df["longitude"] = lon # Adiciona coluna de Longitude
        df["latitude"] = lat # Adiciona coluna de Latitude
        df["geometry"] = df.apply(lambda row: Point(row["longitude"], row["latitude"]), axis=1) # Cria coluna de geometria a partir das duas anteriores
        
        # Converte para geoDataFrame utilizando a coluna geometry e crs = EPSG:4326
        gdf = gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")
        lista_gdfs.append(gdf) # Adiciona o geoDataFrame à lista lista_gdfs
        
        # Adiciona uma nova coluna de cor ao GeoDataFrame
        gdf['color'] = gdf['major_emitter'].map(color_map)
        
        # Substituir valores NaN por uma cor padrão (exemplo: cinza)
        gdf['color'] = gdf['color'].fillna('#A9A9A9')  # Cinza escuro
        
        # Cria o Envelope do Buffer
        gdf['envelope'] = gdf.geometry.buffer(0.0135).envelope  
        
        # Cria um GeoDataFrame apenas com o envelope
        gdf = gpd.GeoDataFrame(gdf, geometry='envelope', crs="EPSG:4326")
             
        for i in range(2):
        
            # Define o tamanho da figura
            fig, ax = plt.subplots(figsize=(12, 12))
            
            # Plotar o GeoDataFrame com as cores definidas
            gdf.plot(color=gdf['color'], ax=ax)
            
            # Criar a legenda manualmente
            legend_patches = [mpatches.Patch(color=color, label=label) for label, color in color_map.items()]
            
            # Adicionar a legenda à direita do plot
            ax.legend(handles=legend_patches, title="Major Emitter", loc="center left", bbox_to_anchor=(1, 0.5))
    
            # Plota o shapefile 'shp' com borda preta e face preta
            shp.plot(ax=ax, edgecolor='black', facecolor='none', alpha=0.5)
            
            if i == 0:
                # Salva a figura utilizando para nomear as variáveis var, freq e i
                fig.savefig(r'C:\BolsaCongonhas\Git\Congonhas\E04\figures\high_emitter_' + var + '_' + freq + '_' + str(i) + '.png')
            else:
                # Definir os limites do gráfico
                ax.set_xlim(-44.35, -43.35)
                ax.set_ylim(-21, -20)
                
                # Salva a figura utilizando para nomear as variáveis var, freq e i
                fig.savefig(r'C:\BolsaCongonhas\Git\Congonhas\E04\figures\high_emitter_Congonhas_' + var + '_' + freq + '_' + str(i) + '.png')
            
#%% Função para identificar maior emissão em tempo

def highTime(dfs, lat, lon, shp, freq, var, ind_val):
    '''
    Identifica o período com a maior emissão para cada pixel (linha) e cada setor (coluna)
    e seus respectivos valores.
    
    Parameters:
    - dfs: list of DataFrames, each corresponding to a different time period (e.g., months, days, hours),
           where each column represents an emission sector.
    - time_labels: list of strings representing the labels for each time period, in the same order as the DataFrames.
    
    Returns:
    - A DataFrame where each cell contains the label of the time period with the highest emission for the respective pixel and sector.

    Parameters
    ----------
    dfs : Dictionary ou DataFrame 
            que possui os valores de emissão para cada poluente em cada pixel, 
            qual o maior emissor e qual a soma de emissão
    lat : Array of float64 (e.g. (21904,)) 
            com os valores de latitude para todos os pixels
    lon : Array of float64 (e.g. (21904,)) 
            com os valores de longitude para todos os pixels
        DESCRIPTION.
    shp : DataFrame
            com coluna geometry que possui um POLYGON
    freq : string
            que indica a frequência temporal (e.g. hourly, weekly, monthly, yearly)
    var : string
            que indica o poluente que está sendo avaliado
    ind_val : string
            indica se fará uma imagem dos índices (e.g. hora com o maior valor) 
            ou então valores (e.g. maior valor dentre as horas)

    Returns
    -------
    None

    '''
    
    # Condição if para verificar se dfs não é um dicionário, caso seja, rodará o código
    if type(dfs) != dict:
        return 'O primeiro argumento não está no formato dict ou há apenas um dataframe'
    
    # Condição if para verificar a frequência e assim gerar o label para a plotagem das figuras
    if freq == 'monthly':
        labels = ['Jan', 'Fev', 'Mar', 'Abr', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']    
    elif freq == 'weekly':
        labels = ['Seg', 'Ter', 'Qua', 'Qui', 'Sex', 'Sab', 'Dom']
    elif freq == 'hourly':
        labels = list(np.arange(24))
    
    # Criar uma lista para armazenar a chave com maior valor para cada linha
    max_keys = []
    
    # Transforma os DataFrames em uma única matriz com as chaves
    df_combined = pd.DataFrame({key: df['total_emissions'] for key, df in dict_dfs.items()})
    
    # Pega a chave correspondente ao maior valor em cada linha
    max_keys = df_combined.idxmax(axis=1)  # Obtém os índices dos maiores valores
    max_keys[df_combined.nunique(axis=1) == 1] = np.nan  # Substitui por NaN se todos os valores forem iguais
    max_keys = max_keys.tolist()  # Converte para lista os índices com maior valor
    max_values = df_combined.max(axis=1).tolist() # Converte para lista os maiores valores

    # Cria um dataframe
    df = pd.DataFrame({
        'longitude': lon, # Adiciona uma coluna com longitude
        'latitude': lat, # Adiciona uma coluna com latitude
        'indice': max_keys, # Adiciona uma coluna com os índices que possuem maior valor
        'valor': max_values}) # Adiciona uma coluna com os maiores valores
    
    # Cria a coluna geometria a partir de longitude e latitude no DataFrame
    df["geometry"] = df.apply(lambda row: Point(row["longitude"], row["latitude"]), axis=1)
    
    # Converte o df para GeoDataFrame
    gdf = gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")
    
    # Cria o Envelope do Buffer
    gdf['envelope'] = gdf.geometry.buffer(0.0135).envelope  
    
    # Cria um GeoDataFrame apenas com o envelope
    gdf = gpd.GeoDataFrame(gdf, geometry='envelope', crs="EPSG:4326")
            
    # Define o tamanho da figura
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Cria um colormap discreto baseado na quantidade de categorias
    cmap = plt.get_cmap('jet', len(labels))
    
    # Condição if para verificar se ind_val é indice ou valor
    if ind_val == 'indice':
        
        # Cria uma escala de 0 até tamanho de labels
        norm = mcolors.Normalize(vmin=0, vmax=len(labels))  
        
        # Cria um mapeador para a legenda
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # Necessário para que funcione com colorbar
        
        # Adiciona a colorbar com rótulos personalizados
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label(freq)  # Nome da legenda
        
        # Adiciona os rótulos
        cbar.set_ticks(range(len(labels)))  
        cbar.set_ticklabels(labels)
       
    elif ind_val == 'valor': 
        # Cria uma escala logarítmica para valor
        norm = mcolors.LogNorm(vmin=gdf['valor'].min()+10**-5, vmax=gdf['valor'].max())
        
    for i in range(2):
        
        # Plota gdf    
        gdf.plot(column=ind_val, cmap=cmap, ax=ax, legend=False, norm = norm)
    
        # Plota o shapefile 'shp' com borda preta e face preta
        shp.plot(ax=ax, edgecolor='black', facecolor='none', alpha=0.5)
        
        if i == 0:
            # Salva a figura
            fig.savefig(r'C:\BolsaCongonhas\Git\Congonhas\E04\figures\high_' + ind_val + '_' + var + '_' + freq + '.png')
        else:
            # Definir os limites do gráfico
            ax.set_xlim(-44.35, -43.35)
            ax.set_ylim(-21, -20)

            # Salva a figura
            fig.savefig(r'C:\BolsaCongonhas\Git\Congonhas\E04\figures\high_Congonhas_' + ind_val + '_' + var + '_' + freq + '.png')

    # Define a pasta de destino para salvar o gdf
    output_folder = r'C:\BolsaCongonhas\Git\Congonhas\E04\geodfs'  # Substitua pelo seu caminho real
    
    # Caminho do arquivo CSV
    output_path = os.path.join(output_folder, f"high_{ind_val}_{var}_{freq}.csv")
    
    # Converte a geometria para WKT e salva em csv
    gdf["geometry"] = gdf["geometry"].apply(lambda geom: geom.wkt)  # Converte para texto
    gdf.drop(columns=["geometry","envelope"], inplace=True)  # Remove as colunas geometry e envelope
    gdf.to_csv(output_path, index=False) 
    
    


'''
'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 
'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Grays', 'Greens', 'Greens_r', 'Greys', 
'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 
'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 
'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 
'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 
'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 
'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 
'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 
'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 
'gist_gray_r', 'gist_grey', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 
'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gist_yerg', 'gnuplot', 'gnuplot2', 
'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'grey', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 
'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 
'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 
'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 
'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 
'viridis_r', 'winter', 'winter_r'
'''

#%% Funções para o PM10

def netCDFEmiswriter(ds,data,sourceID,polu,name,folderOut):
    
    if ~np.isnan(data).all():
        # Get datesTime and removing duplicates
        datesTime, data = tst.getTime(ds,data)
        datesTimeAll = datesTime.copy()
        #datesTimeAll = pd.to_datetime(datesTimeAll, format='%Y%m%d%H')
        # Get coordinates from ioapi
        xv,yv,lon,lat = tst.ioapiCoords(ds)
        '''print(ds)
        print(data)
        print(sourceID)
        print(polu)
        print(name)
        print(folderOut)'''
        # Transforming mercator to latlon/degrees
        xlon, ylat = tst.eqmerc2latlon(ds,xv,yv)
        bu.createNETCDFtemporal(folderOut,name+sourceID+
                             '_'+polu['tag']+'_'+
                             str(datesTimeAll.year[0])+'_'+
                             str(datesTimeAll.month[0]).zfill(2)+'_'+
                             str(datesTimeAll.day[0]).zfill(2)+'_'+
                             str(datesTimeAll.hour[0]).zfill(2)+'_to_'+
                             str(datesTimeAll.iloc[-1].year)+'_'+
                             str(datesTimeAll.iloc[-1].month).zfill(2)+'_'+
                             str(datesTimeAll.iloc[-1].day).zfill(2)+'_'+
                             str(datesTimeAll.iloc[-1].hour).zfill(2)+
                             '.nc',data,ds,polu,xlon,ylat,datesTime)
        
    return data

def intermediateFileRemover(folderOut,intermediateFilePrefix):
    for fname in os.listdir(folderOut):
        if fname.startswith(intermediateFilePrefix):
            os.remove(os.path.join(folderOut, fname))

    return intermediateFilePrefix

def addPM10(name,dataset,path):
    '''
    

    Returns
    -------
    None.

    '''
    
    PM10 = {
       "Pollutant": "$PM_{10}$",
       "Unit": '$\u03BCg.m^{-3}$',
       "tag":'PM10'
     }

    pol = [PM10]

    folderOut=r'C:\BolsaCongonhas\Git\Congonhas\E04\out'

    tflag = dataset['TFLAG'].values
    
    # Extrai o ano e o dia juliano
    dia = tflag[:, 0, :][0, 0] % 1000  # Últimos 3 dígitos são o dia do ano (DDD)
    
    dias_por_mes = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    mes = 1  # Começa em janeiro
    while dia > dias_por_mes[mes - 1]:  # Ajuste no índice (0 a 11)
        dia -= dias_por_mes[mes - 1]
        mes += 1

    if name == 'Unpaved':
        fileType = 'modified_unpaved_Speciated'
        dataFile = '_2023-' + str(f"{mes:02d}") + '-' + str(f"{dia:02d}")
    elif name == 'Industrial':
        fileType = 'IND2CMAQ_2023'
        dataFile = '_' + str(f"{mes:02d}") + '_' + str(f"{dia:02d}")
    elif name == 'Queimadas':
        fileType = 'GLOB_GEOSchem_Con_3km.3D.2023'
        dataFile = '-' + str(f"{mes:02d}") + '-' + str(f"{dia:02d}")

    print(fileType + dataFile)

    # Selecting files and variables
    prefixed = sorted([path + '/' + filename for filename in os.listdir(path) if filename.startswith(fileType + dataFile)])
    

    # Opening netCDF files
    ds = nc.MFDataset(prefixed)
   
    if fileType=='IND2CMAQ_2023':
        intermediateFilePrefix='IND_'
        for pref in prefixed:
            dsi = nc.MFDataset(pref)
            data = np.zeros((dsi['PMC'][:].shape[0], 1,dsi['PMC'][:].shape[2], dsi['PMC'][:].shape[3]))
            data[:,0,:,:] = np.nansum(dsi['CO'][:],axis=1)
            '''polu = {'tag': pol['tag'],
                'Unit': dsi['PMC'].units}'''
            ''''netCDFEmiswriter(dsi,data,fileType,polu,'IND2D_',folderOut)'''
        
        dataset["PM10"] = dataset["PMC"].copy()  # Criando a variável nova
        dataset["PM10"].loc[:, 0, :, :] = data[:, 0, :, :]  # Atribuindo apenas para LAY=0
        
        '''prefixed2 = prefixed = sorted([filename for filename in os.listdir(path) if filename.startswith(fileType)])
        ds = nc.MFDataset(prefixed2)
        data = ds['PMC'][:]
        polu = {'tag': pol['tag'],
            'Unit': dsi['PMC'].units}
        #intermediateFileRemover(folderOut,intermediateFilePrefix)
        ds["PM10"] = (ds['PMC'].dims, ds['PMC'][:])  '''

    #REVISAR ISTO
        #data=ATOTI*1+ATOTJ*1+ATOTK*0.5


    elif fileType == 'GLOB_GEOSchem_Con_3km.3D.2023':

        sources='FINN'
        intermediateFilePrefix = 'FINN2D'
        pspec =['POC','PEC','PSO4','PNO3','PMOTHR']
        for kk,pref in enumerate(prefixed):
            dsi = nc.MFDataset(pref)
            data = np.zeros((dsi[pspec[0]][:].shape[0], 1,dsi[pspec[0]][:].shape[2], dsi[pspec[0]][:].shape[3]))
            for ps in pspec:
                data[:,0,:,:] = data[:,0,:,:]+ np.nansum(dsi[ps][:],axis=1)                    
                '''polu = {'tag': pol['tag'],
                    'Unit': dsi[ps].units}'''
            #netCDFEmiswriter(dsi,data,sources,polu,'FINN2D_',folderOut)
        '''polu = {'tag': pol['tag'],
            'Unit': dsi[ps].units}
        prefixed2 = prefixed = sorted([filename for filename in os.listdir(path) if filename.startswith(fileType)])'''
        
        dataset["PM10"] = dataset[ps].copy()  # Criando a variável nova
        dataset["PM10"].loc[:, 0, :, :] = data[:, 0, :, :]  # Atribuindo apenas para LAY=0
        
        '''
        #prefixed2 =  sorted([filename for filename in os.listdir('/home/artaxo/CMAQ_REPO/PREP/emis/finn2cmaq-master/hourly/2021/09') if filename.startswith('FINN2D_'+sources+'_'+polu['tag'])])
        print(prefixed2)
        ds = nc.MFDataset(prefixed2)
        data = ds[ps][:]
        polu = {'tag': pol['tag'],
            'Unit': ds[ps].units}
        #intermediateFileRemover(folderOut,intermediateFilePrefix)'''

    elif fileType == 'modified_unpaved_Speciated':

        sources='unpaved'

        intermediateFilePrefix = 'unpaved'
        pspec =['PAL','PCA','PFE','PK','PMN','PSI','PTI']
        for kk,pref in enumerate(prefixed):
            dsi = nc.MFDataset(pref)
            data = np.zeros((dsi[pspec[0]][:].shape[0], 1,dsi[pspec[0]][:].shape[2], dsi[pspec[0]][:].shape[3]))
            for ps in pspec:
                data[:,0,:,:] = data[:,0,:,:]+ np.nansum(dsi[ps][:],axis=1)                    
                '''polu = {'tag': pol['tag'],
                    'Unit': dsi[ps].units}'''
            #netCDFEmiswriter(dsi,data,sources,polu,'modified_',folderOut)
            
        dataset["PM10"] = (dataset[ps].dims, data) 
        '''
        polu = {'tag': pol['tag'],
            'Unit': dsi[ps].units}
        prefixed2 = prefixed = sorted([filename for filename in os.listdir(path) if filename.startswith(fileType)])
        #prefixed2 =  sorted([filename for filename in os.listdir('/home/artaxo/CMAQ_REPO/PREP/emis/finn2cmaq-master/hourly/2021/09') if filename.startswith('FINN2D_'+sources+'_'+polu['tag'])])
        print(prefixed2)
        ds = nc.MFDataset(prefixed2)
        data = ds[ps][:]
        polu = {'tag': pol['tag'],
            'Unit': ds[ps].units}
        #intermediateFileRemover(folderOut,intermediateFilePrefix)'''
 

    else:
        data = ds[pol['tag']][:]

    return dataset


#%% Rodar códigos

dir_folder = "C:\BolsaCongonhas\Git\Congonhas\E04\emission_data"
var = 'PM10' # Nome do poluente
op = np.mean
freq = 'weekly' # 'monthly', 'weekly', 'hourly', or 'yearly'

dict_dfs = aggEmis(dir_folder, var, op, freq)

lon, lat = latlon_2d("C:\BolsaCongonhas\Git\Congonhas\E04\emission_data\wbd_updated_PM10_2023-01-01.nc")
shp = gpd.read_file('C:\BolsaCongonhas\Git\Congonhas_LCQAr\shp\Shapefile_Congonhas.shp')

highEmitter(dict_dfs, lat, lon, shp, freq, var)

ind_val = 'indice' # 'indice', 'valor'

highTime(dict_dfs, lat, lon, shp, freq, var, ind_val)

#%% Abrir os arquivos em lista

#lista_arquivos = [arquivo for arquivo in os.listdir("C:\BolsaCongonhas\Git\Congonhas\E04\emission_data") if arquivo.endswith((".nc", ".ncf"))]
lista = []

for file in os.listdir(dir_folder):
    
    if file.endswith(('.nc','.ncf')):
        sector_name = file.split('_')[0]
        
        print(f"Processing sector: {sector_name}")
        
        dir_data = os.path.join(dir_folder, file)
        data = xr.open_dataset(dir_data)
        lista.append(data)
        
