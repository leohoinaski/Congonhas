# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 13:39:42 2024

O código analisa a meteorologia nos dias que tiveram rajadas de vento na cidade de Congonhas - MG

@author: 
    José Henrique Hess
    Researcher in LCQAr - Air Quality Control Laboratory
    Student in Sanitary and Environmental Engineering
    Federal University of Santa Catarina
"""

#%% Célula de importações

import numpy as np
import pandas as pd
import geopandas as gpd
import os

#%% Abrir os arquivos METEO.SFC

def read_METEO_SFC():
    '''

    Esta função busca ler arquivos METEO.SFC nos moldes da pág. 83 do seguinte link:
    https://www.epa.gov/sites/default/files/2020-09/documents/aermet_userguide.pdf
    Gera um dataframe com o nome das colunas

    Parameters
    ----------
    a : pd.DataFrame()
        Dataframe 

    Returns
    -------
    df: pd.Dataframe()
        Dataframe com os dados meteorológicos da estação
    '''
    

#%%

columns = [
    "Year", # Ano
    "Month", # Mês
    "Day", # Dia
    "JulianDay", # Dia Juliana 
    "Hour", # Hora
    "H", # Fluxo de calor sensível (W/m²)
    "u*", # Velocidade de fricção (m/s)
    "w*", # Escala de velocidade convectiva (m/s)
    "VPTG", # Gradiente de temperatura potencial vertical acima de Zic (K/m)
    "Zic", # Altura da camada limite por convecção (m)
    "Zim", # Altura da camada limite por mecânica (m)
    "L", # Altura Monin-Obukhov (m)
    "z0", # Comprimento da rugosidade da superfície (m)
    "B0", # Razão de Bowen
    "r", # Albedo
    "ws", # Velocidade do vento (m/s)
    "wd", # Direção do vento (grau)
    "zref", # Altura do anenometro (m)
    "temp", # Temperatura (K)
    "ztemp", # Altura do termômetro (m)
    "ipcode", # Tipo de chuva (0=sem, 11=líquida, 22=neve, 99=sem dado)
    "pamt", # Quantidade de chuva (mm/h)
    "rh", # Umidade relatica (%)
    "pres", # Pressão (mb)
    "ccvr", # Cobertura de nuvens (décimos)
    "WSADJ", # Ajuste da velocidade do vento e sinalizados da fonte de dados
]

arquivos_sfc = [arquivo for arquivo in os.listdir(r'C:\BolsaCongonhas\Git\Congonhas\E07\data\METEO') if arquivo.endswith('.sfc')]

nomes_estacao = pd.read_csv(r"C:\BolsaCongonhas\Git\Congonhas\E07\data\estacoes_wrf_xy.csv")

# Ler o arquivo como texto

for i in range(0, len(arquivos_sfc)):
    with open(r"C:\BolsaCongonhas\Git\Congonhas\E07\data\METEO/" + arquivos_sfc[i], 'r') as file:
        lines = file.readlines()
    
    # Processar o texto e transformar em uma lista de dicionários
    data = []
    for line in lines:
        if line != lines[0]:
            # Dividir os campos, supondo separador por espaço ou tabulação
            fields = line.strip().split()
            
            fields = [float(field) if field.replace('.', '', 1).isdigit() else field for field in fields]
            data.append(fields)
    
    # Criar o DataFrame
    df = pd.DataFrame(data, columns=columns)
    
    nome = nomes_estacao.loc[(nomes_estacao['wrf_nearest_lat'] == round(float(arquivos_sfc[i].split('_')[1]),6)) & (nomes_estacao['wrf_nearest_lon'] == round(float(arquivos_sfc[i].split('_')[2]),6)), 'Nome da estacao']
    
    print(nome.iloc[0])
    
    # Salvar como CSV
    df.to_csv('C:\BolsaCongonhas\Git\Congonhas\E07\data\METEO/'+ arquivos_sfc[i][:-4] +'_'+nome.iloc[0]+ '.csv', index=False)
