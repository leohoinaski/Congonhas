#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 11:13:56 2024

@author: leohoinaski
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import scipy.stats
import matplotlib.dates as mdates
import pandas as pd
import os
import modelEval_filter as mefil



def modelEvaluation (model,obs):
    trueNaN = np.isnan(obs[:,0])
    obs=obs[~trueNaN,:]
    model = model.iloc[~trueNaN,:].reset_index(drop=True)
    trueNaNmodel = np.isnan(model.iloc[:,0])
    obs=obs[~trueNaNmodel,:]
    model = model[~trueNaNmodel]
    bias = np.nanmean((model)-obs,axis=0)
    spearman, pv_spearman = scipy.stats.spearmanr(model,obs[:,0],axis=0)
    spearman=spearman[-1,0:model.shape[1]]
    pv_spearman=pv_spearman[-1,0:model.shape[1]]
    idBest = spearman.argmax(axis=0)
    stats=[]
    statsAll=[]
    stats.append(bias[idBest])
    stats.append(spearman[idBest])
    stats.append(pv_spearman[idBest])
    rmse = np.sqrt(((model - obs) ** 2).mean())
    stats.append(rmse[idBest])
    mae = np.nansum(abs(model - obs),axis=0)/obs.shape[0]
    stats.append(mae[idBest])
    stats.append(model.columns[idBest])
    stats.append(obs.shape[0])
    stats = pd.DataFrame(stats).transpose()
    stats.columns = ["Bias", "Spearman", "Spearman_pval","RMSE","MAE","Pixel",'n']
    statsAll.append(bias)
    statsAll.append(spearman)
    statsAll.append(pv_spearman)
    statsAll.append(rmse)
    mae = np.nansum(abs(model - obs),axis=0)/obs.shape[0]
    statsAll.append(mae)
    statsAll.append(model.columns)
    statsAll.append(np.matlib.repmat(obs.shape[0],obs.shape[1],1))
    statsAll = pd.DataFrame(statsAll).transpose()
    statsAll.columns = ["Bias", "Spearman", "Spearman_pval","RMSE","MAE","Pixel",'n']
    return stats,statsAll


def modelScatterplot(new_df,file,varType,unidade):
    #Scatterplots
    #select range os stations
    fig, ax = plt.subplots(figsize=(4,4))
    
    xy = np.vstack([new_df.iloc[:,0],new_df.iloc[:,1]])
    new_df = new_df.iloc[~np.any(np.isnan(xy), axis=0).transpose(),:]
    xy = xy[:,~np.any(np.isnan(xy), axis=0)]
    z = gaussian_kde(xy)(xy)
    new_df = new_df.iloc[~np.any(np.isnan(xy), axis=0).transpose(),:]
    ax.scatter(new_df.iloc[:,0],new_df.iloc[:,1],c=z,s=15,alpha=.5)
    
    ###calculate Spearman correlation using new_df
    corr, p_value = scipy.stats.spearmanr(new_df.iloc[:,0], new_df.iloc[:,1])
   
    ###insert text with Spearman correlation
    ax.annotate('ρ = {:.2f}'.format(corr), 
            xy=(0.70, 0.9), xycoords='axes fraction', 
            fontsize=8, ha='left', va='center')
    
    
    ax.set_xlabel(varType +' observada (' + unidade+')'+'\n'+file.split('.')[0].split('_')[-1]
                  ,fontsize=9)
    ax.set_ylabel(varType +' modelada (' + unidade+')'+' \nWRF' ,fontsize=9)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.set_xlim([new_df.min().min(),new_df.max().max()])
    ax.set_ylim([new_df.min().min(),new_df.max().max()])
    ax.set_aspect('equal')
    ax.plot([new_df.min().min(), new_df.max().max()],
             [new_df.min().min(), new_df.max().max()], 'k-', lw=1,dashes=[2, 2])
    ax.fill_between(np.linspace(new_df.min().min(),new_df.max().max(),new_df.shape[0]), 
                    np.linspace(new_df.min().min(),new_df.max().max(),new_df.shape[0])*0.5,
                    alpha=0.2,facecolor='gray',edgecolor=None)
    ax.fill_between(np.linspace(new_df.min().min(),new_df.max().max(),new_df.shape[0]),
                    np.linspace(new_df.max().max(),new_df.max().max(),new_df.shape[0]),
                    np.linspace(new_df.min().min(),new_df.max().max(),new_df.shape[0])*2,
                    alpha=0.2,facecolor='gray',edgecolor=None)
    fig.tight_layout()
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    return fig
    


folder = '/mnt/sdb1/Congonhas/data'
varType = ['Temperatura',
           'Pressão Atmosférica',
           'Umidade específica', 
           'Precipitação',
           'Umidade',
           'Direção do vento',
           'Velocidade do vento']

unidade = ['°C',
           'mBar',
           'kg/kg',
           'mm',
           '%',
           '°',
           'm/s']

varId = ['T2',
         'PATM',
         'Q2',
         'RAIN',
         'RH',
         'WDIR',
         'WIND']


outfolder = '/mnt/sdb1/Congonhas/E06/outputs'

os.makedirs(outfolder,exist_ok=True)
os.makedirs(outfolder+'/figures',exist_ok=True)
os.makedirs(outfolder+'/tables',exist_ok=True)


statsT = pd.DataFrame()
statsTall = pd.DataFrame()

for ii, varI in enumerate(varId): 
    files = [f for f in os.listdir(folder+'/'+varId[ii]) if os.path.isfile(os.path.join(folder+'/'+varId[ii], f))]
    for file in files:
        print(file)
        data = pd.read_csv(folder+'/'+varId[ii]+'/'+file,encoding ='latin-1')
        if data[~data['Valor'].isnull()].shape[0]/8760>0.5:
            model = data.iloc[:,data.columns.str.startswith('POINT')]
            obs = data['Valor']
            bolObs,statsOut=mefil.fixObs(obs)
            obs = obs[bolObs==False]
            model=model[bolObs==False]
            if obs.shape[0]>0:
                if model.shape[1]>0:
                    obs = np.matlib.repmat(obs,model.shape[1],1).transpose()
                    stats,statsAll = modelEvaluation (model,obs)
                    stats['Station'] = file.split('.')[0]
                    statsT = pd.concat([statsT,stats])
                    statsTall = pd.concat([statsTall,stats])
                    #statsTall=  statsTall.groupby(by="Pixel").mean().reset_index()
                    statsTall = statsTall.sort_values(by='Spearman', ascending=False)
                    statsTall = statsTall.drop_duplicates(subset='Pixel', keep="first")
    
                    #select point with best correlation
                    point_col = stats['Pixel'].values[0]
                    point_df = model[point_col]
                    
                    #create new df with point corresponding values
                    new_df = pd.DataFrame({'Valor': obs[:, stats.index[0]]})
                    new_df[point_col] = point_df
                    fig = modelScatterplot(new_df,file,varType[ii],unidade[ii])
                    fig.savefig(outfolder+'/figures/Scatter_'+varType[ii]+'_'+file+'.png',dpi=300)
 
statsTall.to_csv(outfolder+'/tables/stats.csv') 