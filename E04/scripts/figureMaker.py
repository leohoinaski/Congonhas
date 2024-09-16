#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:03:09 2024

@author: leohoinaski
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import geopandas as gpd
import temporalStatistics as tst


                                
def cityEmissTimeSeries(data,xlon,ylat,datesTime,shapeFilePath,folder,
                        pol,IBGE_CODE,source,cmap):
    
    import matplotlib.dates as mdates
    
    cityData,cityDataPoints,cityDataFrame,matData,gdf = tst.timeSeriesByCity(
        data,xlon,ylat,datesTime,shapeFilePath,folder,pol,IBGE_CODE,source)

    if len(matData.shape)==4:
        aveFigData= np.nanmean(matData,axis=0)[0,:,:]
    else:
        aveFigData= np.nanmean(matData,axis=0)

    if (np.nanmax(aveFigData)>0):

        cityArea=gdf[gdf['CD_MUN']==str(IBGE_CODE)]
        #cmap = plt.get_cmap(cmap,5)    
        fig, ax = plt.subplots(1,2,gridspec_kw={'width_ratios': [1, 3]})
        cm = 1/2.54  # centimeters in inches
        fig.set_size_inches(14*cm, 7*cm)
        cmap.set_under('white')

        heatmap = ax[0].pcolor(xlon,ylat,aveFigData,cmap=cmap)
        cbar = fig.colorbar(heatmap,fraction=0.04, pad=0.02,

                            #extend='both',
                            spacing='uniform',
                            orientation='horizontal',
                            ax=ax[0])

        cbar.ax.tick_params(rotation=30)
        #tick_locator = mpl.ticker.MaxNLocator(nbins=5)
        #cbar.locator = tick_locator
        #cbar.ax.set_xscale('log')
        #cbar.update_ticks()
        #cbar.ax.locator_params(axis='both',nbins=5)
        #cbar.ax.set_yscale('log')
        #cbar.update_ticks()
        #cbar.ax.set_xticklabels(['{:.1e}'.format(x) for x in bounds],rotation=30)
        cbar.ax.set_xlabel(cityArea['NM_MUN'].to_string(index=False)+'\nAverage', rotation=0,fontsize=6)
        cbar.ax.get_xaxis().labelpad = 5
        cbar.ax.tick_params(labelsize=6) 


        ax[0].set_xlim([cityArea.boundary.total_bounds[0],cityArea.boundary.total_bounds[2]])
        ax[0].set_ylim([cityArea.boundary.total_bounds[1],cityArea.boundary.total_bounds[3]])
        ax[0].set_frame_on(False)
        ax[0].set_xticks([])
        ax[0].set_yticks([])
        cityArea.boundary.plot(edgecolor='black',linewidth=0.5,ax=ax[0])
        gdf[gdf['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='gray',linewidth=0.3,ax=ax[0])

        ax[1].fill_between(cityDataFrame.mean(axis=1).index,cityDataFrame.max(axis=1), cityDataFrame.min(axis=1),
                         color=cmap(0.8),       # The outline color
                         facecolor=cmap(0.8),
                         edgecolor=None,
                         alpha=0.2,label='Min-Max')          # Transparency of the fill
        ax[1].plot(cityDataFrame.mean(axis=1).index,cityDataFrame.mean(axis=1),
                   color=cmap(0.8),linewidth=1,label='Average')
        ax[1].xaxis.set_tick_params(labelsize=6)
        ax[1].yaxis.set_tick_params(labelsize=6)
        ax[1].set_ylim([np.nanmin(matData)*0.95,np.nanmax(matData)*1.05])
        ax[1].set_xlim([cityDataFrame.index.min(),cityDataFrame.index.max()])
        ax[1].xaxis.set_major_locator(mdates.MonthLocator(interval=1))
        # set formatter
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
        for label in ax[1].get_xticklabels(which='major'):
            label.set(rotation=30, horizontalalignment='right')
        ax[1].legend(prop={'size': 6})
        ax[1].set_ylabel(cityArea['NM_MUN'].to_string(index=False)+'\n'+source,fontsize=6)
        fig.tight_layout()
        fig.savefig(folder+'/cityTimeSeries_'+pol+'_'+source+'.png', format="png",
                   bbox_inches='tight')
        return matData.shape


def spatialFigure(data,xlon,ylat,legend,cmap,borderShapePath,folder,pol,IBGE_CODE,source):
    fig, ax = plt.subplots()
    cm = 1/2.54  # centimeters in inches
    fig.set_size_inches(15*cm, 10*cm)
    #cmap = plt.get_cmap(cmap, 6)
    bounds = np.array([np.percentile(data[data>0],1),
                       np.percentile(data[data>0],5),
                       np.percentile(data[data>0],10),
                        np.percentile(data[data>0],25),
                        np.percentile(data[data>0],50),
                        np.percentile(data[data>0],75),
                        np.percentile(data[data>0],90),
                        np.percentile(data[data>0],95),
                        np.percentile(data[data>0],99),
                        np.percentile(data[data>0],100)])
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    heatmap = ax.pcolor(xlon,ylat,data,cmap=cmap,norm=norm)
    cbar = fig.colorbar(heatmap,fraction=0.04, pad=0.02,
                        #extend='both',
                        ticks=bounds,
                        spacing='uniform',
                        orientation='horizontal',
                        norm=norm,
                        ax=ax)

    cbar.ax.tick_params(rotation=30)
    #tick_locator = mpl.ticker.MaxNLocator(nbins=5)
    #cbar.locator = tick_locator
    #cbar.ax.set_xscale('log')
    #cbar.update_ticks()
    
    cbar.ax.set_xlabel(legend, rotation=0,fontsize=6)
    cbar.ax.get_xaxis().labelpad = 2
    cbar.ax.tick_params(labelsize=6)
    #cbar.ax.locator_params(axis='both',nbins=5)
    cbar.ax.minorticks_off()
    br = gpd.read_file(borderShapePath)
    br[br['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='blacK',linewidth=0.5,ax=ax)
    br.boundary.plot(edgecolor='black',linewidth=0.3,ax=ax,alpha=.6)
    ax.set_xlim([xlon.min(), xlon.max()])
    ax.set_ylim([ylat.min(), ylat.max()]) 
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(folder+'/spatialFigure_'+pol+'_'+source+'.png',
                format="png",bbox_inches='tight')
    
def maxPixelFigure(data,xlon,ylat,legend,cmap,borderShapePath,folder,pol,IBGE_CODE,
                   source,aveTime):
    fig, ax = plt.subplots()
    cm = 1/2.54  # centimeters in inches
    fig.set_size_inches(15*cm, 10*cm)
    #cmap = plt.get_cmap(cmap, 6)
    heatmap = ax.pcolor(xlon,ylat,data)
    cbar = fig.colorbar(heatmap,fraction=0.04, pad=0.02,
                        #extend='both',
                        spacing='uniform',
                        orientation='horizontal',
                        ax=ax)
    cbar.ax.tick_params(rotation=30)
    #tick_locator = mpl.ticker.MaxNLocator(nbins=5)
    #cbar.locator = tick_locator
    #cbar.ax.set_xscale('log')
    #cbar.update_ticks()
    
    cbar.ax.set_xlabel(legend, rotation=0,fontsize=6)
    cbar.ax.get_xaxis().labelpad = 2
    cbar.ax.tick_params(labelsize=6)
    #cbar.ax.locator_params(axis='both',nbins=5)
    cbar.ax.minorticks_off()
    br = gpd.read_file(borderShapePath)
    br[br['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='blacK',linewidth=0.5,ax=ax)
    br.boundary.plot(edgecolor='black',linewidth=0.3,ax=ax,alpha=.6)
    ax.set_xlim([xlon.min(), xlon.max()])
    ax.set_ylim([ylat.min(), ylat.max()]) 
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(folder+'/maxSpatialFigure_'+aveTime+'_'+pol+'_'+source+'.png',
                format="png",bbox_inches='tight')