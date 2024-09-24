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
import contextily as cx
import matplotlib.patches as mpatches
import string
import matplotlib.dates as mdates
                                
def cityEmissTimeSeries(data,xlon,ylat,datesTime,shapeFilePath,folder,
                        pol,IBGE_CODE,source,cmap):
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
        fig.set_size_inches(19*cm, 12*cm)
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


def spatialFigure(data,xlon,ylat,legend,cmap,borderShapePath,folder,pol,
                  IBGE_CODE,source):
    
    fig, ax = plt.subplots(1,2)
    cm = 1/2.54  # centimeters in inches
    fig.set_size_inches(19*cm, 12*cm)
    #cmap = plt.get_cmap(cmap, 6)
    # bounds = np.array([np.percentile(data[data>0],1),
    #                    np.percentile(data[data>0],5),
    #                    np.percentile(data[data>0],10),
    #                     np.percentile(data[data>0],25),
    #                     np.percentile(data[data>0],50),
    #                     np.percentile(data[data>0],75),
    #                     np.percentile(data[data>0],90),
    #                     np.percentile(data[data>0],95),
    #                     np.percentile(data[data>0],99),
    #                     np.percentile(data[data>0],100)])
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    heatmap = ax[0].pcolor(xlon,ylat,data,cmap=cmap,
                        norm=mpl.colors.LogNorm(),alpha=0.6,
                        edgecolors=None)
    cbar = fig.colorbar(heatmap,fraction=0.04, pad=0.02,
                        #extend='both',
                        #ticks=bounds,
                        spacing='uniform',
                        orientation='horizontal',
                        #norm=norm,
                        ax=ax[0])

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
    br[br['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='black',linewidth=0.7,
                                                   ax=ax[0])
    #br.boundary.plot(edgecolor='black',linewidth=0.3,ax=ax,alpha=.6)
    ax[0].set_xlim([xlon.min(), xlon.max()])
    ax[0].set_ylim([ylat.min(), ylat.max()]) 
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    cx.add_basemap(ax[0], crs=br.crs, source=cx.providers.OpenStreetMap.Mapnik)
    
    s,cityMat,cityBuffer = tst.citiesBufferINdomain(xlon,ylat,br,IBGE_CODE,'CD_MUN')
    
    matData = data[:,:].copy()
    matData[np.isnan(cityMat)]=np.nan
    
    if (np.nanmin(matData)>0) and (np.nanmax(matData)):
        heatmap = ax[1].pcolor(xlon,ylat,matData,cmap=cmap,
                            norm=mpl.colors.LogNorm(),alpha=0.6,
                            edgecolors=None)
    else:
        heatmap = ax[1].pcolor(xlon,ylat,matData,cmap=cmap,
                            alpha=0.6,
                            edgecolors=None)
    print(np.nanmin(matData))
    print(np.nanmax(matData))
    cbar = fig.colorbar(heatmap,fraction=0.04, pad=0.02,
                         #extend='both',
                         #ticks=bounds,
                         spacing='uniform',
                         orientation='horizontal',
                         #norm=norm,
                         ax=ax[1])
    
    cbar.ax.tick_params(rotation=30)
    #tick_locator = mpl.ticker.MaxNLocator(nbins=5)
    #cbar.locator = tick_locator
    #cbar.ax.set_xscale('log')
    #cbar.update_ticks()
    
    cbar.ax.set_xlabel(legend, rotation=0,fontsize=6)
    cbar.ax.get_xaxis().labelpad = 0
    cbar.ax.tick_params(labelsize=6)
    #cbar.ax.locator_params(axis='both',nbins=5)
    cbar.ax.minorticks_off()
    
    br[br['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='blacK',
                                                   linewidth=0.7,ax=ax[1])
    #br.boundary.plot(edgecolor='black',linewidth=0.3,ax=ax,alpha=.6)
    ax[1].set_xlim([br[br['CD_MUN']==str(IBGE_CODE)].bounds.minx.values,
                    br[br['CD_MUN']==str(IBGE_CODE)].bounds.maxx.values])
    ax[1].set_ylim([br[br['CD_MUN']==str(IBGE_CODE)].bounds.miny.values,
                    br[br['CD_MUN']==str(IBGE_CODE)].bounds.maxy.values]) 
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[1].set_anchor('C')
    cx.add_basemap(ax[1], crs=br.crs, source=cx.providers.OpenStreetMap.Mapnik)
    for n, axs in enumerate(ax):
        axs.text(0.01, 0.95, string.ascii_lowercase[n]+')', transform=axs.transAxes, 
                size=10, weight='normal')
    fig.tight_layout()
    fig.savefig(folder+'/spatialFigure_'+pol+'_'+source+'.png',
                format="png",bbox_inches='tight')
    
def maxPixelFigure(data,xlon,ylat,legend,cmap,borderShapePath,folder,pol,IBGE_CODE,
                   source,aveTime):


    fig, ax = plt.subplots(1,2)
    cm = 1/2.54  # centimeters in inches
    fig.set_size_inches(19*cm, 12*cm)
    #cmap = plt.get_cmap(cmap, 6)
    
    if aveTime =='Month': 
        cmap = plt.cm.get_cmap(cmap, 13) 
        bound = np.arange(1,13)
        title='Mês'
    elif aveTime =='DayOfWeek': 
        cmap = plt.cm.get_cmap(cmap, 8) 
        bound = np.arange(1,8)
        title='Dia da semana'
    elif aveTime =='Hour':
        cmap = plt.cm.get_cmap(cmap, 24) 
        bound = np.arange(0,24)
        title='Hora'
    else:
        cmap = plt.cm.get_cmap(cmap, 13) 
        bound = np.arange(0,13)    
            
    #cmap.set_bad(color='white')
    cmap.set_over('black')
    heatmap = ax[0].pcolormesh(xlon,ylat,data, cmap=cmap,alpha=0.5,
                            vmin=bound.min(),vmax=bound.max())
    
    ax[1].legend([mpatches.Patch(color=cmap(b)) for b in bound],
               ['{} '.format(bound[i-1]) for i in bound], ncol=4,
               prop={'size': 6},frameon=False,
               loc='upper center', bbox_to_anchor=(0.5, -0.05),
               title=title,columnspacing=0.6)
    
    br = gpd.read_file(borderShapePath)
    br[br['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='blacK',
                                                   linewidth=0.8,ax=ax[0])
    #br.boundary.plot(edgecolor='black',linewidth=0.3,ax=ax,alpha=.6)
    ax[0].set_xlim([xlon.min(), xlon.max()])
    ax[0].set_ylim([ylat.min(), ylat.max()]) 
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    cx.add_basemap(ax[0], crs=br.crs, source=cx.providers.OpenStreetMap.Mapnik)

    heatmap = ax[1].pcolormesh(xlon,ylat,data, cmap=cmap,alpha=0.5,
                            vmin=bound.min(),vmax=bound.max())
    
    br = gpd.read_file(borderShapePath)
    br[br['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='blacK',
                                                   linewidth=0.8,ax=ax[1])
    #br.boundary.plot(edgecolor='black',linewidth=0.3,ax=ax,alpha=.6)
    ax[1].set_xlim([br[br['CD_MUN']==str(IBGE_CODE)].bounds.minx.values,
                    br[br['CD_MUN']==str(IBGE_CODE)].bounds.maxx.values])
    ax[1].set_ylim([br[br['CD_MUN']==str(IBGE_CODE)].bounds.miny.values,
                    br[br['CD_MUN']==str(IBGE_CODE)].bounds.maxy.values]) 
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    cx.add_basemap(ax[1], crs=br.crs, source=cx.providers.OpenStreetMap.Mapnik)
    ax[1].set_anchor('C')
    fig.tight_layout()
    for n, axs in enumerate(ax):
        axs.text(0.01, 0.95, string.ascii_lowercase[n]+')', transform=axs.transAxes, 
                size=10, weight='normal')
        
    fig.savefig(folder+'/maxSpatialFigure_'+aveTime+'_'+pol+'_'+source+'.png',
                format="png",bbox_inches='tight',dpi=300)
    return fig


def maxPixelFigureAll(data,xlon,ylat,legend,SOURCES,borderShapePath,folder,pol,IBGE_CODE,
                   source,aveTime,sources):
    import matplotlib.patches as mpatches
    import ismember
    fig, ax = plt.subplots(1,2)
    cm = 1/2.54  # centimeters in inches
    fig.set_size_inches(19*cm, 12*cm)
    #cmap = plt.get_cmap(cmap, 6)
    
    lia, loc = ismember.ismember(sources,SOURCES['file'])
    cmap = mpl.colors.ListedColormap(np.array(SOURCES['color'])[loc])
    bound = np.arange(0,len(sources))    
            
    #cmap.set_bad(color='white')
    #cmap.set_over('white')
    heatmap = ax[0].pcolormesh(xlon,ylat,data, cmap=cmap,vmin=0, vmax=len(sources), 
                               alpha=0.5)
    
    ax[1].legend([mpatches.Patch(color=cmap(b)) for b in bound],
               ['{} '.format(np.array(SOURCES['source'])[loc][i]) for i in bound], ncol=1,
               prop={'size': 6},frameon=False,
               loc='center left', bbox_to_anchor=(1, 0.5))
    
    br = gpd.read_file(borderShapePath)
    br[br['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='black',
                                                   linewidth=0.8,ax=ax[0])
    #br.boundary.plot(edgecolor='black',linewidth=0.3,ax=ax,alpha=.6)
    ax[0].set_xlim([xlon.min(), xlon.max()])
    ax[0].set_ylim([ylat.min(), ylat.max()]) 
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    cx.add_basemap(ax[0], crs=br.crs, source=cx.providers.OpenStreetMap.Mapnik)

    heatmap = ax[1].pcolormesh(xlon,ylat,data,vmin=0,vmax=len(sources),
                               cmap=cmap,alpha=0.5)

    br = gpd.read_file(borderShapePath)
    br[br['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='black',
                                                   linewidth=0.8,ax=ax[1])
    #br.boundary.plot(edgecolor='black',linewidth=0.3,ax=ax,alpha=.6)
    ax[1].set_xlim([br[br['CD_MUN']==str(IBGE_CODE)].bounds.minx.values,
                    br[br['CD_MUN']==str(IBGE_CODE)].bounds.maxx.values])
    ax[1].set_ylim([br[br['CD_MUN']==str(IBGE_CODE)].bounds.miny.values,
                    br[br['CD_MUN']==str(IBGE_CODE)].bounds.maxy.values]) 
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    cx.add_basemap(ax[1], crs=br.crs, source=cx.providers.OpenStreetMap.Mapnik)
    ax[1].set_anchor('C')
    fig.tight_layout()
    fig.savefig(folder+'/maxSpatialFigure_'+aveTime+'_'+pol+'_'+source+'.png',
                format="png",bbox_inches='tight',dpi=300)
    return fig

def maxPixelFigureAllbyPeriod(data,xlon,ylat,legend,cmap,borderShapePath,
                              folder,pol,IBGE_CODE,
                   source,aveTime,sources):
    import math
    import matplotlib.patches as mpatches
    nrows = math.ceil(data.shape[0]/2)
    
    fig, ax = plt.subplots(1,2)
    cm = 1/2.54  # centimeters in inches
    fig.set_size_inches(19*cm, 12*cm)
    #cmap = plt.get_cmap(cmap, 6)

    cmap = plt.cm.get_cmap(cmap, len(sources)) 
    bound = np.arange(0,len(sources))    
            
    #cmap.set_bad(color='white')
    cmap.set_over('white')
    heatmap = ax[0].pcolormesh(xlon,ylat,data, cmap=cmap,alpha=0.5,
                            vmin=1,vmax=bound.max())
    
    ax[1].legend([mpatches.Patch(color=cmap(b)) for b in bound],
               ['{} '.format(sources[i]) for i in bound], ncol=1,
               prop={'size': 6},frameon=False,
               loc='center left', bbox_to_anchor=(1, 0.5))
    
    br = gpd.read_file(borderShapePath)
    br[br['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='blacK',
                                                   linewidth=0.8,ax=ax[0])
    #br.boundary.plot(edgecolor='black',linewidth=0.3,ax=ax,alpha=.6)
    ax[0].set_xlim([xlon.min(), xlon.max()])
    ax[0].set_ylim([ylat.min(), ylat.max()]) 
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    cx.add_basemap(ax[0], crs=br.crs, source=cx.providers.OpenStreetMap.Mapnik)

    heatmap = ax[1].pcolormesh(xlon,ylat,data, cmap=cmap,alpha=0.5,
                            vmin=1,vmax=bound.max())

    
    br = gpd.read_file(borderShapePath)
    br[br['CD_MUN']==str(IBGE_CODE)].boundary.plot(edgecolor='blacK',
                                                   linewidth=0.8,ax=ax[1])
    #br.boundary.plot(edgecolor='black',linewidth=0.3,ax=ax,alpha=.6)
    ax[1].set_xlim([br[br['CD_MUN']==str(IBGE_CODE)].bounds.minx.values,
                    br[br['CD_MUN']==str(IBGE_CODE)].bounds.maxx.values])
    ax[1].set_ylim([br[br['CD_MUN']==str(IBGE_CODE)].bounds.miny.values,
                    br[br['CD_MUN']==str(IBGE_CODE)].bounds.maxy.values]) 
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    cx.add_basemap(ax[1], crs=br.crs, source=cx.providers.OpenStreetMap.Mapnik)
    ax[1].set_anchor('C')
    fig.tight_layout()
    fig.savefig(folder+'/maxSpatialFigure_'+aveTime+'_'+pol+'_'+source+'.png',
                format="png",bbox_inches='tight',dpi=300)
    return fig