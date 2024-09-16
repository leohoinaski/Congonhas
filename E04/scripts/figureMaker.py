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