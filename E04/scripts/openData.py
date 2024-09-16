#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 13:39:11 2024



@author: leohoinaski
"""



import os
import netCDF4 as nc


NO2 = {
  "Pollutant": "$NO_{2}$",
  "Unit": '$\u03BCg.m^{-3}$',
  "conv": 1880,
  "tag":'NO2',
  #"Criteria": 260, # 260, 240, 220, 200
}

CO = {
  "Pollutant": "CO",
  "Unit": 'ppb',
  "conv": 1000, # Conversão de ppm para ppb
  "tag":'CO',
}

O3 = {
  "Pollutant": "$O_{3}$",
  "Unit": 'ppm',
  "conv": 1,
  "tag":'O3'
}

SO2 = {
  "Pollutant": "$SO_{2}$",
  "Unit": '$\u03BCg.m^{-3}$',
  "conv": 2620,
  "tag":'SO2'
}

PM10 = {
  "Pollutant": "$PM_{10}$",
  "Unit": '$\u03BCg.m^{-3}$',
  "conv": 1,
  "tag":'PM10',
}

PM25 = {
  "Pollutant": "$PM_{2.5}$",
  "Unit": '$\u03BCg.m^{-3}$',
  "conv": 1,
  "tag":'PM25',
}


pollutants = [CO]
pol = 'PM10'

GDNAM = 'Con_3km'


cwdPath = os.path.abspath(os.getcwd())
rootPath = os.path.dirname(os.path.dirname(cwdPath))
dataPath = rootPath + '/data/'+GDNAM

prefixed = [filename for filename in os.listdir(dataPath) if 
            filename.startswith("BRAIN_BASEMIS_"+GDNAM)]

matching = [s for s in prefixed if pol in s]





