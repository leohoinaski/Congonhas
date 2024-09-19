#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 09:59:08 2024

@author: leohoinaski
"""
import os 
import numpy as np


def listDatasets(dataPath,GDNAM,pol):
    
    prefixed = [filename for filename in os.listdir(dataPath) if 
                filename.startswith("BRAIN_BASEMIS_"+GDNAM)]
    matching = [s for s in prefixed if pol in s]
    return matching

def all_equal(sequence):
    return len(set(sequence)) == 1
        
def checkMatEquals(monthlySum):
    if len(monthlySum.shape)==4:
        test = np.empty((monthlySum.shape[2],monthlySum.shape[3])).astype(bool)
        for ii in range(0,monthlySum.shape[2]):
            for jj in range(0,monthlySum.shape[3]):
                test[ii,jj] = all_equal(monthlySum[:,0,ii,jj])
                
    if len(monthlySum.shape)==3:
        test = np.empty((monthlySum.shape[1],monthlySum.shape[2])).astype(bool)
        for ii in range(0,monthlySum.shape[1]):
            for jj in range(0,monthlySum.shape[2]):
                test[ii,jj] = all_equal(monthlySum[:,ii,jj])
    return test

def agrmaxArray(data):
    if len(data.shape)==4:
        test = np.empty((data.shape[2],data.shape[3]))
        for ii in range(0,data.shape[2]):
            for jj in range(0,data.shape[3]):
                test[ii,jj] = np.nanargmax(data[:,0,ii,jj])
                
    if len(data.shape)==3:
        test = np.empty((data.shape[1],data.shape[2]))
        for ii in range(0,data.shape[1]):
            for jj in range(0,data.shape[2]):
                test[ii,jj] = np.nanargmax(data[:,ii,jj])
    return test