#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 23:28:24 2023

@author: doan
"""

# download sounding data
    
    
###
import pandas as pd
import os, glob, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap, PowerNorm
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import wget



if True: # download by country

    ctry = 'VM'
    ctry = 'NG'
    
    # https://www.ncei.noaa.gov/pub/data/igra/
    #df = pd.read_csv('igra2-station-list.txt', error_bad_lines=False, index_col=0, header=None, delim_whitespace=True)
    
    ll = open('igra2-station-list.txt').read().split('\n')[:-1]
    l1 = [ [l[:11], l[11:20], l[20:30], l[30:37], l[37:72], l[72:76], l[76:82 ], l[82:]] for l in ll]
    l4 = [ [l2.strip() for l2 in l2] for l2 in l1]
    df = pd.DataFrame(l4, columns=['IND', 'LAT', 'LON', 'ELEV', 'NAME', 'STY', 'ENY', 'X'] )
    

    df = df.set_index('IND')
    for c in ['LAT', 'LON', 'ELEV', 'STY', 'ENY', 'X']:
        df.loc[:,c] = pd.to_numeric( df.loc[:,c]) #, errors='coerse')
    
    lat, lon = 42.6977, 23.3219
    lat, lon = 6.5244, 3.3792
    b = df[ (df.LAT < lat+1) & (df.LAT > lat-1)  & 
            (df.LON < lon+1) & (df.LON > lon-1)  
            ]
    
    url0 = 'https://www1.ncdc.noaa.gov/pub/data/igra/data/data-por/' 
    url0 = 'ftp://ftp.ncdc.noaa.gov/pub/data/igra/data/data-por/'
    
    b['f'] = b.index+'-data.txt.zip'


    for i, r in b[:].iterrows():
        odir = 'download/'
        if not os.path.exists(odir): os.makedirs(odir)
        url = url0 + r['f']
        print(url)
        try:
            os.system('wget '+url)
            #filename = wget.download(url, out=odir)
        except:
            print('no exisit')
            
            
            
            
            
            
            
            
            
            
            
            