#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 18:46:13 2021

@author: doan
"""
import pandas as pd
import glob, os, sys
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap, PowerNorm
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import numpy as np
from sklearn.linear_model import LinearRegression
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

def grid(ax,st):
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = False
    gl.ylines = False
    gl.xlocator = mticker.FixedLocator(np.arange(0,361,st))
    gl.ylocator = mticker.FixedLocator(np.arange(-90.,90,st))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':9, 'color':'gray'}
    gl.ylabel_style = {'size':9, 'rotation':90, 'va':'center', 'color':'gray'}
    return ax



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import cartopy.feature as cfeature
#land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
#                                        edgecolor='face',
#                                        facecolor=cfeature.COLORS['land'])

bodr = cfeature.NaturalEarthFeature(category='cultural', 
    name='admin_0_boundary_lines_land', scale='10m', facecolor='none', alpha=0.7)


from matplotlib.image import imread


def plot_basemap(hr=False, size=(10, 6)):
    
    proj =  ccrs.PlateCarree() #ccrs.Robinson(central_longitude=145)

    fig = plt.figure(figsize=size)
    ax = plt.axes(projection= proj )
    #ax.set_extent([123,150,25, 49])
    ax.coastlines(resolution='10m',lw=.5)
    #ax.gridlines()
    ax.stock_img()
    
    fname = os.path.join('station_info/NE1_50M_SR_W', 'NE1_50M_SR_W.tif')
    if hr:
        ax.imshow(imread(fname), origin='upper', transform=proj, 
              extent=[-180, 180, -90, 90])
    #ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5,lw=.5)
    ax.add_feature(bodr, linestyle='-', edgecolor='k', alpha=.5, lw=.5)
    ax.outline_patch.set_linewidth(.5)
    
    return fig, ax
    grid(ax,1)
    
    
    
    

if __name__ == "__main__":
    fig, ax = plot_basemap()
    
    
    sys.exit()
    df = pd.read_csv('Amedas_list.csv', index_col=0).set_index('station_id')
    df = df[~df.index.duplicated()]

    #df['val'] = np.log2(x.values)


    mm = ['^', 'v']
    cc = ['r', 'b']
    
    
    #for icond, cond in enumerate([df['val']>0, df['val']<=0]):
    #    d1 = df[cond]
        
    x = df.longitude.values
    y = df.latitude.values
    #c = d1['val'].values
    #a = (1 + d1['val'].values ) * 100
        
    plt.scatter(x, y, c='r', s = 1, 
                    alpha=1, transform=ccrs.Geodetic(),marker='o')
        
    


































