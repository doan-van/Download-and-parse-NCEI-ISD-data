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
import xarray as xr


# proccess station informations

# ========
# convert txt to csv for list of stations
if True:
    ll = open('data/isd-history.txt').readlines()
    dd = np.array([ [l[:6], l[6:13], l[13:43], l[43:48], l[48:51], l[51:57], l[57:65], l[65:74], l[74:82], l[82:91], l[91:99] ]  for l in ll])
    df = pd.DataFrame(dd[22:], columns=[ c.strip() for c in dd[20] ] )
    for c in ['LAT', 'LON', 'ELEV(M)']: df.loc[:,c] =  pd.to_numeric(df[c], errors='coerce')
    for c in ['BEGIN', 'END']: df.loc[:,c] =  pd.to_datetime(df[c]) #, errors='coerce')

    df.to_csv('data/list-isd-history.csv')
    
    #do = df[ (df.BEGIN < '2010-01-01') & (df.END > '2021-05-31') ]
    #do.to_csv('list_station_10y.csv')

#=======










# Plot stations 
if False:
    df = pd.read_csv('data/list-isd-history.csv')
    out_filename='map_all_station.png'
    #proj = ccrs.PlateCarree(central_longitude=145)
    #proj._threshold /= 10
    proj = ccrs.Robinson(central_longitude=145)
    
    plt.figure(figsize=(10, 6))
    ax = plt.axes(projection= proj )
    ax.coastlines(resolution='10m',lw=.1)
    #ax.gridlines()
    ax.stock_img()
    
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5,lw=.3)
    ax.outline_patch.set_linewidth(.5)

    ax.scatter(df['LON'],df['LAT'],color='r',transform=ccrs.Geodetic(),s=.1)
    plt.savefig(out_filename, format='png', bbox_inches='tight', dpi = 200)
    
    
    
    
    
# find stations within 100 km from cities
# 
if False:
    ci = pd.read_excel('../WORLD-CITIES/out/cities_list_llcoast_check_cities.xlsx',index_col=0,header=0) #.sort_values(2020, ascending=False)
    df = pd.read_csv('list-isd-history.csv', index_col=0)  
    
    
    for i, r in ci[:].iterrows():
        print(r['Urban Agglomeration'])
        la, lo = r['Latitude'], r['Longitude']
    
        bu = 2
        d1 = df[ (df.LAT < la+bu) & (df.LAT > la-bu)  & (df.LON < lo+bu) & (df.LON > lo-bu)  ]    
        y, x = d1.LAT.values, d1.LON.values
        y = np.where( (y > 90) | (y < -90), np.NaN, y)
        x = np.where( (x > 180) | (x < -180), np.NaN, x)
    
        p = 0.017453292519943295
        b = 0.5 - np.cos((la-y)*p)/2 + np.cos(la*p)*np.cos(y*p) * (1-np.cos((x-lo)*p)) / 2
        a = 12742 * np.arcsin(np.sqrt(b))    
    
        d2 = d1[a < 100]
        ci.loc[i, 'IDS'] = ','.join((d2.USAF.astype(str) + d2.WBAN.astype(str)).values)
    
    ci.to_csv('data/world_cities_list_with_ids_station.csv')
    
    
   
    
    
    
    
#=====================================================   
# calculate distance to the coastline for each station
if False:
    
    df = pd.read_csv('list-isd-history.csv', index_col=0)    
    dc = xr.open_dataset('../WORLD-CITIES/d2cl.nc')['d2cl']
    slats = []
    slons = []
    d2s = []
        
    for i, r in df[:].iterrows():
         
        print(r['STATION NAME'])
        
        la, lo = r['LAT'], r['LON']
        
        
        if np.isnan(la) | np.isnan(lo) | (la > 90) | (la < -90) | (lo > 180) | (lo < -180) :
            slats.append(np.NaN)
            slons.append(np.NaN)
            d2s.append( np.NaN )               
            continue
        
        d2sea = dc.sel(lat=la, lon = lo, method='nearest').values
        # if located over ocean then get box 1 deg
        bu = np.ceil( - d2sea / 100)*1.5
        if bu <= 0: bu = 1
        
        print(d2sea, bu)
        d1 = dc.sel(lat = slice(la-bu, la+bu), lon=slice(lo-bu, lo+bu) )
        
        print(d2sea)
        
        inan = (d1 > 0).values.flatten()
        if (~inan).all(): # all ocean
            slats.append(lo)
            slons.append(la)
            d2s.append( 0. )   
            continue
            
        
        lon, lat = d1.lon.values, d1.lat.values
        x0, y0 = np.meshgrid(lon, lat)
        
        x1, y1 = x0.flatten(), y0.flatten()
        x, y = np.where(inan, x1, np.nan), np.where(inan, y1, np.nan)
        
        p = 0.017453292519943295
        b = 0.5 - np.cos((la-y)*p)/2 + np.cos(la*p)*np.cos(y*p) * (1-np.cos((x-lo)*p)) / 2
        a = 12742 * np.arcsin(np.sqrt(b))
        
        distance = np.nanmin(a, axis=None)
        ind = np.unravel_index(np.nanargmin(a, axis=None), x0.shape)
        
        print(ind)
        slat, slon = lat[ind[0]], lon[ind[1]]
        slats.append(slat)
        slons.append(slon)
        d2s.append( distance )    
    
    
        if False:
            
            sys.path.append('../../')
            from common import *
            fig, ax = plot_basemap(hr=True)
            ax.set_extent([d1.lon.min(),d1.lon.max(),d1.lat.min(), d1.lat.max()])
            plt.scatter(lo, la, color='r',s=40, alpha=1, marker='s')
            plt.scatter(slon, slat, color='k',s=20, alpha=1)
            plt.plot([lo,slon], [la,slat])
            
            
            
    df['lat_cl'] = slats
    df['lon_cl'] = slons
    df['d2s'] = d2s
    df.to_csv('data/list-isd-history_llcoast_all.csv')
#=======================================
# need to do it again ==================    
    
    
    
    
    
# nead to calculate angle to coastline
    
if False:
    df = pd.read_csv('data/list-isd-history_llcoast_all.csv', index_col=0)    

    
    y1, x1 = df.lat_cl.values, df.lon_cl.values
    y2, x2 = df.LAT.values, df.LON.values
    
    if True:
        #x1, y1, x2, y2 = slon, slat, lon2, lat2
        dx = x2 -x1
        phi0 = np.where(dx > 0, 90, 270)
                
        #--- cal direction ---#
        term1 = np.sin( np.deg2rad( dx ) )
        term1 = np.where(term1 == 0, 1.0*10**(-5), term1)
                
        term2 = np.cos( np.deg2rad( y1 ) ) * np.tan(np.deg2rad( y2 ) )\
                        - np.sin( np.deg2rad( y1 ) ) * np.cos(np.deg2rad( dx ) )
                
        phi = phi0 - np.rad2deg ( np.arctan( term2 / term1 ) )
                
    df['angle'] = phi
    df.to_csv('data/list-isd-history_llcoast.csv')
    
    
    
    
    
    
# Plot
if False:
    df = pd.read_csv('data/list-isd-history_llcoast.csv', index_col=0)   
    df = df[df.d2s > 20]
    # Plot stations 
    if False:
        out_filename='map_all_station_c.png'
        #proj = ccrs.PlateCarree(central_longitude=145)
        #proj._threshold /= 10
        proj = ccrs.Robinson(central_longitude=145)
        
        plt.figure(figsize=(10, 6))
        ax = plt.axes(projection= proj )
        ax.coastlines(resolution='10m',lw=.1)
        #ax.gridlines()
        ax.stock_img()
        
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5,lw=.3)
        ax.outline_patch.set_linewidth(.5)
    
        ax.scatter(df['LON'],df['LAT'],color='r',transform=ccrs.Geodetic(),s=.1)
        plt.savefig(out_filename, format='png', bbox_inches='tight', dpi = 200)
    
    
    if True:
        
        import random
        ir = random.randint(0, len(df)-1)
        d = df.iloc[ir]
        fig = plt.figure(figsize=(4.5,4.5))
        ax = plt.axes([.1,.2,.8,.75], projection=ccrs.PlateCarree())
        la, lo = d['LAT'], d['LON']
        sla, slo = d.lat_cl, d.lon_cl
        bu = 1.5
        ax.set_extent([lo-bu,lo+bu,la-bu,la+bu])
        ax.coastlines(resolution='10m',lw=.5)
        #ax.stock_img()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5,lw=.3)
        
        ax.scatter(lo,la, c='r', marker='s',s=20 )
        ax.scatter(slo,sla, c='g' )
        ax.plot([lo,slo], [la,sla], lw=1,c='b')
        plt.axhline(y=la, color='k', linestyle='-', lw=.5)
        plt.axvline(x=lo, color='k', linestyle='-', lw=0.5)

        #scale_bar(ax, ccrs.PlateCarree(), 50) 
        #add_grid(ax,1)
            
        ax.text(.99,.99,'$\Phi$%.0f'%d.angle, fontsize=12,
                        ha='right', va='top',transform=ax.transAxes)



    
    
    
    
    
    
    
    
    
    
    
    