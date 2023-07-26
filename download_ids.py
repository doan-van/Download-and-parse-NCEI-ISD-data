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
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import cartopy.feature as cfeature
from matplotlib.image import imread

#land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
#                                        edgecolor='face',
#                                        facecolor=cfeature.COLORS['land'])


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
#==============================================================================




#==============================================================================
#==============================================================================
#==============================================================================
def plot_basemap(hr=False, size=(10, 6)):
    bodr = cfeature.NaturalEarthFeature(category='cultural', 
        name='admin_0_boundary_lines_land', scale='10m', facecolor='none', alpha=0.7)
    
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
    grid(ax,1)
    return fig, ax
#==============================================================================
#==============================================================================
#==============================================================================




# https://www.ncei.noaa.gov/data/global-summary-of-the-day/doc/



#==============================================================================
#==============================================================================
#==============================================================================
# proccess station informations
def process_data(fil, odir, ofile=None, rm=True):
    # process original data
    d = pd.read_csv(fil, index_col=1, parse_dates = True)
    print(d)
    stn = fil.split('.')[0]
    d1 = d[d.STATION.astype(int) == int(stn) ]
    
    d2 = d1.drop(['STATION', 'LATITUDE', 'LONGITUDE', 'ELEVATION', 'NAME']+[c for c in d1.columns if 'ATTRI' in c], axis=1)
    
    miss = {'TEMP':9999.9,  # temperature
            'DEWP':9999.9,  # dew-point temperature
            'SLP':9999.9,   #
            'STP':9999.9, 
            'VISIB':999.9, 
            'WDSP':999.9, 
            'MXSPD':999.9 , 
            'GUST':999.9 , 
            'MAX':9999.9,
            'MIN':9999.9, 
            'PRCP':99.99, 
            'SNDP':999.9, 
            'FRSHTT':99999}
   
    for c in d2.columns: 
        if c == 'FRSHTT': continue
        d2.loc[:,c] = pd.to_numeric(d2.loc[:,c] ,errors='coerce') 
        d2.loc[:,c] = np.where(d2.loc[:,c].values== miss[c], np.NaN, d2.loc[:,c].values)
 

    # F to C
    for c in ['TEMP', 'DEWP', 'MAX', 'MIN']:
        d2.loc[:,c] = (d2.loc[:,c] - 32) * 5./9
    # miles to km
    d2.loc[:,'VISIB'] = d2.loc[:,'VISIB'] * 1.60934
    # inches to cms
    for c in ['PRCP', 'SNDP']: d2.loc[:,c] = d2.loc[:,c]  * 2.54
    # knots to m/s
    for c in ['WDSP', 'MXSPD', 'GUST']: d2.loc[:,c] = d2.loc[:,c]  * 0.514444
    if ofile is None:
        d2.to_csv(odir+os.path.basename(fil), float_format='%.2f')
    else:
        d2.to_csv(odir+ofile, float_format='%.2f')
        
    if rm: os.system('rm '+fil)    
    print(d2)
    return d2
#==============================================================================
#==============================================================================
#==============================================================================



    

#==============================================================================
#==============================================================================
#==============================================================================
def get_stations_info(meta_info, list_isd_history, ofile = 'list_of_stn.csv'):
    # Preprocessing
    la1, lo1, la2, lo2 = meta_info['lat_low'], meta_info['lon_low'], meta_info['lat_up'], meta_info['lon_up']
    sdate, edate = meta_info['st_date'], meta_info['en_date']
    if pd.isnull(sdate): sdate = pd.to_datetime('2100-01-01')
    if pd.isnull(edate): sdate = pd.to_datetime('1900-01-01')
    
    #'station_info/list-isd-history.csv'
    df = pd.read_csv(list_isd_history, index_col=0, parse_dates=[10,11])   
    cond = (df.LAT >= la1)&(df.LAT <= la2)&(df.LON >= lo1)&(df.LON <= lo2) 
    cond1 = cond
    ds = df[cond]
    ds = ds[ (ds.END > edate) & (ds.BEGIN < sdate) ]
    ds['f'] = ds.USAF+ds.WBAN.astype(str)    
    odir = os.path.dirname(ofile)
    print(odir)
    if not os.path.exists(odir): os.makedirs(odir)
    ds.to_csv(ofile)
    
    return ds
#==============================================================================
#==============================================================================
#==============================================================================




#==============================================================================
#==============================================================================
#==============================================================================
def download_daily(station_list_file, output_dir):
    import wget
    url0 = 'https://www.ncei.noaa.gov/data/global-summary-of-the-day/access/'
    ds = pd.read_csv(station_list_file, index_col=0)
    
    # Download from here
    for i1, r1 in ds[:].iterrows():
        
        print(i1, r1)        
        st, en = pd.to_datetime(r1['BEGIN']), pd.to_datetime(r1['END'])
        
        # important
        stn = str(r1.USAF) + str(r1.WBAN )  
        
        if len(stn) != 11: 
            print('error in stn index')
            sys.exit()
        # --------
            
            
        # donwload and process data
        for y in range(st.year, en.year+1)[:]:
            print(y)
                        
            fil = stn+'.csv'
            
            url = url0 + '/'+str(y)+'/'+fil
            # download to directory
            ddir = 'tmp/daily/'+str(y)+'/'
            if not os.path.exists(ddir): os.makedirs(ddir)
            #os.system('wget -O '+ ddir + fil + ' '+url)
            
            #url = ddir + fil + ' '+url
            try:
                
                filename = wget.download(url)
            
                odir = output_dir + '/'+stn+'/'
                if not os.path.exists(odir): os.makedirs(odir)
                    
                if os.path.exists(fil): 
                    ofile = str(y)+'.csv'
                    d2 = process_data(fil, odir, ofile, rm=True)
                    
                else:
                    print(y, 'no exisit')        
            except:
                print(y, 'no exisit')   
    
        
    return
#==============================================================================
#==============================================================================
#==============================================================================




#========================
# check data availability
def check_year_available(output_dir, station_list_file, ofile = 'list_of_data_available_yy.csv'):
    # check download results
    ds = pd.read_csv(station_list_file, index_col=0)
    #print(sorted(glob.glob(od0+'/*/')))
    xx = {}
    for i1, r1 in ds[:].iterrows():
        stn = str(r1['f'])
        odir = output_dir + '/'+ stn +'/'
        ifiles = sorted(glob.glob(odir+'*.csv'))
        print(len(ifiles))
        x = []
        for f in ifiles[:]:
            #print(f)
            d = pd.read_csv(f, index_col=0, parse_dates=True) 
            year = int(f.split('/')[-1].split('.')[0])
            doy = pd.Timestamp(year, 12, 31).dayofyear
            count = d.count()/doy
            count.name = year
            x.append(count)
            #print(count)
            
        if len(x) > 0: 
            ca = pd.concat(x, axis=1)
            print(len(ca))
            if not len(ca)==13: 
                print('error')
            xx[stn] = ca
        #d = pd.concat([ pd.read_csv(f, index_col=0, parse_dates=True) for f in ifiles[:]])
    yy = pd.concat(xx)
    odir = os.path.dirname(station_list_file)
    yy.to_csv(odir+'/'+ofile)
    return yy
#==============================================================================









if __name__ == "__main__":
    #=======
    #ds = xr.open_dataset('/Users/doan/working/share_DATA/DBSH/geo/lu_duong/geo_em.d03.nc')
    #la1, la2, lo1, lo2 = ds.XLAT_M.values.min(), ds.XLAT_M.values.max(), ds.XLONG_M.values.min(), ds.XLONG_M.values.max()
    
    '''
    xx = {'DBSH': {'lat_low': la1,
                       'lat_up': la2,
                       'lon_low': lo1,
                       'lon_up': lo2,
                       'st_date': pd.to_datetime('2023-01-01 00:00:00') , 
                       'en_date': pd.to_datetime('1900-01-01 00:00:00') 
                       }  }
    
    '''
    
    xx = {'Hanoi': {'lat_low': 20,
                   'lat_up': 22,
                   'lon_low': 104.5,
                   'lon_up': 106.5,
                   'st_date': pd.to_datetime('2023-01-01 00:00:00') , 
                   'en_date': pd.to_datetime('1900-01-01 00:00:00') 
                   }}
    
    kk = list(xx.keys())
    for code in kk:
        
        print(code)
        
        
        meta_info = xx[code]
        output_dir = 'ddata/'
        station_list_file = output_dir + '/list_of_stn.csv'
        
        if 1:
            
            list_isd_history = 'station_info/list-isd-history.csv'
            ds = get_stations_info(meta_info, list_isd_history, ofile = station_list_file)
        
        # to plot
        if 0:
            # plot
            fig, ax = plot_basemap(hr=False)
            bu = 1.5
            extent = meta_info['lon_low'], meta_info['lon_up'], meta_info['lat_low'], meta_info['lat_up']
            ax.set_extent(extent)
            #plt.scatter(slon, slat, color='k',s=20, alpha=1)
            #plt.plot([-180,180], [la, la], lw=0.5, c='k')    
            #plt.plot([lo, lo], [-180,180], lw=0.5, c='k')    
            #circle1 = plt.Circle((lo, la), .9, color='r', lw=2, alpha=0.1)
            #ax.add_patch(circle1)
    
            plt.scatter(ds.LON, ds.LAT, color='r',s=20, alpha=1, marker='s')
    
        
        # Download:        
        
        
        if 0: download_daily(station_list_file, output_dir)
        
        
    
    
            
        check_year_available(output_dir, station_list_file)
            
        # 
        if 1:
            f_ay = output_dir + '/' + 'list_of_data_availability.csv'
            da = pd.read_csv(f_ay, index_col=[0,1])
            var = da.index.get_level_values(1).unique()
            
            for v in var[:]:
                
                if not v == 'MAX': continue
                print(v)
                x = da[da.index.get_level_values(1) ==  v ]
                x1 = x.where(x>0.75)
                odir = output_dir+'/data_availability/'
                if not os.path.exists(odir): os.makedirs(odir)
                x1.to_csv(odir+'/'+v+'.csv')
                
                ds = pd.read_csv(output_dir+'/list_of_stn.csv', index_col=0)
                y = x1[~x1.isnull().all(axis=1)].index.get_level_values(0)
                ds1 = ds.set_index('f').loc[y]            
                
                if 1:
                    
                    
                    lo, la = ds.LON.values, ds.LAT.values
                    lo1 = lo.min() - ( lo.max() - lo.min() )*.05
                    lo2 = lo.max() + ( lo.max() - lo.min() )*.05
                    la1 = la.min() - ( la.max() - la.min() )*.05
                    la2 = la.max() + ( la.max() - la.min() )*.05
                    # plot
                    sys.path.append('../../')
                    from common import *
                    fig, ax = plot_basemap(hr=False)
                    bu = 1.5
                    extent = lo1, lo2, la1, la2
                    ax.set_extent(extent)
                    #plt.scatter(slon, slat, color='k',s=20, alpha=1)
                    #plt.plot([-180,180], [la, la], lw=0.5, c='k')    
                    #plt.plot([lo, lo], [-180,180], lw=0.5, c='k')    
                    #circle1 = plt.Circle((lo, la), .9, color='r', lw=2, alpha=0.1)
                    #ax.add_patch(circle1)
            
                    plt.scatter(ds.LON, ds.LAT, color='r',s=20, alpha=1, marker='s')
                    #plt.scatter(lo, la, color='r',s=40, alpha=1, marker='o')   
                    
    
                    plt.scatter(ds1.LON, ds1.LAT, color='b',s=40, alpha=1, marker='o')
                    
                if True:
                    
                    for i1, r1 in ds1[:].iterrows():
                        
                        sname = r1['STATION NAME'].strip()
                        #if not sname == 'SENO': continue   
                        print(sname)                 
                        stn = str(i1)
                        inul = x1.loc[i1].T.isnull()
                        ind = np.argwhere(inul.values[:,0])[-2,0]+1
                        yy = inul.index[ind:-1]
                        if len(yy) == 0: continue
                        #yy = x2.index[~x2.values[:,0]]
                        #yy = x2.index[x2.values[:,0]]
                        #yy = x2.loc[yy[-2]:x2.index[-2]]
                        odir = output_dir + '/'+ stn +'/'
                        #ifiles = sorted(glob.glob(odir+'*.csv'))
                        ifiles = [ odir+y+'.csv' for y in yy if os.path.exists(odir+y+'.csv' )]
                        
                        
    
                        fig = plt.figure(figsize=[6,2.5])
                        ax = plt.axes([.1,.1,.8,.8])                   
                        ax.text(.0,1.05,sname, transform = ax.transAxes)
                    
                        x = []
                        
                        for f in ifiles[:]:
                            #print(f)
                            d = pd.read_csv(f, index_col=0, parse_dates=True) 
                            #d[v].plot(color='gray', lw=.5)
                            d[v].groupby(d.index.month).mean().plot(color='gray', lw=.5)
                        plt.show()
                    
                    
                
                
                
                
            
            
            
            
            
            
            
    
                
                
            

        
        






    
    
    
  
      
                

        
        
        
        
    sys.exit()
        
        
    
    
    
    # Download hourly data
    if False:              
        #!/usr/bin/env python3
        # -*- coding: utf-8 -*-
        """
        Created on Mon Jun  7 16:23:41 2021
        
        @author: doan
        """
        
        
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
            
            df = pd.read_csv('list-isd-history.csv', index_col = 0)[:]
            df.loc[:,'CTRY'] = [c.strip() for c in df.loc[:,'CTRY']]
            b = df[df.CTRY == ctry]
            lat, lon = 42.6977, 23.3219
            lat, lon = 6.5244, 3.3792
            b = df[ (df.LAT < lat+1) & (df.LAT > lat-1)  & 
                    (df.LON < lon+1) & (df.LON > lon-1)  
                    ]
            b['f'] = b.USAF+b.WBAN.astype(str)+'.csv'
            #bo = b[ (b.BEGIN < '2010-01-01') & (b.END > '2015-05-31')]
            url0 = 'https://www.ncei.noaa.gov/data/global-hourly/access'
            
            #sys.exit()
            for i, r in b[:].iterrows():
                print(i, r)
                st, en = pd.to_datetime(r['BEGIN']), pd.to_datetime(r['END'])
                
                for y in range(st.year, en.year+1):
                    print(y)
                    url = url0 + '/'+str(y)+'/'+r['f']
                    
                    odir = ctry+'/'+str(y)
                    if not os.path.exists(odir): os.makedirs(odir)
                    try:
                        filename = wget.download(url, out=odir)
                    except:
                        print('no exisit')
        
        
        
    # download sounding data
        
    if False:
        #!/usr/bin/env python3
        # -*- coding: utf-8 -*-
        """
        Created on Mon Jun  7 16:23:41 2021
        
        @author: doan
        """
        
        
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
                
        
        
    
        
    
        
        
        
        
        
        
    
    