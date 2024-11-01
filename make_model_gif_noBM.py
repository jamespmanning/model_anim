#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Apr 17 11:29:14 2019
@author: leizhao
Modified by JiM in Fall 2019
Modified by JiM in Summer 2020 with Jack Polentes help
Modified by JiM in Fall 2020 by adding drifter tracks and, in late Oct 2020, NCEP wind vector
Modified by JiM in Spring 2021 to add all model options including FVCOM
Modified by JiM in Fall 2021 to add OOI variables 
Modified by JiM in Winter 2023-24 to add emolt obs and ptown wind, and rid of basemap
Stored in Github as "model_anim" repository
See hardcodes at top of code where there is, for example, "start_date" and "ndays" 
Note: You may need to adjust the "clevs" to get good range of colors.
Note: You may want to adjust the Basemap resolution to "c" for crude in experiments.
Note: You may need to "conda install -c conda-forge basemap-data-hires" in order to use the higher resolution coastlines
Note: 
"""
#hardcodes########################
#area='NEC'#'SNW'#'NorthShore'#'SNW'#'GBANK_RING'#'Gloucester'
area='inside_CCBAY'
start_date='2023-12-22'#'2013-04-01'
#clevs=[65.,80.,.5]#gb ring surf June
clevs=[50.,72.,.5]#ns bottom June
clevs=[52.,78.,.5]#gb ring June
clevs=[58.,74.,.5]#SNE-W in July
clevs=[58.,80.,.5]#SNE in August
#clevs=[56.,75.,.5]#NEC in December
clevs=[43.,52.,.2]#CCBay in Dec
clevs=[34.,48.,.2]#CCBay in late Dec
dtime=[]
units='degF'
ndays=25 # number of days to run
detide='n'# for FVCOM hindast only now
include_temp_obs='yes' # 'yes' overlays positions of observed bottom temps
include_temp_obs_LFA='no' # 'yes' overlays LFA positions of observed bottom temp
include_wind='yes'
ptownwind_file='CCS_12-1-23_12-00_AM_3_Month_1705326014_v2.csv'
include_ooi='no'
include_weather_balloon='no'
include_miniboats='no'
#cluster='ep_2022_1' # code of drifter cluster to include
cluster='fhs_2023_2'
surf_or_bot=0#0 for most put -1 for surface (opposite for FVCOM)
lat_w,lon_w=41.9,-70.3 # base of wind vector legend (actual vector appears mid plot)
model='FVCOM'#'DOPPIO'# FVCOM, DOPPIO, or GOMOFS ...
subsample=5
maxvel = 1.0 # for quiver legend
#########
import os,imageio
import conda
import pandas as pd
from pylab import *
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
from datetime import datetime,timedelta,timezone
import time
import zlconversions as zl
from conversions import mps2knots
import sys
import warnings
warnings.filterwarnings("ignore") # gets rid of warnings at runtime but you may want to comment this out to see things
import netCDF4 # for wind access
from math import sqrt
try:
    import cPickle as pickle
except ImportError:
    import pickle
import glob
import ssl # added this in Oct 2021 when I got a "certificate_verify_failed" message
ssl._create_default_https_context = ssl._create_unverified_context

def getgbox(area):
  # gets geographic box based on area
  if area=='SNE':
    gbox=[-71.,-67.,39.5,42.] # for SNE shelf east
  elif area=='SNW':
    gbox=[-71.5,-69.5,40.,41.75] # for SNw shelf west
  elif area=='MABN':
    gbox=[-73.,-68.,39.,42.] # for SNw shelf west  
  elif area=='OOI':
    gbox=[-71.5,-69.,39.5,41.6] # for OOI
  elif area=='GBANK':
    gbox=[-71.,-66.,40.,42.] # for GBANK
  elif area=='GBANK_RING':
    gbox=[-71.,-65.,39.,42.5] # for typical GBANK Ring 
  elif area=='GS':           
    gbox=[-71.,-63.,38.,42.5] # for Gulf Stream
  elif area=='NorthShore':
    gbox=[-71.,-69.5,41.75,43.25] # for north shore
  elif area=='Gloucester':
    gbox=[-71.,-70.,42.25,43.] # for north shore
  elif area=='IpswichBay':
    gbox=[-71.,-70.,42.5,43.] # for IpswitchBay
  elif area=='CCBAY':
    gbox=[-70.75,-69.8,41.5,42.23] # CCBAY
  elif area=='inside_CCBAY':
    gbox=[-70.75,-70.,41.7,42.15] # inside CCBAY
  elif area=='NEC':
    gbox=[-68.,-63.,38.,43.5] # NE Channel
  elif area=='NE':
    gbox=[-76.,-66.,35.,44.5] # NE Shelf 
  return gbox

def get_wind_ncep_2(starttime,endtime,lat,lon): # added Sep 2023
    #function get a time series of u & v wind m/s from ncep returns knots using "mps2knots" function
    #added option for both reanalysis (historical) and forecast in Dec 2023
    st=starttime.strftime('%Y-%m-%dT%H:%M:%SZ')
    et=endtime.strftime('%Y-%m-%dT%H:%M:%SZ')
    
    if starttime<datetime.now()-timedelta(days=60):# get reanalysis
        #expects times to be datetimes and lon to be degrees east (~290)
        print('extracting wind time series from CoastWatch reanalysis erddap...')
        url='https://coastwatch.pfeg.noaa.gov/erddap/griddap/esrlNcepRe.csvp?uwnd%5B('+st+'):1:('+et+')%5D%5B('+str(lat)+'):1:('+str(lat)+')%5D%5B('+str(360+lon)+'):1:('+str(360+lon)+')%5D,vwnd%5B('+st+'):1:('+et+')%5D%5B('+str(lat)+'):1:('+str(lat)+')%5D%5B('+str(360+lon)+'):1:('+str(360+lon)+')%5D'
        df=pd.read_csv(url)       
        t=pd.to_datetime(df['time (UTC)'],format='%Y-%m-%dT%H:%M:%S').values
        u_wind=df['uwnd (m/s)'].values
        v_wind=df['vwnd (m/s)'].values
        df_w=pd.DataFrame([t,u_wind,v_wind])
    else: # get GFS forecast gets 42N and 71W
        print('extracting wind time series from CoastWatch GFS forecast erddap...')
        lon=360+lon
        #df=pd.read_csv('ncep_global_ca54_a976_3c39.csv')# download this once to save time
        #df=pd.read_csv('ncep_global_48dc_3bd9_1414.csv') # startes 12/23
        #df=pd.read_csv('https://coastwatch.pfeg.noaa.gov/erddap/griddap/NCEP_Global_Best.csvp?ugrd10m%5B('+st+'):1:('+et+')%5D%5B(42.0):1:(42.0)%5D%5B(289.0):1:(289.0)%5D,vgrd10m%5B('+st+'):1:('+et+')%5D%5B(42.0):1:(42.0)%5D%5B(289.0):1:(289.0)%5D')
        df=pd.read_csv('https://coastwatch.pfeg.noaa.gov/erddap/griddap/NCEP_Global_Best.csvp?ugrd10m%5B('+st+'):1:('+et+')%5D%5B('+str(lat)+'):1:('+str(lat)+')%5D%5B('+str(lon)+'):1:('+str(lon)+')%5D,vgrd10m%5B('+st+'):1:('+et+')%5D%5B('+str(lat)+'):1:('+str(lat)+')%5D%5B('+str(lon)+'):1:('+str(lon)+')%5D')
        t=pd.to_datetime(df['time (UTC)'],format='%Y-%m-%dT%H:%M:%S').values
        u_wind,v_wind=[],[]
        for k in range(len(t)):
            u_wind.append(mps2knots(df['ugrd10m (m s-1)'].values[k]))
            v_wind.append(mps2knots(df['vgrd10m (m s-1)'].values[k]))
        #u_wind=df['ugrd10m (m s-1)'].values
        #v_wind=df['vgrd10m (m s-1)'].values
        df_w=pd.DataFrame([t,u_wind,v_wind])
    print('got wind!')    
    return df_w
    
def plotit(lons,lats,slons,slats,stemp,temp,depth,time_str,path_save,dpi=80,area='OOI',clevs=[39.,44.,0.5],lat_w=42.9,lon_w=-70.3,wind=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]]),dtime=dtime,model='DOPPIO',detide='n'):
    '''dtimes,u,v are wind in m/s'''
    if clevs[1]<32.:
        units='degC'
    else:
        units='degF'
    fig = plt.figure(figsize=(10,6))
    #ax = fig.add_axes([0.05,0.05,0.9,0.82])
    ax = subplot(111,aspect=(1.0/cos(42.*pi/180.0)))
    gb=getgbox(area)
    if include_wind=='yes':

       # first plot the NCEP wind
       xw,yw=lon_w,lat_w 
       ax.quiver(xw,yw,10.0,0.0,color='teal',scale=50.,zorder=10)#wind legend in middle
       ax.text(xw,yw,'~20 knot wind',fontsize=12,color='teal',zorder=10,va='bottom')
       wind=wind.T
       wind=wind.set_index(0)
       idex=wind.index.get_loc(datetime.strptime(time_str[0:15],"%Y-%m-%d %H%M"),method='nearest')
       ax.quiver(xw,yw,wind[1][idex],wind[2][idex],color='teal',scale=50.0,zorder=9) # plot an arrow
       
       # here we load the ptown wind as downloaded from weatherlink with particular encoding found using "import chartdet;endoding=chardet.detect(filename.read())  
       dfptw=pd.read_csv(ptownwind_file,skiprows=5, encoding='ISO-8859-1')
       dfptw.drop(dfptw[dfptw['Wind Speed - mph']=='--'].index, inplace=True)
       dfptw.drop(dfptw[dfptw['Wind Direction']=='--'].index, inplace=True)
       dfptw['datet']=pd.to_datetime(dfptw['Date & Time'])
       dfptw['datet'] = dfptw['datet'].dt.tz_localize(timezone.utc)
       plottime=pd.to_datetime(time_str,utc=True)
       idx = abs(plottime - dfptw['datet']).idxmin()
       d={'':0, '--':0,'N':0, 'NNE':22.5,"NE":45,"ENE":67.5, 'E':90,'ESE':112.5, 'SE':135,'SSE':157.5, 'S':180,'SSW':202.5, 'SW':225,'WSW':247.5, 'W':270,'WNW':292.5,'NW':315,'NNW':337.5, 'N':0,'North':0,'East':90,'West':270,'South':180}
       dfptw['wdir']=dfptw['Wind Direction'].str.strip().map(d)
       [up,vp]=zl.sd2uv(dfptw['Wind Speed - mph'][idx],dfptw['wdir'][idx])
       mptwx,mptwy=-70.1864,42.0629
       Q=ax.quiver(mptwx,mptwy,-1*up,-1*vp,scale=50.0,color='teal',zorder=8)
       #maxstr='%3.1f m/s (~2 knots)' % maxvel
       #qk = quiverkey(Q,0.93,0.08,maxvel,maxstr,labelpos='W')
       
    if len(slons)!=0:
        x1,y1=slons,slats
        for jj in range(len(stemp)):
            ax.text(x1[jj],y1[jj],str.format('{0:.1f}',stemp[jj]*1.8+32),color='w',fontsize=10,fontweight='bold',horizontalalignment='center',verticalalignment='center')

    if (area=='IpswichBay') or (area=='Gloucester'):
        labint=0.2
        dept_clevs=[30,50,100, 150]
    elif area[0:3]=='CCB':
        labint=0.5
        dept_clevs=[30,50,100]
    elif (area=='NE') or (area=='SNW'):
        labint=1.0
        dept_clevs=[50,100,1000]
        x,y=-69.5,40.5    
        plt.text(x,y,'Great South Channel',fontsize=12, rotation=0) 
    elif (area=='NorthShore'):
        labint=.50
        dept_clevs=[50,100,150]  
    elif (area=='GBANK') or (area=='GBANK_RING'):
        labint=1.0
        dept_clevs=[50,100,150]
        x,y=-68.3,40.65    
        plt.text(x,y,' Georges Bank',fontsize=16, rotation=30) 
    elif (area=='inside_CCBAY'):
        labint=.2
        dept_clevs=[30,50,100,150]
        x,y=-70.32,41.84    
        ax.plot(x,y,'o',color='k',markersize=10)
        ax.text(x,y,'  mooring',fontsize=12)#, rotation=30)
    else:
        labint=1.0
        dept_clevs=[30,50,100, 150,300,1000]
    x, y = lons, lats # compute map proj coordinates.
    dtthis=datetime(int(time_str[0:4]),int(time_str[5:7]),int(time_str[8:10]),int(time_str[11:13]),int(time_str[13:15]))
    clevs=np.arange(clevs[0],clevs[1],clevs[2])  #for all year:np.arange(34,84,1) or np.arange(34,68,1)
    if (model=='DOPPIO') or (model=='GOMOFS'):
            # draw filled contours.
            cs = m.contourf(x,y,temp,clevs,cmap=plt.get_cmap('rainbow'))
            # draw depth contours.
            dept_cs=m.contour(x,y,depth,dept_clevs,colors='black')
            plt.clabel(dept_cs, inline = True, fontsize =12,fmt="%1.0f")

    elif model=='FVCOM':
        
        if dtthis<datetime(2020,6,1,0,0,0):
            url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
        elif dtthis>datetime.now()-timedelta(days=2):    
            url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
            url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_GOM7_FORECAST.nc'
        else:
            url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_NORTHEAST_FORECAST.nc'
        #print(url)
        nc = netCDF4.Dataset(url).variables
        time_var = nc['time']
        itime = netCDF4.date2index(dtthis,time_var,select='nearest')
        # Get lon,lat coordinates for nodes (depth)
        lats = nc['lat'][:]
        lons = nc['lon'][:]
        # Get lon,lat coordinates for cell centers (depth)
        latc = nc['latc'][:]
        lonc = nc['lonc'][:]
        # Get Connectivity array
        nv = nc['nv'][:].T - 1 
        # Get depth
        depth = nc['h'][:]  # depth
        #print('most nc commands done')
        #figure(figsize=(12,8))
        #subplot(111,aspect=(1.0/cos(mean(lats)*pi/180.0)))

        dtime = netCDF4.num2date(time_var[itime],time_var.units)
        daystr = dtime.strftime('%Y-%b-%d %H:%M')
        
        #x, y = m(lons, lats) # compute map proj coordinates.
        
        if surf_or_bot==3:# vertically averaged
            temp = nc['temp'][itime,0,:]#if len(clevs)==0:
            u = nc['ua'][itime,:]# surface fields
            v = nc['va'][itime,:]
        else: 
            #print('geting temp converted, u, and v ..')
            temp = nc['temp'][itime, surf_or_bot, :]*1.8+32
            u = nc['u'][itime,surf_or_bot,:]
            v = nc['v'][itime,surf_or_bot,:]
            #print('done')
        if len(clevs)==0:
        	clevs=np.arange(int(min(temp)),int(max(temp)+1),1)
        # find velocity points in bounding box
        idv = np.argwhere((lonc >= gb[0]) & (lonc <= gb[1]) & (latc >= gb[2]) & (latc <= gb[3]))
        idv=idv[0::subsample]
        tri = Tri.Triangulation(lons,lats, triangles=nv)
        plt.tricontourf(tri, temp,levels=clevs,shading='faceted',cmap=plt.cm.rainbow,zorder=0)
        axis(gb)
        
        gca().patch.set_facecolor('0.5')
        cbar=colorbar()
        #cbar.set_label('SST (degF)', rotation=-90)
        Q = ax.quiver(lonc[idv],latc[idv],u[idv],v[idv],scale=20)
        maxstr='%3.1f m/s (~2 knots)' % maxvel
        qk = quiverkey(Q,0.93,0.08,maxvel,maxstr,labelpos='W')


    # add colorbar.
    #cbar = colorbar(cs,ax=cs,location='right',pad="2%",size="5%")
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    #cbar.set_label(units,fontsize=25)

    df=pd.read_csv('/home/user/drift/drift_'+cluster+'.csv')
    #df1=pd.read_csv('http://apps-nefsc.fisheries.noaa.gov/drifter/drift_ep_2022_1_ap3.csv')
    #df2=pd.read_csv('http://apps-nefsc.fisheries.noaa.gov/drifter/drift_wms_2022_1.csv')
    #df=pd.concat([df1,df2],axis=0)
    ids=np.unique(df['ID'])
    #print(ids)
    for k in ids:
        #print(k)
        datett=[]
        df1=df[df['ID']==k]
        for kk in range(len(df1)):
            datett.append(datetime(int(start_date[0:4]),1,1)+timedelta(df1['YEARDAY'].values[kk]-1))
        df1['datet']=datett
        df1=df1[df1['datet']<dtthis]# only the track less than this time
        df1=df1[df1['datet']>dtthis-timedelta(2.0)]
        #make a plot for this hour and every hour previous with decreasing thickness
        if k==226400691:
            print(len(df1))
        for kk in range(len(df1)-1):
            x1,y1=df1['LON'].values[kk],df1['LAT'].values[kk]
            x2,y2=df1['LON'].values[kk+1],df1['LAT'].values[kk+1]
            plt.plot([x1,x2],[y1,y2],'m',linewidth=kk/(len(df1)/4))#markersize=kk)#,linewidth=kk)
    
    #if sys.argv[1]=='drift_whs_2022_1': # here is where we are adding miniboat track from Cassie
    #include_miniboats='yes'
    df=pd.read_csv('https://educationalpassages.org/wp-content/uploads/csv/LadyLance_1.csv')
    df1=pd.DataFrame()
    df1['datet']=pd.to_datetime(df['moment_date'])
    df1['LON']=df['longitude']
    df1['LAT']=df['latitude']
    df1=df1[df1['datet']<dtthis]# only the track less than this time
    df1=df1[df1['datet']>dtthis-timedelta(2.0)]
    for kk in range(len(df1)-1):
            x1,y1=df1['LON'].values[kk],df1['LAT'].values[kk]
            x2,y2=df1['LON'].values[kk+1],df1['LAT'].values[kk+1]
            plt.plot([x1,x2],[y1,y2],'w',linewidth=kk/(len(df1)/4))#markersize=kk)#,linewidth=kk)

   
    # add title
    clayer=''# default current layer
    if (surf_or_bot==-1) & (model[1]=='O'): #case of both DOPPIO and GOMOFS
        layer='surface'
    elif (surf_or_bot==0) & (model[1]=='O'):
        layer='bottom'
    elif (surf_or_bot==-1) & (model=='FVCOM'):
        layer='bottom'
    elif (surf_or_bot==0) & (model=='FVCOM'):
        layer='surface'
    elif surf_or_bot==3:
        clayer='VA_'
        layer='surface'
    #plt.title(model+' '+layer+' temps (color),'+clayer+'current (black) & depth (m)',fontsize=12,fontweight='bold')
    #plt.title('eMOLT bottom temps (white#s) '+model+' '+layer+' temps (color) '+clayer+' & depth (meters)',fontsize=12,fontweight='bold')
    #plt.title('eMOLT temps (white#s), '+model+' '+layer+' temps (color) '+clayer+', depth(meters), drifter (purple), miniboat (white)',fontsize=10,fontweight='bold')
    plt.title('drifters (purple), '+model+' '+layer+' model temp (color) & current (black) '+clayer+', NCEP daily wind (teal arrow), depth (m) ',fontsize=10,fontweight='bold')
    if detide=='y':
        time_str=time_str[0:10]
    #plt.suptitle(model+' at '+time_str, fontsize=24) 
    plt.suptitle(time_str, fontsize=24) 
    if not os.path.exists(path_save):
        os.makedirs(path_save)
    plt.savefig(os.path.join(path_save,time_str.replace(' ','t')+'.png'),dpi=dpi)
    plt.close()

def mean_temp(temps): #makes daily averages if needed for multi-week animations
    mean_temp=temps[0,0]
    for i in range(1,24):
        mean_temp+=temps[i,0]# note we are using the bottom level 0
    return mean_temp/24.0
        
def make_images(dpath,path,dt=datetime(2019,5,1,0,0,0),interval=31,area='OOI',clevs=[39.,44.,0.5],lat_w=42.5,lon_w=-70.3,wind=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]]),ooi=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]]),model=model,detide='n',ooiw=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]])):
    '''dpath: the path of dictionary, use to store telemetered data
        path: use to store images
        dt: start time
        interval: how many days we need make 
        dtimes,u,v are wind in m/s
    '''
    with open(dpath,'rb') as fp:
         telemetered_dict=pickle.load(fp)
    if detide=='n':     
        interval=interval*24
        tdint=timedelta(hours=1)
    else:
        interval=interval
        tdint=timedelta(days=1)
        #tdint=timedelta(hours=1)
    for j in range(interval):
        #dtime=dt+timedelta(days=j)
        dtime=dt+tdint*j
        print(dtime)

        if model=='DOPPIO':
            #url=get_doppio_url(dtime)
            url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best'
            while True:
                if zl.isConnected(address=url):
                    break
                print('check if website is well or internet is connected?')
            time.sleep(5)
            skip=0
            #while True: 
            #        try:
            nc = NetCDFFile(url)
            lons=nc.variables['lon_rho'][:]
            lats=nc.variables['lat_rho'][:]
            temps=nc.variables['temp']
            depth=nc.variables['h'][:]
            time_var = nc['time']
            itime = netCDF4.date2index(dtime,time_var,select='nearest')
            m_temp=temps[itime,surf_or_bot]#-1 for surface?
            temp=m_temp*1.8+32
            if skip==1:
                continue
            
        elif model=='GOMOFS':
            url=gomofs_modules.get_gomofs_url(dtime)
            nc = NetCDFFile(url)
            lons=nc.variables['lon_rho'][:]
            lats=nc.variables['lat_rho'][:]
            temps=nc.variables['temp']
            depth=nc.variables['h'][:]
            #time_var = nc['time']
            #itime = netCDF4.date2index(dtime,time_var,select='nearest')
            m_temp=temps[0,surf_or_bot,:,:]#-1 for surface
            temp=m_temp*1.8+32

        elif model=="FVCOM":
            #print('fvcom')
            lons=[]
            lats=[]
            depth=[]
            temp=[]
                
        ntime=dtime
        time_str=ntime.strftime('%Y-%m-%d %H%MUTC')
        
        # this is where we get eMOLT data
        Year=str(ntime.year)
        Month=str(ntime.month)
        Day=str(ntime.day)
        slons,slats=[],[]
        try:
            slons,slats,stemp=[],[],[]
            for i in telemetered_dict[Year][Month][Day].index:
                slons.append(telemetered_dict[Year][Month][Day]['lon'].iloc[i])
                slats.append(telemetered_dict[Year][Month][Day]['lat'].iloc[i])
                stemp.append(telemetered_dict[Year][Month][Day]['temp'].iloc[i])
        except:
            slons,slats,stemp=[],[],[]
        #print(slons,temp)            
        plotit(lons,lats,slons,slats,stemp,temp,depth,time_str,path,dpi=80,area=area,clevs=clevs,lat_w=lat_w,lon_w=lon_w,wind=wind,dtime=dtime,model=model,detide=detide)
            
        
def read_telemetry(path):
    """read the telemetered data and fix a standard format, the return the standard data"""
    tele_df=pd.read_csv(path,sep='\s+',names=['vessel_n','esn','month','day','Hours','minutes','fracyrday',\
                                          'lon','lat','dum1','dum2','depth','rangedepth','timerange','temp','stdtemp','year'])
    if len(tele_df)<6000:
        print('Warning: the emolt.dat file is not complete at this time.')
        #sys.exit()
        
    return tele_df

def make_dict(filepathsave,ptelemetered='http://emolt.org/emoltdata/emolt.dat'):   
    '''create a dictionary use to store the data from telemetered, index series is year, month, day and hour
    ptelemetered: the path of telemetered
    '''
    dfdict={}
    df=read_telemetry(ptelemetered)
    
    # JiM added the following 6/16/2021 to get mean value per vessel per day
    df=df[df['depth']>2.0]
    df = df.groupby(['vessel_n','year','month','day'], as_index=False).agg({'temp': 'mean', 'lat': 'mean', 'lon': 'mean','Hours':'first','minutes':'first'})
    for i in df.index:
        #if df['depth'][i]<2.0:
        #    continue
        #if df['minutes'].iloc[i]<=30:
        Ctime=datetime.strptime(str(df['year'].iloc[i])+'-'+str(df['month'].iloc[i])+'-'+str(df['day'].iloc[i])+' '+\
                                         str(df['Hours'].iloc[i])+':'+str(df['minutes'].iloc[i])+':'+'00','%Y-%m-%d %H:%M:%S')
        #else:
        #    Ctime=datetime.strptime(str(df['year'].iloc[i])+'-'+str(df['month'].iloc[i])+'-'+str(df['day'].iloc[i])+' '+\
        #                                 str(df['Hours'].iloc[i])+':'+str(df['minutes'].iloc[i])+':'+'00','%Y-%m-%d %H:%M:%S')+timedelta(seconds=1800)
        Year=str(Ctime.year)
        Month=str(Ctime.month)
        Day=str(Ctime.day)
        #Vessel=df['vessel_n'][i]
        if not Year in dfdict:
            dfdict[Year]={}
        if not Month in dfdict[Year]:
            dfdict[Year][Month]={}
        if not Day in dfdict[Year][Month]:
            dfdict[Year][Month][Day]={}
        #if not Vessel in dfdict[Year][Month][Day]:
        #    dfdict[Year][Month][Day][Vessel]={}

        if len(dfdict[Year][Month][Day])!=0:
            dfdict[Year][Month][Day]=dfdict[Year][Month][Day].append(pd.DataFrame(data=[[df['lat'].iloc[i],df['lon'].iloc[i],df['temp'].iloc[i]]],columns=['lat','lon','temp']).iloc[0])
            dfdict[Year][Month][Day].index=range(len(dfdict[Year][Month][Day]))
        else:
            dfdict[Year][Month][Day]=pd.DataFrame(data=[[df['lat'].iloc[i],df['lon'].iloc[i],df['temp'].iloc[i]]],columns=['lat','lon','temp'])
    with open(filepathsave,'wb') as fp:
        pickle.dump(dfdict,fp,protocol=pickle.HIGHEST_PROTOCOL)


def make_gif(gif_name,png_dir,start_time=False,end_time=False,frame_length = 0.2,end_pause = 4 ):
    '''use images to make the gif
    frame_length: seconds between frames
    end_pause: seconds to stay on last frame
    the format of start_time and end time is string, for example: %Y-%m-%d(YYYY-MM-DD)'''
    
    if not os.path.exists(os.path.dirname(gif_name)):
        os.makedirs(os.path.dirname(gif_name))
    allfile_list = glob.glob(os.path.join(png_dir,'*.png')) # Get all the pngs in the current directory
    #print(allfile_list)
    file_list=[]
    if start_time:    
        for file in allfile_list:
            if start_time<=os.path.basename(file).split('.')[0]<=end_time:
                file_list.append(file)
    else:
        file_list=allfile_list
    #list.sort(file_list, key=lambda x: x.split('/')[-1].split('t')[0]) # Sort the images by time, this may need to be tweaked for your use case
    file_list.sort()
    print(file_list)
    images=[]
    # loop through files, join them to image array, and write to GIF called 'wind_turbine_dist.gif'
    for ii in range(0,len(file_list)):       
        file_path = os.path.join(png_dir, file_list[ii])
        if ii==len(file_list)-1:
            for jj in range(0,int(end_pause/frame_length)):
                images.append(imageio.imread(file_path))
        else:
            images.append(imageio.imread(file_path))
    # the duration is the time spent on each image (1/duration is frame rate)
    imageio.mimsave(gif_name, images,'GIF',duration=frame_length)
    
#MAINCODE###########################
start_date_datetime=datetime(int(start_date[0:4]),int(start_date[5:7]),int(start_date[8:10]),0,0,0)
end_date_datetime=datetime(int(start_date[0:4]),int(start_date[5:7]),int(start_date[8:10]),0,0,0)+timedelta(days=ndays)
end_date=str(end_date_datetime.year)+'-'+str(end_date_datetime.month).zfill(2)+'-'+str(end_date_datetime.day).zfill(2)
realpath=os.path.dirname(os.path.abspath(__file__))
dpath=realpath[::-1].replace('py'[::-1],'result/Doppio'[::-1],1)[::-1]  # the directory of the result
if not os.path.exists(dpath):
    os.makedirs(dpath)
dictionary=os.path.join(dpath,'dictionary_emolt.p')
gif_path=os.path.join(dpath,'gif')
map_save=os.path.join(dpath,'map')
gif_name =os.path.join(gif_path,start_date+area+'_'+model+'_'+str(surf_or_bot)+'.gif')
if include_wind=='yes': # get wind time series at one location
    #df_w=get_wind_ncep(int(start_date[0:4]),lat_w,lon_w)# returns a dataframe of wind
    df_w=get_wind_ncep_2(start_date_datetime,end_date_datetime,lat_w,lon_w)# returns a dataframe of wind
else:
    df_w=np.nan

#############################
#run functions
if include_temp_obs=='yes':
    make_dict(filepathsave=dictionary)
make_images(dpath=dictionary,path=map_save,dt=start_date_datetime,interval=ndays,area=area,clevs=clevs,lat_w=lat_w,lon_w=lon_w,wind=df_w,model=model,detide=detide)
make_gif(gif_name,map_save,start_time=start_date,end_time=end_date)# sending datetimes
    


