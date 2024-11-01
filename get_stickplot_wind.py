# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 16:41:47 2021

@author: James.manning

extracted from Liu et al 2019 Fig 8 code which reads in a csv file downloaded from NERACOOS "graphing and download"site and generates a stick plot
but this version gets it data from an ERDDAP site
"""
### Import modules
import numpy as np
import netCDF4
import os
from datetime import datetime as dt
from datetime import timedelta
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.dates import date2num,DateFormatter,WeekdayLocator, MONDAY
import math

### HARDCODES #####
datasource='ncep'
if datasource=='neracoos':
    buoy=r'A01'#name of NERACOOS buoy
    timespan=r'mar2020'# assumes you have renamed "download.csv" to "wind_<timespan>_<buoy>"
    folder=os.path.join('C:'+os.sep, 'Users', 'Joann','Downloads'+os.sep)# provides back or forward slashes dependent on OS
elif datasource=='ncep':    
    lat=40.5
    lon=-67.5# degrees east
    starttime=dt(2023,9,16,0,0,0)
    endtime=dt(2023,9,17,0,0,0)
how_to_average='D' #R for raw, W for weekly, D for daily
#hta='Daily' # written out for title in case of raw depending on data source

### Functions
def get_wind_ncep(starttime,endtime,lat,lon):
        #function get a time series of u & v wind m/s from ncep 
        #expects times to be datetimes and lon to be degrees east (~290)
        st=starttime.strftime('%Y-%m-%dT%H:%M:%SZ')
        et=endtime.strftime('%Y-%m-%dT%H:%M:%SZ')
        url='https://coastwatch.pfeg.noaa.gov/erddap/griddap/esrlNcepRe.csvp?uwnd%5B('+st+'):1:('+et+')%5D%5B('+str(lat)+'):1:('+str(lat)+')%5D%5B('+str(lon)+'):1:('+str(lon)+')%5D,vwnd%5B('+st+'):1:('+et+')%5D%5B('+str(lat)+'):1:('+str(lat)+')%5D%5B('+str(lon)+'):1:('+str(lon)+')%5D'
        df=pd.read_csv(url)
        
        t=pd.to_datetime(df['time (UTC)'],format='%Y-%m-%dT%H:%M:%S').values
        u_wind=df['uwnd (m/s)'].values
        v_wind=df['vwnd (m/s)'].values
        #LAT=df['latitude (degrees_north)'][:]
        #LON=df['longitude (degress_east)'][:]
        #for i in range(len(LON)):#transfer lon from (0,360) to (-180.180)
        #    if(LON[i]>180):
        #        LON[i]=-360+LON[i]
        return t,u_wind,v_wind
def sd2uv(s,d): # converts speed and direction to eastward and northward
    u = float(s)*math.sin(math.pi*float(d)/180.)
    v = float(s)*math.cos(math.pi*float(d)/180.)
    return u,v
def stick_plot(time, u, v, **kw):
    width = kw.pop('width', 0.004) 
    headwidth = kw.pop('headwidth', 4)
    headlength = kw.pop('headlength', 6)
    headaxislength = kw.pop('headaxislength', 6)
    angles = kw.pop('angles', 'uv')
    ax = kw.pop('ax', None)
    if angles != 'uv':
        raise AssertionError("Stickplot angles must be 'uv' so that"
                             "if *U*==*V* the angle of the arrow on"
                             "the plot is 45 degrees CCW from the *x*-axis.")
    time, u, v = map(np.asanyarray, (time, u, v))
    if not ax:
        ax = ax1
    q=ax.quiver(date2num(time), [[0]*len(time)], u, v,color='red',scale=20.,width=width, headwidth=headwidth,headlength=headlength, headaxislength=headaxislength) #This Warning says the variable is assigned to but never used but this is needed to plot the wind stress direction  
    ref=10
    ax.quiverkey(q, 0.3, 0.85, ref,
                  "%s m/s" % ref,
                  labelpos='N', coordinates='axes',fontproperties={"weight": "bold","size":40})
    ax.axes.get_yaxis().set_visible(False)
    plt.xticks()
    ax.tick_params(axis="x", labelsize=40)

### read and parse wind data ###
if datasource=='neracoos':
    infile=folder+r"wind_"+timespan+"_"+buoy+".csv"# wind spd & dir as downloaded from NERACOOS and renamed accordingly
    dateparse = lambda x: pd.datetime.strptime(x[0:-4], '%Y-%m-%d %H:%M:%S')#removes the last 4 characters ' UTC'
    df = pd.read_csv(infile, index_col='Time-UTC',parse_dates=['Time-UTC'], date_parser=dateparse)
    time=list(df.index)
    year=df.index[0].year # used later when labeling xaxis
    wspd=df['A01-Hourly-Wind_Speed_m/s'].values
    wdir=df['A01-Hourly-Wind_Direction_degrees'].values
    #convert spd & dir to eastward and northward components of the wind
    uw,vw=[],[]
    for k in range(len(df)):
        [uw1,vw1]=sd2uv(wspd[k],wdir[k])
        uw.append(-1*uw1)# converts to direction towards multiplying by -1
        vw.append(-1*vw1)
    df['uw']=uw
    df['vw']=vw
elif datasource=='ncep':
    print('extracting NCEP wind time series ...')
    [time,uw,vw]=get_wind_ncep(starttime,endtime,lat,360+lon)
    print('done.')
#average wind to daily or weekly    
if how_to_average!='R':
    df=pd.DataFrame({'time':time,'uw':uw,'vw':vw})
    df=df.set_index('time')
    time=df.uw.resample(how_to_average,loffset=-1*timedelta(days=15)).mean().index
    uw=df.uw.resample(how_to_average,loffset=-1*timedelta(days=15)).mean().values
    vw=df.vw.resample(how_to_average,loffset=-1*timedelta(days=15)).mean().values

# START FIGURE, CALL STICK_PLOT, and ANNOTATE PLOT
width=0.4 # not sure if this is needed
fig=plt.figure(figsize=(15,12))
ax1=fig.add_subplot(111)
if how_to_average=='D':
    hta='Daily'
elif how_to_average=='W':
    hta='Weekly'
elif how_to_average=='M':
    hta='Monthly'
if datasource=='neracoos':
    ax1.set_title('Buoy '+buoy+' '+hta+' Wind Stickplot',fontsize=30)     
elif datasource=='ncep':
    ax1.set_title(str(lat)+'N, '+str(lon)+'W '+hta+' Averaged NCEP Wind ',fontsize=30)    
stick_plot(time,uw,vw,color='red')
#mondays = WeekdayLocator(MONDAY)
#ax1.xaxis.set_major_locator(mondays)
#weekFormatter = DateFormatter('%b %d %Y')
fig.autofmt_xdate()
ax1.xaxis_date()
ax1.tick_params(axis="x", labelsize=20)
#ax1.xaxis.set_major_formatter(weekFormatter)
#ax1.text(starttime+timedelta(days=1),0,'yearday 123',rotation=90)
ax1.text(endtime-timedelta(days=19),0,'yearday 195',rotation=90)
ax1.set_xlim(time[0]-timedelta(days=2),time[-1]+timedelta(days=20))
plt.show()
if datasource=='neracoos':
    plt.savefig("wind_"+timespan+"_"+buoy+"_"+hta+".png")
elif datasource=='ncep':
    datet=st=starttime.strftime('%Y_%m_%d')
    plt.savefig("wind_"+str(lat)+'N_'+str(abs(lon))+'W_'+hta+"_starts_"+datet+".png")