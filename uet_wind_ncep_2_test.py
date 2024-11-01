#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 15:10:31 2024

@author: user
"""
st=starttime.strftime('%Y-%m-%dT%H:%M:%SZ')
et=endtime.strftime('%Y-%m-%dT%H:%M:%SZ')

if starttime<datetime.now()-timedelta(days=4):# get reanalysis
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
    df=pd.read_csv('ncep_global_ca54_a976_3c39.csv')# download this once to save time
    #df=pd.read_csv('https://coastwatch.pfeg.noaa.gov/erddap/griddap/NCEP_Global_Best.csvp?ugrd10m%5B('+st+'):1:('+et+')%5D%5B(42.0):1:(42.0)%5D%5B(289.0):1:(289.0)%5D,vgrd10m%5B('+st+'):1:('+et+')%5D%5B(42.0):1:(42.0)%5D%5B(289.0):1:(289.0)%5D')
    #df=pd.read_csv('https://coastwatch.pfeg.noaa.gov/erddap/griddap/NCEP_Global_Best.csvp?ugrd10m%5B('+st+'):1:('+et+')%5D%5B('+str(lat)+'):1:('+str(lat)+')%5D%5B('+str(lon)+'):1:('+str(lon)+')%5D,vgrd10m%5B('+st+'):1:('+et+')%5D%5B('+str(lat)+'):1:('+str(lat)+')%5D%5B('+str(lon)+'):1:('+str(lon)+')%5D')
    t=pd.to_datetime(df['time (UTC)'],format='%Y-%m-%dT%H:%M:%S').values
    #u_wind=mps2knots(df['ugrd10m (m s-1)'].values)
    #v_wind=mps2knots(df['vgrd10m (m s-1)'].values)
    u_wind=df['ugrd10m (m s-1)'].values
    v_wind=df['vgrd10m (m s-1)'].values
    df_w=pd.DataFrame([t,u_wind,v_wind])
print('got wind!')