#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 07:26:56 2023

@author: JiM
looking at the Cape Cod Bay mooring
"""
import pandas as pd
from netCDF4 import Dataset
import datetime
from matplotlib import pylab as plt
import numpy as np
from conversions import c2f

#get mooring data
df=Dataset('https://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/221p1_rt.nc?waveTime[0:1:25041],waveHs[0:1:25041],sstTime[0:1:149326],sstSeaSurfaceTemperature[0:1:149326]')
time = df.variables['sstTime'][:].data
dtime=[]
for k in range(len(time)):
    dtime.append(datetime.datetime.fromtimestamp(time[k]))
sst  = df.variables['sstSeaSurfaceTemperature'][:].data

#get miniboat data
df=pd.read_csv('https://educationalpassages.org/wp-content/uploads/csv/sensor/Rock_Star_4s.csv')
x5=df['longitude'].values
y5=df['latitude'].values
t5=c2f(df['water_temp'].values)
dt=pd.to_datetime(df['moment_date'])

fig=plt.figure(figsize=(10,8))
ax=fig.add_subplot(111)
ax.plot(dt,t5,color='r',linewidth=3)
ax.plot(dtime,sst,color='b',linewidth=3)
ax.set_ylabel('fahrenheit',fontsize=16)
ax.set_ylim(np.nanmin(t5),np.nanmax(t5))
ax.set_xlim(np.nanmin(dt),np.nanmax(dt))