#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 06:39:08 2023

@author: JiM
"""
import pandas as pd
from matplotlib import pylab as plt

df1=pd.read_csv('/home/user/drift/data/EPMBD_084_1.csv')
x=df1['LON']
y=df1['LAT']
t=df1['buoy_temp_F']

fig = plt.figure()
plot = plt.scatter(x, y, s= 10, c = t, cmap='coolwarm')
fig.colorbar(plot)
plt.grid(True, 'both')

# add another scatterplot
'''
x_line = np.linspace(np.min(x), np.max(x), num=1000)
y_line = x_line + np.sin(np.pi * x_line)
z_line = 5 * x_line
plt.scatter(x_line, y_line, c=z_line, s=0.1, cmap='coolwarm')
'''

plt.show()