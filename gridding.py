#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 11:15:14 2021

@author: leohoinaski
"""

import geopandas as gpd
import numpy as np
from shapely.geometry import MultiLineString
from shapely.ops import polygonize
#from temporalDisagregation import temporalDisagVehicular
import numpy.matlib

#%% Gridding and populatingGrid functions

def gridding(lon,lat):
    xv, yv = np.meshgrid(lon, lat)
    hlines = [((x1, yi), (x2, yi)) for x1, x2 in zip(lon[:-1], lon[1:]) for yi in lat]
    vlines = [((xi, y1), (xi, y2)) for y1, y2 in zip(lat[:-1], lat[1:]) for xi in lon]
    grids = list(polygonize(MultiLineString(hlines + vlines)))
    grids = gpd.GeoDataFrame(grids) 
    grids.columns =['geometry'] 
    grids['geometry'] = grids['geometry']
    grids.crs = "EPSG:4326"  
    grids['X'] = grids.geometry.centroid.x
    grids['Y'] = grids.geometry.centroid.y
    xX = np.array(grids['X']).reshape((lon.shape[0]-1,lat.shape[0]-1)).transpose()
    yY = np.array(grids['Y']).reshape((lon.shape[0]-1,lat.shape[0]-1)).transpose()
    return grids,xv,yv,xX,yY

def populatingGrid(dataEmiss,center,xX,yY,xv,yv):   
    data = np.zeros([1,dataEmiss.shape[1],np.size(yv,0)-1, np.size(xv,1)-1])
    xcenter = center.geometry.centroid.x
    ycenter = center.geometry.centroid.y
   
    for ii in range(0,dataEmiss.shape[0]):
        dist = ((xcenter[ii]-xX)**2 + (ycenter[ii]-yY)**2)**(1/2)
        mindist = np.where(dist == np.amin(dist))
        print('Fire number '+str(ii)+' from '+str(dataEmiss.shape[0]))
        for kk in range (0,dataEmiss.shape[1]):
            data[0,kk,mindist[0][0],mindist[1][0]]= data[0,kk,mindist[0][0],mindist[1][0]]+dataEmiss.iloc[ii,kk]     
    return data

def populatingGridMat(dataMat,center,xv,yv,xX,yY):   
    dataTempo = np.zeros([dataMat.shape[2],dataMat.shape[1],np.size(yv,0)-1, np.size(xv,1)-1])
    xcenter = center.geometry.centroid.x
    ycenter = center.geometry.centroid.y
   
    for ii in range(0,dataMat.shape[0]):
        dist = ((xcenter[ii]-xX)**2 + (ycenter[ii]-yY)**2)**(1/2)
        mindist = np.where(dist == np.amin(dist))
        print('cell number = ' + str(ii))
        for kk in range (0,dataMat.shape[1]):
            dataTempo[:,kk,mindist[0][0],mindist[1][0]]= np.nansum([dataTempo[:,kk,mindist[0][0],mindist[1][0]],dataMat[ii,kk,:]],0)        
    return dataTempo