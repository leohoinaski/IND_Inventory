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
from verticalProfile import ptVerticalProfile

#%% Gridding and populatingGrid functions

def gridding(lon,lat):
    print('Calling gridding function')
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
    print('Calling populatingGrid function')
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

def griddingMCIP(GRIDCRO2D,GRIDDOT2D,METCRO3D):
    print('Calling griddingMCIP function for MCIP grid')
    #xv, yv = np.meshgrid(lon, lat)
    yv = GRIDDOT2D['LATD'][0,0,:,:]
    xv = GRIDDOT2D['LOND'][0,0,:,:]
    yY = GRIDCRO2D['LAT'][0,0,:,:]
    xX = GRIDCRO2D['LON'][0,0,:,:]

    return xv,yv,xX,yY

def populatingGridTemp3D(dataEmiss,center,emisPar,xv,yv,xX,yY,METCRO3D,METCRO2D,GRIDCRO2D):   
    print('Calling populatingGridTemp3D function for MCIP grid')
    # METCRO3D = nc.Dataset(METCRO3Dfolder)
    # GRIDCRO2D = nc.Dataset(GRIDCRO2Dfolder)
    dataTempo = np.zeros([dataEmiss.shape[1],METCRO3D['TFLAG'][:].shape[0],METCRO3D.NLAYS,
                          METCRO3D.NROWS, METCRO3D.NCOLS],dtype=numpy.single)
    xcenter = center.geometry.centroid.x
    ycenter = center.geometry.centroid.y
    P = METCRO3D['PRES'][:]
    T = METCRO3D['TA'][:]
    TEMP2 = METCRO2D['TEMP2'][:]
    HT = GRIDCRO2D['HT'][:]
    Uas = METCRO2D['WSPD10'][:]
    for ii in range(0,dataEmiss.shape[0]):
        dist = ((xcenter[ii]-xX)**2 + (ycenter[ii]-yY)**2)**(1/2)
        mindist = np.where(dist == np.amin(dist))
        print('------------Source number = ' + str(ii))
        for jj in range(0,dataTempo.shape[1]):
            print('TSTEP = ' + str(jj))
            Pu =  P[jj,:,mindist[0][0],mindist[1][0]]
            Tu = T[jj,:,mindist[0][0],mindist[1][0]]
            TEMP2u = TEMP2[:][jj,0,mindist[0][0],mindist[1][0]]
            Tu = np.append(TEMP2u,Tu)
            HTu = HT[0,0,mindist[0][0],mindist[1][0]]  
            Uasu = Uas[0,0,mindist[0][0],mindist[1][0]] 
            factor = ptVerticalProfile(Pu,Tu,Uasu,HTu,
                                       emisPar.Ts[ii],emisPar.Vs[ii],
                                       emisPar.Ds[ii],emisPar.Hs[ii])
            
            # for zl in range(0,dataTempo.shape[2]):
            #     print('Layer = ' + str(zl))
            for kk in range (0,dataEmiss.shape[1]):
                dataTempo[kk,jj,:,mindist[0][0],mindist[1][0]]= \
                    np.nansum([dataTempo[kk,jj,:,mindist[0][0],mindist[1][0]],
                               dataEmiss.iloc[ii,kk]*factor],0)        
    return dataTempo