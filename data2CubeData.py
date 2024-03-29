#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import geopandas as gpd
import pandas as pd
import numpy as np
#from netCDFcreator import createNETCDFtemporal
#from temporalDisagregation import temporalDisagVehicular
import numpy.matlib
from shapely.geometry import Polygon
from gridding import griddingMCIP, populatingGridTemp3D
from IND_speciate import IndSpeciate
#import datetime
import netCDF4 as nc


#%% MCIPDomain
def MCIPDomain(GRIDDOT2D):   
    print('Extracting MCIP GRIDDOT2D coordinates')
    
    #dataVar = list(ds.variables.keys())
    x = np.unique(GRIDDOT2D['LOND'][:])
    y = np.unique(GRIDDOT2D['LATD'][:])
    #-------------------- Reading industrial emissions-----------------------------   
    # Creating domain window
    domain = Polygon(zip([np.min(x),np.min(x),np.max(x),np.max(x)],
                         [np.min(y),np.max(y),np.max(y),np.min(y)])) 
    domain = gpd.GeoDataFrame(index=[0],geometry=[domain])
    domain.crs = "EPSG:4326"
    #Reading
    return domain,x,y

#%% emissReader
def emissReader(rootPath,domain):
    convKg2g = 1000 # Conversion from kg/s to g/s
    print('Reading indutrial emissions from ' + rootPath+'/Inputs/BR_Ind.xlsx')
    dfind = pd.read_excel(rootPath+'/Inputs/BR_Ind.xlsx')
    
    # Converting to geodataframe
    ind = gpd.GeoDataFrame(
        dfind, geometry=gpd.points_from_xy(dfind.Long, dfind.Lat))
    ind.crs = "EPSG:4326" # Setting EPSG
    ind = ind[ind['Type']=='POINT'].copy()
    
    # Cliping industries inside domain
    pip_mask = ind.within(domain.iloc[0,0])
    ind = ind.loc[pip_mask]
    ind=ind.reset_index(drop=True)

    # Selecting data from industrial inventory - emissions in kg/s
    # Measured data
    dataEmissIND = ind[['PMemis', 
                        'COemis', 
                        'NOxemis',
                        'SOxemis',
                        'VOCemis']].copy() 
    
    # Converting from kg/s to g/s
    dataEmissIND=dataEmissIND*convKg2g
    
    # Estimated data
    dataEmissINDestmeas = ind[['PMestmeas', 
                        'COestmeas', 
                        'NOxestmeas',
                        'SOxestmeas',
                        'VOCestmeas']].copy()
    
    # Converting from kg/s to g/s
    dataEmissINDestmeas=dataEmissINDestmeas*convKg2g 

    # Substituting nan by estimated data
    isna = dataEmissIND.isna()
    dataEmissIND[isna] = dataEmissINDestmeas[isna]
    dataEmissIND['ID'] = ind['ID']
    
    # Getting emissions centroid
    centerIND = ind.geometry.centroid
    centerIND.to_crs("EPSG:4326")
    
    emisPar = ind[['Hs', 'Ts', 'Ds', 'Vs']].copy()
    isna = emisPar.Hs.isna()
    emisPar.Hs[isna] = 50.0
    isna = emisPar.Ts.isna()
    emisPar.Ts[isna] = 30.0
    emisPar.Ts = emisPar.Ts+273
    isna = emisPar.Ds.isna()
    emisPar.Ds[isna] = 1.0
    isna = emisPar.Vs.isna()
    emisPar.Vs[isna] = 5.0
    
    
    return dataEmissIND, centerIND, emisPar
    

#%% gridSpec
def gridSpecMCIP (rootPath,outPath,
                  mcipMETCRO3Dpath,mcipMETCRO2Dpath,
                  mcipGRIDCRO2DPath,mcipGRIDDOT2DPath):
    
    print('Calling gridSpecMCIP')
    METCRO3D = nc.Dataset(mcipMETCRO3Dpath)
    METCRO2D = nc.Dataset(mcipMETCRO2Dpath)
    GRIDCRO2D = nc.Dataset(mcipGRIDCRO2DPath)        
    GRIDDOT2D = nc.Dataset(mcipGRIDDOT2DPath) 
   
    
    domain,x,y = MCIPDomain(GRIDDOT2D)
    
    dataEmissIND, centerIND, emisPar = emissReader(rootPath,domain)
    
    #Loop over each cel in x direction
    polygons=[]
    for ii in range(1,x.shape[0]):
        #Loop over each cel in y direction
        for jj in range(1,y.shape[0]):
            #roadClip=[]
            lat_point_list = [y[jj-1], y[jj], y[jj], y[jj-1]]
            lon_point_list = [x[ii-1], x[ii-1], x[ii], x[ii]]
            cel = Polygon(zip(lon_point_list, lat_point_list))
            polygons.append(cel)

    # Creating basegridfile
    baseGrid = gpd.GeoDataFrame({'geometry':polygons})
    baseGrid.to_csv(outPath+'/baseGrid'+'.csv')
    baseGrid.crs = "EPSG:4326" 
    print('baseGrid.csv was created at ' + outPath ) 
    
    # Calling speciation function
    print('Speciating emissions')
    dataEmissX = IndSpeciate(rootPath,dataEmissIND)
    
    xv,yv,xX,yY = griddingMCIP(GRIDCRO2D,GRIDDOT2D,METCRO3D)
    
    dataTempo = populatingGridTemp3D(dataEmissX,centerIND,emisPar,xv,yv,xX,yY,
                                     METCRO3D,METCRO2D,GRIDCRO2D)
    
    return dataTempo
