#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                             IND2netCDF.py
                             
                             
This is the main script to convert industrial emission from .csv file 
to chemical speciated  and temporal disaggregated netCDF files.


Inputs: 
    
    rootPath: Path to functions
    
    outPath: Path to output folder
    
    lati: Initial latitude (lower-left)
    
    latf: Final latitude (upper-right)
    
    loni: Initial longitude (lower-left)
    
    lonf: Final longitude (upper-right)
    
    deltaX: Grid resolution/spacing in x direction
    
    deltaY: Grig resolution/spacing in y direction
    
    year: Base-year for your simulation
            
    fileId = identification of your output files
    
    
Input file:
    
    BR_Ind.xlsx - excel file with industrial characteristics, emissions, 
        and coordinates. Emissions are in kg/s.
       

Outputs:
    
    'INDannualEmiss_'+fileId+'_'+str(deltaX)+'x'+str(deltaY)+'_'+str(year)+'.nc'
        (unit in g/year)
    
External functions:
    gridding, populatingGrid, netCDFcreator
    
    
    
Last update = 29/10/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br

-------------------------------------------------------------------------------
"""

import geopandas as gpd
import pandas as pd
import numpy as np
#from netCDFcreator import createNETCDFtemporal
#from temporalDisagregation import temporalDisagVehicular
import numpy.matlib
from shapely.geometry import Polygon
from gridding import gridding, populatingGrid, griddingMCIP, populatingGridTemp3D
from IND_speciate import IndSpeciate
#import datetime
import netCDF4 as nc


#%% userDomain
def userDomain(lati,latf,loni,lonf,deltaX,deltaY):
    # Creating domain window
    domain = Polygon(zip([loni,loni,lonf,lonf],[lati,latf,latf,lati])) 
    domain = gpd.GeoDataFrame(index=[0],geometry=[domain])
    domain.crs = "EPSG:4326"
    print('Setting domain borders')
    x = np.arange(loni, lonf+2*deltaX, deltaX)
    y = np.arange(lati, latf+2*deltaY, deltaY)
    return domain,x,y

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
def gridSpec (rootPath,outPath,dataEmissIND,centerIND,emisPar,x,y):
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
    baseGrid.to_csv(outPath+'/baseGrid.csv')
    baseGrid.crs = "EPSG:4326" 
    print('baseGrid.csv was created at ' + outPath )

   
    
    # Calling speciation function
    dataEmissX = IndSpeciate(rootPath,dataEmissIND)
    
    # Calling gridding function
    grids,xv,yv,xX,yY = gridding(x,y)
    # Calling populatingGrid function
    dataIND = populatingGrid(dataEmissX.fillna(0),centerIND,xX,yY,xv,yv)
    return dataIND

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
# # Setting output file's name
# name = 'INDannualEmiss_'+fileId+'_'+str(deltaX)+'x'+str(deltaY)+'_'+str(year)+'.nc'

# # Calling creatNETCDF function
# createNETCDF(outPath,name,dataIND,xv,yv,y,x,centerIND,year,month)

# name = 'INDannualSpecEmiss_'+fileId+'_'+str(deltaX)+'x'+str(deltaY)+'_'+str(year)+'.nc'
# # Calling createNETCDFtemporal - ANNUAL EMISSIONS
# startDate = datetime.datetime(year, month, 1, 0, 0)
# endDate = datetime.datetime(year, month, 1, 1, 0)
# datePfct = np.arange(np.datetime64(startDate),np.datetime64(endDate),3600000000)
# dates = pd.DataFrame(datePfct)   
# dates['year'] = year
# dates['month'] = 1
# dates['day'] = 1
# dates['hour'] = 00
# # Getting number of days in a year
# fullDate = np.arange(np.datetime64(datetime.datetime(year, 1, 1, 0, 0)),
#                      np.datetime64(datetime.datetime(year, 12, 31, 0, 0)),
#                      3600000000*24)
# # Conversion from second to year
# convSec2year = fullDate.shape[0]*24*60*60 
# createNETCDFtemporal(outPath,name,dataIND*convSec2year,xv,yv,y,x,dates,month,'Annual')

#%% Calling netcdf creator

# # Creating hourly basis inventory - results in grams per seconds
# name = 'INDtemporalEmiss_'+str(year)+'_'+str(month)+'.nc'
# dataMat,datePfct,disvec = temporalDisagIndustrial(dataEmissIND, year, month)
# dataTempo = populatingGridMat(dataMat,centerIND,xX,yY)
# datePfct = disvec
# createNETCDFtemporalIND(folder,name,dataTempo,xv,yv,lat,lon,centerIND,disvec,month)
