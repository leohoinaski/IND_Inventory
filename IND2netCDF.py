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
import os 
import numpy as np
from netCDFcreator import createNETCDFtemporal
#from temporalDisagregation import temporalDisagVehicular
import numpy.matlib
from shapely.geometry import Polygon
from gridding import gridding, populatingGrid
from IND_speciate import IndSpeciate
import datetime


#%%============================= INPUTS =========================================

rootPath = '/media/leohoinaski/HDD/IND_inventory'

outPath = rootPath +"/Outputs"

#-------------------------Setting grid resolution------------------------------

# Users can change the domain and resolution here.
lati =-30 #lati = int(round(bound.miny)) # Initial latitude (lower-left)

latf = -24 #latf = int(round(bound.maxy)) # Final latitude (upper-right)

loni = -54 #loni = int(round(bound.minx)) # Initial longitude (lower-left)

lonf = -47 #lonf = int(round(bound.maxx)) # Final longitude (upper-right)

deltaX = 0.05 # Grid resolution/spacing in x direction

deltaY = 0.05 # Grig resolution/spacing in y direction

fileId = 'SC' # Code to identify your output files

prefix = str(deltaX)+'x'+str(deltaY) # grid definition identification

year = 2019

month = 1

convKg2g = 1000 # Conversion from kg/s to g/s


#%%------------------------- PROCESSING ---------------------------------------

#-------------------- Reading industrial emissions-----------------------------

print('Reading indutrial emissions from ' + rootPath+'/Inputs/BR_Ind.xlsx')

dfind = pd.read_excel(rootPath+'/Inputs/BR_Ind.xlsx')

# Converting to geodataframe
ind = gpd.GeoDataFrame(
    dfind, geometry=gpd.points_from_xy(dfind.Long, dfind.Lat))
ind.crs = "EPSG:4326" # Setting EPSG

domain = Polygon(zip([loni,loni,lonf,lonf],[lati,latf,latf,lati]))
domain = gpd.GeoDataFrame(index=[0],geometry=[domain])
domain.crs = "EPSG:4326"
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

# Calling speciation function
dataEmissX = IndSpeciate(rootPath,dataEmissIND)

# -----------------------Gridding and populating-------------------------------
# cd to the main folder
os.chdir(rootPath)

# Creating output directory
if os.path.isdir(outPath)==0:
    os.mkdir(outPath)
     
print('Setting domain borders')
x = np.linspace(loni, lonf, int((lonf-loni)/deltaX))
y = np.linspace(lati, latf, int((latf-lati)/deltaY))

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
baseGrid.to_csv(outPath+'/baseGrid_'+prefix+'.csv')
baseGrid.crs = "EPSG:4326" 
print('baseGrid_'+prefix+'.csv was created at ' + outPath )


# Calling gridding function
grids,xv,yv,xX,yY = gridding(x,y)

# Calling populatingGrid function
dataIND = populatingGrid(dataEmissX.fillna(0),centerIND,xX,yY,xv,yv)

# # Setting output file's name
# name = 'INDannualEmiss_'+fileId+'_'+str(deltaX)+'x'+str(deltaY)+'_'+str(year)+'.nc'

# # Calling creatNETCDF function
# createNETCDF(outPath,name,dataIND,xv,yv,y,x,centerIND,year,month)

name = 'INDannualSpecEmiss_'+fileId+'_'+str(deltaX)+'x'+str(deltaY)+'_'+str(year)+'.nc'
# Calling createNETCDFtemporal - ANNUAL EMISSIONS
startDate = datetime.datetime(year, month, 1, 0, 0)
endDate = datetime.datetime(year, month, 1, 1, 0)
datePfct = np.arange(np.datetime64(startDate),np.datetime64(endDate),3600000000)
dates = pd.DataFrame(datePfct)   
dates['year'] = year
dates['month'] = 1
dates['day'] = 1
dates['hour'] = 00
# Getting number of days in a year
fullDate = np.arange(np.datetime64(datetime.datetime(year, 1, 1, 0, 0)),
                     np.datetime64(datetime.datetime(year, 12, 31, 0, 0)),
                     3600000000*24)
# Conversion from second to year
convSec2year = fullDate.shape[0]*24*60*60 
createNETCDFtemporal(outPath,name,dataIND*convSec2year,xv,yv,y,x,dates,month,'Annual')

#%% Calling netcdf creator

# # Creating hourly basis inventory - results in grams per seconds
# name = 'INDtemporalEmiss_'+str(year)+'_'+str(month)+'.nc'
# dataMat,datePfct,disvec = temporalDisagIndustrial(dataEmissIND, year, month)
# dataTempo = populatingGridMat(dataMat,centerIND,xX,yY)
# datePfct = disvec
# createNETCDFtemporalIND(folder,name,dataTempo,xv,yv,lat,lon,centerIND,disvec,month)
