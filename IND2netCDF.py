#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------IND2netCDF-----------------------------
Este scrip converte os dados das emissões do industriais provenientes
do BRAIN em matricial para posterior criação do arquivo
netCDF pela função netCDFcreator.py
A resolução temporal é horária durante um mês de dados. 
O usuário precisa indicar o mês e ano. A resolução espacial 
é a mesma do arquivo do BRAVES. 

Inputs:
    folder: pasta com os arquivos do BRAIN
    Inputs das emissões da queima da biomassa em kg/hora
    FALTA ATUALIZAR A GRADE PARA O ARQUIVO DO MCIP

Outputs:
    Arquivos netCDF - INDtemporalEmiss_ e INDannualEmiss.nc
    Outputs em g/s ou mol/s
    

@author: leohoinaski - leonardo.hoinaski@ufsc.br
Atualizado em 04/05/2021
---------------------------------------------------------------
"""

import geopandas as gpd
import pandas as pd
import os 
import numpy as np
from netCDFcreator import createNETCDF
#from temporalDisagregation import temporalDisagVehicular
import numpy.matlib
from shapely.geometry import Polygon
from gridding import gridding, populatingGrid


#%%============================= INPUTS =========================================

rootPath = '/media/leohoinaski/HDD/IND_inventory'
outPath = rootPath +"/Outputs"

#-------------------------Setting grid resolution------------------------------

# Users can change the domain and resolution here.
lati =-30 #lati = int(round(bound.miny)) # Initial latitude (lower-left)

latf = -24 #latf = int(round(bound.maxy)) # Final latitude (upper-right)

loni = -54 #loni = int(round(bound.minx)) # Initial longitude (lower-left)

lonf = -47 #lonf = int(round(bound.maxx)) # Final longitude (upper-right)

deltaX = 0.01 # Grid resolution/spacing in x direction

deltaY = 0.01 # Grig resolution/spacing in y direction

fileId = 'SC' # Code to identify your output files

prefix = str(deltaX)+'x'+str(deltaY) # grid definition identification

year = 2019

month = 1

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

# Estimated data
dataEmissINDestmeas = ind[['PMestmeas', 
                    'COestmeas', 
                    'NOxestmeas',
                    'SOxestmeas',
                    'VOCestmeas']].copy() 

# Substituting nan by estimated data
isna = dataEmissIND.isna()
dataEmissIND[isna] = dataEmissINDestmeas[isna]

# Getting emissions centroid
centerIND = ind.geometry.centroid
centerIND.to_crs("EPSG:4326")   


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
dataIND = populatingGrid(dataEmissIND,centerIND,xv,yv,xX,yY)

# Setting output file's name
name = 'INDannualEmiss_'+fileId+'_'+str(deltaX)+'x'+str(deltaY)+'_'+str(year)+'.nc'

# Calling creatNETCDF function
createNETCDF(outPath,name,dataIND,xv,yv,y,x,centerIND,year,month)


#%% Calling netcdf creator

# # Creating hourly basis inventory - results in grams per seconds
# name = 'INDtemporalEmiss_'+str(year)+'_'+str(month)+'.nc'
# dataMat,datePfct,disvec = temporalDisagIndustrial(dataEmissIND, year, month)
# dataTempo = populatingGridMat(dataMat,centerIND,xX,yY)
# datePfct = disvec
# createNETCDFtemporalIND(folder,name,dataTempo,xv,yv,lat,lon,centerIND,disvec,month)
