#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                             IND_speciate.py
                             
                             
This function speciates the industrial emission in chemical species. Conversion
from NOx to NO (0.495) and NO2 (0.505)

Inputs: 
    
    rootPath: Path to functions
    
    dataEmissIND: dataframe with industrial emissions and coordinates
    
    
Input files:
    
    CMAQ_species.csv: list of CMAQ species in IndustriaSpeciation folder
        
    IND_PROFILES_speciate.csv: Speciation profiles by industrial ID in 
            IndustriaSpeciation folder
        
    BR_Ind_Profiles.csv: speciate profiles to be used in each industrial ID
   
    CMAQ_speciesMW: Molecular Weight of chemical species    
       

Outputs:
    
    dataEmissX: Dataframe with speciated emissions

    

Last update = 29/10/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br

-------------------------------------------------------------------------------
"""

#%% Importing packages

import pandas as pd
import numpy as np


#%% ------------------------- PROCESSING---------------------------------------
def IndSpeciate (rootPath,dataEmissIND):
    
    # Opening CMAQ_species.csv
    file_path_CMAQ_species = "/IndustrialSpeciation/CMAQ_species.csv"
    dfCMAQspc = pd.read_csv(rootPath+file_path_CMAQ_species)
    
    # Opening SPEC_names.csv
    file_path_SPEC_names = "/IndustrialSpeciation/SPEC_names.csv"
    dfspec = pd.read_csv(rootPath+file_path_SPEC_names,index_col=0)
    
    # Opening SPEC_names.csv
    file_path_IND_PROFILES_speciate = "/Inputs/IND_PROFILES_speciate.xlsx"
    proConv = pd.read_excel(rootPath+file_path_IND_PROFILES_speciate,index_col=0)
    proConv['proID'] = proConv.index.astype(str)
    
    # Opening SPEC_names.csv
    file_path_BR_Ind_Profiles = "/Inputs/BR_Ind_Profiles.xlsx"
    IndProfiles = pd.read_excel(rootPath+file_path_BR_Ind_Profiles,index_col=0)
    
    file_path_CMAQ_speciesMW = '/IndustrialSpeciation/CMAQ_speciesMW.csv'
    smm = pd.read_csv(rootPath+file_path_CMAQ_speciesMW)
    
    dfCMAQspc2 = pd.DataFrame()
    for idx in IndProfiles.index:
        if isinstance(IndProfiles['PROFILE'][idx], str):
            proPM = IndProfiles['PROFILE'][idx].split(";")
            for idxPM in proPM:
                dfPM = pd.read_csv(rootPath+'/IndustrialSpeciation/'+idxPM.replace(" ", "")+'.csv',index_col=0)
                dfPJoin = dfPM.join(dfspec)
                dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
                dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
                dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
                dfCMAQspcX['IDX'] = idx
                
                # Converting TOG to VOC
                if np.isnan(proConv[proConv['proID'].str.replace(" ", "") == 
                            idxPM.replace(" ", "")]['TOG_to_VOC RATIO'].values)==False:
                    dfCMAQspcX['WEIGHT_PERCENT'] = dfCMAQspcX['WEIGHT_PERCENT']*\
                        proConv[proConv['proID'].str.replace(" ", "") == 
                                idxPM.replace(" ", "")]['TOG_to_VOC RATIO'].values
                   
                dfCMAQspc2 = pd.concat([dfCMAQspc2,dfCMAQspcX])
                              
            
    dfCMAQspc2['Formula'] = dfCMAQspc2.index.get_level_values(0)
    dfCMAQspc2['ID'] = dfCMAQspc2['IDX']
    dfCMAQspc2 = dfCMAQspc2.groupby(by=['IDX','Formula']).mean()
    
    
    df3 = dfCMAQspc2.unstack(level='IDX')
    df3.columns = df3.columns.droplevel()
    cmaqSpecies = pd.DataFrame()
    cmaqSpecies['Formula']=df3.index
    
    for idx in IndProfiles.index:
        if any(dfCMAQspc2['ID']==idx):
            cmaqSpecies['ID_'+str(idx)] = np.array(dfCMAQspc2[dfCMAQspc2['ID']==idx]['WEIGHT_PERCENT'])
        else:
            cmaqSpecies['ID_'+str(idx)] = np.zeros((cmaqSpecies.shape[0],1))


    # Creating the .csv file   
    cmaqSpecies = dfCMAQspc.set_index('Formula').merge(cmaqSpecies.set_index('Formula'),how='left', on='Formula')
    cmaqSpecies=cmaqSpecies.fillna(0)   
    cmaqSpecies.iloc[:,3:]=cmaqSpecies.iloc[:,3:]/100 # converting from percentage to factor
    cmaqSpecies.to_csv(rootPath+'/IndustrialSpeciation/IND_speciation.csv')
    
    
    dataEmissIND.columns = dataEmissIND.columns.str.replace(' ', '')
    # Creating dataEmiss matrix
    dataEmissX=pd.DataFrame()

    vocSpecs = ['ACET','ACROLEIN','ALD2','ALD2_PRIMARY','BENZ','BUTADIENE13','CH4',
                'CH4_INV','ETH','ETHA','ETHY','ETOH','FORM','FORM_PRIMARY','ISO',
                'MEOH','NAPH','PRPA','TOL','XYLMN']
                

    # Filling VOC emissions
    for ii in range(0,dfCMAQspc.shape[0]):
        for jj in range(0,dataEmissIND['ID'].shape[0]):
            col = 'ID_'+str(int(dataEmissIND['ID'][jj]))
            dataEmissX[dfCMAQspc.ID[ii]]= np.array(cmaqSpecies[col][ii])*dataEmissIND['PMemis']/smm.iloc[ii,1] 
    
    for ii in range(0,len(vocSpecs)):
        for jj in range(0,dataEmissIND['ID'].shape[0]):
            col = 'ID_'+str(int(dataEmissIND['ID'][jj]))
            dataEmissX[vocSpecs[ii]][jj]= \
                np.array(cmaqSpecies[col][cmaqSpecies.ID==vocSpecs[ii]])*\
                    dataEmissIND['VOCemis'][jj]/smm[smm.ID==vocSpecs[ii]].MM

    dataEmissX['CO'] = dataEmissIND['COemis']/np.array(smm[smm.iloc[:,0]=='CO'].MM)
    dataEmissX['SO2'] = dataEmissIND['SOxemis']/np.array(smm[smm.iloc[:,0]=='SO2'].MM)
    dataEmissX['VOC_INV'] = dataEmissIND['VOCemis']
    dataEmissX['PMC'] = dataEmissIND['PMemis']
    dataEmissX['NO'] = dataEmissIND['NOxemis']*np.array(0.495/smm[smm.iloc[:,0]=='NO'].MM)
    dataEmissX['NO2'] = dataEmissIND['NOxemis']*np.array(0.505/smm[smm.iloc[:,0]=='NO2'].MM)
    
    return dataEmissX
