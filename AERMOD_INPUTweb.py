# -*- coding: utf-8 -*-
"""
--------------------------AERMOD_INPUT------------------------------

This function is used to create AERMOD input files for point sources

inputs: 
    folder = folder of excel file with point sources data

Update: 25/03/2020

Developed by: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br
Laboratório de Controle da Qualidade do Ar - LCQAr
Universidade Federal de Santa Catarina

"""

#%% Import section
import pandas as pd
#import numpy as np
from django.http import HttpResponse

#%% Input Folder and files
def pointSourceInput(data,user, pollutant, standard, start_year, start_month, start_day, end_year, end_month, end_day,folder):
    data = data.reset_index(drop=True)
    # Opening excel file with input data
    #data = pd.read_excel(file)
    # Creating a aermod.inp file
    #folder = "media/PointSO/"
    
    file2 = open(folder + "/aermod.inp","w") 
    #file1 = HttpResponse()
    #file1.write('')
    file2.write('')
    file2 = open(folder + "/aermod.inp","a") 
    #file1 = HttpResponse()
    
    # Adding control text
    str_CO = 'CO STARTING \nCO TITLEONE '+str(user)+ '\nCO MODELOPT DFAULT CONC \
         \nCO AVERTIME 1 8 24 \nCO POLLUTID  '+pollutant+ ' \nCO FLAGPOLE 1.50 \
             \nCO RUNORNOT RUN \nCO FINISHED \nSO STARTING \nSO ELEVUNIT METERS \n'
    #file1.write(str_CO)
    file2.write(str_CO)

    # Loop over each source ID
    for ii in range (0,data.shape[0]):
        
        # Loop over each pollutant
        
        for jj in range(9, data.shape[1]):
            # Adding source text for each pollutant
            str_SO = 'SO LOCATION ' + data.fonte[ii] +'p'+str(jj-8) + ' ' + data.tipo[ii] + ' ' + \
                str(data.Long[ii]) + ' ' +  str(data.Lat[ii]) + ' ' + \
                    str(data.Altitude[ii]) + '\n'
            #file1.write(str_SO)
            file2.write(str_SO)

        
            str_SOparam = 'SO SRCPARAM ' + data.fonte[ii] +'p'+str(jj-8)+ \
                ' ' + str(data[data.columns[jj]][ii]) + \
                    ' ' + str(data['Altura da chaminé'][ii]) + \
                        ' ' + str(data['Temperatura de saída'][ii]) + \
                            ' ' + str(data['Velocidade de saída'][ii]) + \
                                ' ' + str(data['Diâmetro'][ii]) + '\n'
            #file1.write(str_SOparam)
            file2.write(str_SOparam)

            
     # Loop for grouping sources         
    for jj in range(9, data.shape[1]):
        
        str_sogr=[]
        str_sogrN = ' '
        for ii in range (0,data.shape[0]):
            str_sogr = data.fonte[ii] +'p'+str(jj-8) +' '
            str_sogrN = str_sogrN + str_sogr
            
        str_gr = 'SO SRCGROUP ' + 'POL' + str(jj-8) + str_sogrN  +'\n'
        #file1.write(str_gr)
        file2.write(str_gr)

    
    # Receptor and meteorological information    
    str_rm = '\nSO FINISHED \nRE STARTING \nRE INCLUDED RECEPT.ROU \nRE FINISHED \nME STARTING \
        \nME SURFFILE METEO.SFC \nME PROFFILE METEO.PFL \nME SURFDATA 66666 ' + str(start_year)+ \
            '\nME UAIRDATA 66666 '+ str(start_year)+ '\nME PROFBASE 60 METERS \nME STARTEND ' \
            + str(start_year) + ' ' + str(start_month) + ' ' + str(start_day) + ' ' \
            + str(end_year) + ' '+ str(end_month) + ' '+ str(end_day) + '\nME FINISHED \n'
    #file1.write(str_rm)
    file2.write(str_rm)
    
    # Outputs
    str_OU = 'OU STARTING \nOU RECTABLE ALLAVE 1ST SECOND \nOU MAXTABLE ALLAVE 100 \
        \nOU RANKFILE 1 50 RANK1.RNK \nOU RANKFILE 8 50 RANK8.RNK \
            \nOU RANKFILE 24 50 RANK24.RNK \n'       
    #file1.write(str_OU)  
    file2.write(str_OU)  
    
    for jj in range(9, data.shape[1]):
        str_OU2 = 'OU MAXIFILE 1 ' +  'POL' + str(jj-8) + ' '+str(standard)+' MAX1H_'+'POL' + str(jj-8)+'.OUT \n' +\
            'OU MAXIFILE 8 ' +  'POL' + str(jj-8)  + ' '+str(standard)+' MAX8H_'+'POL' + str(jj-8)+'.OUT \n'+ \
                'OU MAXIFILE 24 ' +  'POL' + str(jj-8)  + ' '+str(standard)+' MAX24H_'+'POL' + str(jj-8)+'.OUT \n'
        
        str_OU3 = 'OU PLOTFILE 1 ' +  'POL' + str(jj-8) + ' FIRST PLOT1H_'+  'POL' + str(jj-8) +'.PLT \n' +\
            'OU PLOTFILE 8 ' +  'POL' + str(jj-8) + ' FIRST PLOT8H_'+  'POL' + str(jj-8) +'.PLT \n' +\
                'OU PLOTFILE 24 ' +  'POL' + str(jj-8) + ' FIRST PLOT24H_'+  'POL' + str(jj-8) +'.PLT \n' 
                 #+\'OU PLOTFILE ANNUAL ' +  'POL' + str(jj-8) + ' FIRST PLOTANNUAL_'+  'POL' + str(jj-8) +'.PLT \n'
        
        #file1.write(str_OU2)
        file2.write(str_OU2)
        #file1.write(str_OU3)
        file2.write(str_OU3)

        
    #file1.write('OU FINISHED')
    file2.write('OU FINISHED')

    file2.close()


    return data

def pointNewSourceInput(file,user, pollutant, standard, start_year, start_month, start_day, end_year, end_month, end_day,folder):
    
    # Opening excel file with input data
    data = pd.read_excel(file)
    # Creating a aermod.inp file
    #folder = "media/PointSO/"
    file2 = open(folder + "aermod.inp","w") 
    file1 = HttpResponse()
    file1.write('')
    file2.write('')
    file2 = open(folder + "aermod.inp","a") 
    file1 = HttpResponse()
    
    # Adding control text
    str_CO = 'CO STARTING \nCO TITLEONE '+str(user)+ '\nCO MODELOPT DFAULT CONC \
         \nCO AVERTIME 1 8 24 \nCO POLLUTID  '+pollutant+ ' \nCO FLAGPOLE 1.50 \
             \nCO RUNORNOT RUN \nCO FINISHED \nSO STARTING \nSO ELEVUNIT METERS \n'
    file1.write(str_CO)
    file2.write(str_CO)

    # Loop over each source ID
    for ii in range (0,data.shape[0]):
        
        # Loop over each pollutant
        
        for jj in range(9, data.shape[1]):
            # Adding source text for each pollutant
            str_SO = 'SO LOCATION ' + data.fonte[ii] +'p'+str(jj-8) + ' ' + data.tipo[ii] + ' ' + \
                str(data.Long[ii]) + ' ' +  str(data.Lat[ii]) + ' ' + \
                    str(data.Altitude[ii]) + '\n'
            file1.write(str_SO)
            file2.write(str_SO)

        
            str_SOparam = 'SO SRCPARAM ' + data.fonte[ii] +'p'+str(jj-8)+ \
                ' ' + str(data[data.columns[jj]][ii]) + \
                    ' ' + str(data['Altura da chaminé'][ii]) + \
                        ' ' + str(data['Temperatura de saída'][ii]) + \
                            ' ' + str(data['Velocidade de saída'][ii]) + \
                                ' ' + str(data['Diâmetro'][ii]) + '\n'
            file1.write(str_SOparam)
            file2.write(str_SOparam)

            
     # Loop for grouping sources         
    for jj in range(9, data.shape[1]):
        
        str_sogr=[]
        str_sogrN = ' '
        for ii in range (0,data.shape[0]):
            str_sogr = data.fonte[ii] +'p'+str(jj-8) +' '
            str_sogrN = str_sogrN + str_sogr
            
        str_gr = 'SO SRCGROUP ' + 'POL' + str(jj-8) + str_sogrN  +'\n'
        file1.write(str_gr)
        file2.write(str_gr)

    
    # Receptor and meteorological information    
    str_rm = '\nSO FINISHED \nRE STARTING \nRE INCLUDED RECEPT.ROU \nRE FINISHED \nME STARTING \
        \nME SURFFILE METEO.SFC \nME PROFFILE METEO.PFL \nME SURFDATA 66666 ' + str(start_year)+ \
            '\nME UAIRDATA 66666 '+ str(start_year)+ '\nME PROFBASE 60 METERS \nME STARTEND ' \
            + str(start_year) + ' ' + str(start_month) + ' ' + str(start_day) + ' ' \
            + str(end_year) + ' '+ str(end_month) + ' '+ str(end_day) + '\nME FINISHED \n'
    file1.write(str_rm)
    file2.write(str_rm)
    
    # Outputs
    str_OU = 'OU STARTING \nOU RECTABLE ALLAVE 1ST SECOND \nOU MAXTABLE ALLAVE 100 \
        \nOU RANKFILE 1 50 RANK1.RNK \nOU RANKFILE 8 50 RANK8.RNK \
            \nOU RANKFILE 24 50 RANK24.RNK \n'       
    file1.write(str_OU)  
    file2.write(str_OU)  
    
    for jj in range(9, data.shape[1]):
        str_OU2 = 'OU MAXIFILE 1 ' +  'POL' + str(jj-8) + ' '+str(standard)+' MAX1H_'+'POL' + str(jj-8)+'.OUT \n' +\
            'OU MAXIFILE 8 ' +  'POL' + str(jj-8)  + ' '+str(standard)+' MAX8H_'+'POL' + str(jj-8)+'.OUT \n'+ \
                'OU MAXIFILE 24 ' +  'POL' + str(jj-8)  + ' '+str(standard)+' MAX24H_'+'POL' + str(jj-8)+'.OUT \n'
        
        str_OU3 = 'OU PLOTFILE 1 ' +  'POL' + str(jj-8) + ' FIRST PLOT1H_'+  'POL' + str(jj-8) +'.PLT \n' +\
            'OU PLOTFILE 8 ' +  'POL' + str(jj-8) + ' FIRST PLOT8H_'+  'POL' + str(jj-8) +'.PLT \n' +\
                'OU PLOTFILE 24 ' +  'POL' + str(jj-8) + ' FIRST PLOT24H_'+  'POL' + str(jj-8) +'.PLT \n' 
                 #+\'OU PLOTFILE ANNUAL ' +  'POL' + str(jj-8) + ' FIRST PLOTANNUAL_'+  'POL' + str(jj-8) +'.PLT \n'
        
        file1.write(str_OU2)
        file2.write(str_OU2)
        file1.write(str_OU3)
        file2.write(str_OU3)

        
    file1.write('OU FINISHED')
    file2.write('OU FINISHED')

    file2.close()


    return data, file1

def pointNewAndOldSourceInput(file,dataOld,user, pollutant, standard, start_year, start_month, start_day, end_year, end_month, end_day,folder):
    
    # Opening excel file with input data
    
    data = pd.read_excel(file)
    data = pd.concat([data,dataOld])
    # Creating a aermod.inp file
    #folder = "media/PointSO/"
    file2 = open(folder + "aermod.inp","w") 
    file1 = HttpResponse()
    file1.write('')
    file2.write('')
    file2 = open(folder + "aermod.inp","a") 
    file1 = HttpResponse()
    
    # Adding control text
    str_CO = 'CO STARTING \nCO TITLEONE '+str(user)+ '\nCO MODELOPT DFAULT CONC \
         \nCO AVERTIME 1 8 24 \nCO POLLUTID  '+pollutant+ ' \nCO FLAGPOLE 1.50 \
             \nCO RUNORNOT RUN \nCO FINISHED \nSO STARTING \nSO ELEVUNIT METERS \n'
    file1.write(str_CO)
    file2.write(str_CO)

    # Loop over each source ID
    for ii in range (0,data.shape[0]):
        
        # Loop over each pollutant
        
        for jj in range(9, data.shape[1]):
            # Adding source text for each pollutant
            str_SO = 'SO LOCATION ' + data.fonte[ii] +'p'+str(jj-8) + ' ' + data.tipo[ii] + ' ' + \
                str(data.Long[ii]) + ' ' +  str(data.Lat[ii]) + ' ' + \
                    str(data.Altitude[ii]) + '\n'
            file1.write(str_SO)
            file2.write(str_SO)

        
            str_SOparam = 'SO SRCPARAM ' + data.fonte[ii] +'p'+str(jj-8)+ \
                ' ' + str(data[data.columns[jj]][ii]) + \
                    ' ' + str(data['Altura da chaminé'][ii]) + \
                        ' ' + str(data['Temperatura de saída'][ii]) + \
                            ' ' + str(data['Velocidade de saída'][ii]) + \
                                ' ' + str(data['Diâmetro'][ii]) + '\n'
            file1.write(str_SOparam)
            file2.write(str_SOparam)

            
     # Loop for grouping sources         
    for jj in range(9, data.shape[1]):
        
        str_sogr=[]
        str_sogrN = ' '
        for ii in range (0,data.shape[0]):
            str_sogr = data.fonte[ii] +'p'+str(jj-8) +' '
            str_sogrN = str_sogrN + str_sogr
            
        str_gr = 'SO SRCGROUP ' + 'POL' + str(jj-8) + str_sogrN  +'\n'
        file1.write(str_gr)
        file2.write(str_gr)

    
    # Receptor and meteorological information    
    str_rm = '\nSO FINISHED \nRE STARTING \nRE INCLUDED RECEPT.ROU \nRE FINISHED \nME STARTING \
        \nME SURFFILE METEO.SFC \nME PROFFILE METEO.PFL \nME SURFDATA 66666 ' + str(start_year)+ \
            '\nME UAIRDATA 66666 '+ str(start_year)+ '\nME PROFBASE 60 METERS \nME STARTEND ' \
            + str(start_year) + ' ' + str(start_month) + ' ' + str(start_day) + ' ' \
            + str(end_year) + ' '+ str(end_month) + ' '+ str(end_day) + '\nME FINISHED \n'
    file1.write(str_rm)
    file2.write(str_rm)
    
    # Outputs
    str_OU = 'OU STARTING \nOU RECTABLE ALLAVE 1ST SECOND \nOU MAXTABLE ALLAVE 100 \
        \nOU RANKFILE 1 50 RANK1.RNK \nOU RANKFILE 8 50 RANK8.RNK \
            \nOU RANKFILE 24 50 RANK24.RNK \n'       
    file1.write(str_OU)  
    file2.write(str_OU)  
    
    for jj in range(9, data.shape[1]):
        str_OU2 = 'OU MAXIFILE 1 ' +  'POL' + str(jj-8) + ' '+str(standard)+' MAX1H_'+'POL' + str(jj-8)+'.OUT \n' +\
            'OU MAXIFILE 8 ' +  'POL' + str(jj-8)  + ' '+str(standard)+' MAX8H_'+'POL' + str(jj-8)+'.OUT \n'+ \
                'OU MAXIFILE 24 ' +  'POL' + str(jj-8)  + ' '+str(standard)+' MAX24H_'+'POL' + str(jj-8)+'.OUT \n'
        
        str_OU3 = 'OU PLOTFILE 1 ' +  'POL' + str(jj-8) + ' FIRST PLOT1H_'+  'POL' + str(jj-8) +'.PLT \n' +\
            'OU PLOTFILE 8 ' +  'POL' + str(jj-8) + ' FIRST PLOT8H_'+  'POL' + str(jj-8) +'.PLT \n' +\
                'OU PLOTFILE 24 ' +  'POL' + str(jj-8) + ' FIRST PLOT24H_'+  'POL' + str(jj-8) +'.PLT \n' 
                 #+\'OU PLOTFILE ANNUAL ' +  'POL' + str(jj-8) + ' FIRST PLOTANNUAL_'+  'POL' + str(jj-8) +'.PLT \n'
        
        file1.write(str_OU2)
        file2.write(str_OU2)
        file1.write(str_OU3)
        file2.write(str_OU3)

        
    file1.write('OU FINISHED')
    file2.write('OU FINISHED')

    file2.close()


    return data, file1