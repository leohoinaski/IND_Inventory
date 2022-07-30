#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 12:06:47 2022

@author: leohoinaski
"""

import numpy as np

def SMOKEplumeRise(Ts,Vs,Ds,Hs):
    G=9.80665 #Mean gravitational acceleration
    U=2.0 #Default wind speed
    T=293.0 #Default ambient air temperature
    Ts = np.array(Ts)
    Vs = np.array(Vs)
    Ds = np.array(Ds, dtype=float)
    Hs = np.array(Hs)
    
    F=0.25*G*Vs*(Ds*Ds)*((Ts-T)/Ts)
    if F<0:
        F=0
    if F<55:
        Hef = Hs + 21.31311057*(F**0.75)/U
    else:
        Hef = Hs + 38.87776061*(F**0.6)/U 
        
    return Hef
               

def BRIGGSplumeRise(Ts,Vs,Ds,Hs,Tas,Uas,dTdZ):
    G=9.80665 #Mean gravitational acceleration
    Ts = np.array(Ts)
    Vs = np.array(Vs)
    Ds = np.array(Ds, dtype=float)
    Hs = np.array(Hs)
    Uas = np.array(Uas)
    
    Fb=G*Vs*(Ds**2)*((Ts-Tas)/Ts)
    Fb = np.max(Fb,0)
    Fm = (Vs**2)*(Ds**2)*Tas/(4*Ts) 
    
    if dTdZ>-0.0098:
        dTethadZ = dTdZ+0.0098
        s = (G/Tas)*dTethadZ
        deltaTc = 0.019582*Ts*Vs*np.sqrt(s)
        
        if deltaTc>dTdZ:
            dH1 = 1.5*(Fm/(Uas*np.sqrt(s)))**(1/3)
            dH2 = 3*Ds*Vs/Uas
            Hef = Hs + np.min([dH1,dH2])
        else:
            Hef = Hs + 2.6*(Fb/(Uas*s))**(1/3)
    
    else:          
        if Fb<55:
            deltaTc = 0.0297*Ts*(Vs**(1/3))/(Ds**(2/3))
            
            if deltaTc>dTdZ:
                Hef = Hs+3*Ds/(Vs*Uas)
            else:
                Hef = Hs + 21.31311057*(Fb**0.75)/Uas
        else:
            deltaTc = 0.00575*Ts*(Vs**(1/3))/(Ds**(1/3))
            
            if deltaTc>dTdZ:
                Hef = Hs+3*Ds/(Vs*Uas)
            else:
                Hef = Hs + 38.87776061*(Fb**0.6)/Uas 
        
    return Hef