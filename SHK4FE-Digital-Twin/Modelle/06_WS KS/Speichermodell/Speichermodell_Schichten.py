# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 10:15:05 2023

@author: leroytomas
"""

import numpy as np
from scipy.integrate import odeint # für die Integration!!!

#%%
""" Speicherparameter berechnen: """

zi = Par_TES_htank / Par_TES_N
Ai = (math.pi*(Par_TES_dtank - 2*Par_TES_thtank)**2) / 4
Aexti = math.pi*Par_TES_dtank*zi # exterior heat transfer surface area of a layer
mi = Ai*zi*Prop_Water_Rho

#%%

m_dot = m_dot_HP - m_dot_load # when positive, then TES is being charged by the HP. When negative the TES is being discharged. Charge and discharged dont happen at the same time. It is just the result of mdot what determines what is happening.

def TES(Par_TES_N, Tfs, Trl, m_dot, TES_top, TES_bottom, Tamb):
                      
    Ti = np.linspace(TES_bottom, TES_top, Par_TES_N).tolist() # set intial top and bottom temps using last value
    t = np.linspace(0, 900)  # Time Grid (900 s)

    """ calculate the bottom layer first--------------------------------------------------"""
    def energybalance_1(y, t): # Function for the bottom layer
        Ti[i] = y

        if m_dot>0:                
            dTdt = (
              ((Ai*lambda_eff/zi) * (Ti[i+1] - Ti[i])) \
            + (m_dot * Prop_Water_Cp * (Ti[i+1] - Ti[i])) \
            - (k * Aexti * (Ti[i] - Tamb))) / (mi*Prop_Water_Cp)
                                               
        else: #The model is still not working properly. If a +X is not added to T[i] then the bottom temperature is stucked in one value and does not allow to discharge the tank. Leroy will check it.           
            dTdt = ((
              ((Ai*lambda_eff/zi)*(Ti[i] - Ti[i+1])) \
            + (m_dot * Prop_Water_Cp * (Ti[i] - Trl)) \
            - (k*Aexti*(Ti[i] - Tamb))) / (mi*Prop_Water_Cp))
                
        return dTdt

    i = 0 #Bottom layer index
    y0 = Ti[i]  # initial condition  
    Ti[i] = float(odeint(energybalance_1, y0, t)[-1]) # solve ODE

    """ calculate and iterate through the middle layers------------------------------------- """
    
    def energybalance_2(y, t): # Function for the middle layer
        Ti[i] = y
        
        if m_dot>0:
            dTdt = ((((Ai*lambda_eff)/zi) * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) \
                    + (m_dot*Prop_Water_Cp*(Ti[i+1] - Ti[i])) \
                    - (k*Aexti*(Ti[i] - Tamb))) / (mi*Prop_Water_Cp)
            
        else:
            dTdt = ((((Ai*lambda_eff)/zi) * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) \
                    + (m_dot*Prop_Water_Cp*(Ti[i] - Ti[i-1])) \
                    - (k*Aexti*(Ti[i] - Tamb))) / (mi*Prop_Water_Cp)
        
        return dTdt

    for i in range(1, Par_TES_N - 1): #All middle layers index
        if  m_dot>0:
            y0 = Ti[i+1] # take layer above the current layer as initial condition          
        else:
            y0 = Ti[i-1] # take layer under the current layer as initial condition          

        Ti[i] = float(odeint(energybalance_2, y0, t)[-1]) # solve ODE
        
    """ calculate the top layer --------------------------------------------------------------"""
    def energybalance_3(y, t): # Function for the top layer
        Ti[i] = y
        
        if m_dot>0:                
            dTdt = (
              ((Ai*lambda_eff/zi)*(Ti[i-1] - Ti[i])) \
            + (m_dot * Prop_Water_Cp * (Tfs - Ti[i])) \
            - (k*Aexti*(Ti[i] - Tamb))) / (mi*Prop_Water_Cp)                                 
        else:       
            dTdt = (
              ((Ai*lambda_eff/zi)*(Ti[i] - Ti[i-1])) \
            + (m_dot * Prop_Water_Cp * (Ti[i] - Ti[i-1])) \
            - (k*Aexti*(Ti[i] - Tamb))) / (mi*Prop_Water_Cp)
                            
        return dTdt
    
    i = Par_TES_N - 1 # Top layer index
    y0 = Ti[i] #initial condition 
    Ti[i] = float(odeint(energybalance_3, y0, t)[-1]) # solve ODE
        
    return Ti  


""" so kann das Aufrufen von der Funktion dann aussehen: """

TES_theta_list = TES(Par_TES_N = Par_TES_N,
                     Tfs = HP_theta_Con_Out_list[-1],         # feed temperature from source (HP). 
                     TES_top = TES_theta_Top_list[-1],        # feed temperature to load (TABS) 
                     TES_bottom = TES_theta_Bottom_list[-1],  # return temperature to source (HP)      
                     Trl = TABS_thetaRL_Last_list[-1]-6,      # return temperature from load (TABS). It has to be lower than Vorlauftemperature otherwise the model with just stagnate in one temperature
                     m_dot = m_dot_list[-1],
                     Tamb = theta_op_list_EG[i]) #The heat pump is located in the Ergeschoss

