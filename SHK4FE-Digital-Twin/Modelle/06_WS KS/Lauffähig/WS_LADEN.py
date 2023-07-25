"""
@author: MoBueh
"""

#%% Libraries
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np
import copy
import math
from scipy.integrate import odeint
plt.rcParams["figure.dpi"] = 400
plt.rc('axes', axisbelow=True)

#%% Parameter
prop = {}
prop['wasser'] = {
    'rho'       : 980,                                                          # Density in kg/m³
    'cp'        : 4180                                                          # Specific heat capacity in J/kg/K
    }
prop['Antifrogen'] = {
    'rho'       : 1026,                                                         # Density in kg/m³
    'cp'        : 3500                                                          # Specific heat capacity in J/kg/K   
    }

param = {}
param['WS'] = {
    'H'             : 1.2,      # m         Speicherhöhe
    'N'             : 7,        # -         Anzahl Schichten
    'D'             : 0.6,      # m         Speicherdurchmesser
    'th'            : 0.02,      # m²        Wanddicke
    'k'             : 0.02,    # W/m²K     Wärmetransferkoeffizient
    'lambda_eff'    : 0.15    # W/mK      vertikale Wärmeleitfähigkeit von Wasser     
    }




#%%
def WS(Vpkt_prim, T_prim_in, Vpkt_sek, T_sek_in, T_vektor_prev, Tamb):
    # "Development and experimental evaluation of grey-box models of a microscale polygeneration system for application in optimal control" 
    # https://www.sciencedirect.com/science/article/abs/pii/S0378778819331342
    
    #                        _______________
    # Vpkt_prim, T_prim_in  |               | T_sek_out
    # --------------------->|_______________|--------------------->
    #                       |               |
    #                       |_______________|
    #                       |               |
    #                       |_______________|
    #            T_prim_out |               | Vpkt_sek, T_sek_in
    # <---------------------|_______________|<---------------------   
    
    # -------------------Parameter-Massenströme-Randbedingungen----------------------------------------------------
    z_i = param['WS']['H']/param['WS']['N']                                     # (25) Höhe einer Schicht
    A_ext = np.pi*param['WS']['D']*z_i                                          # (26) Außenfläche einer Schicht
    A_i = np.pi*(param['WS']['D']-2*param['WS']['th'])**2/4                     # (27) Fläche zwischen zwei Layern
    m_i = A_i*z_i*prop['wasser']['rho']                                         # (28) Masse pro Schicht
    mpkt_prim = Vpkt_prim*prop['wasser']['rho']
    mpkt_sek = Vpkt_sek*prop['wasser']['rho']
    mpkt_i = mpkt_prim-mpkt_sek    
    Ti = copy.deepcopy(T_vektor_prev)
    # -------------------------------------------------------------------------------------------------------------
    # Die folgenden Gleichungen entsprechen Gleichung (29) 
    # -------------------------------------------------------------------------------------------------------------
    # ------------------------------BOTTOM-Layer-------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------
    def EB_BOTTOM(y, t):
    # -----------------------------------DGL-----------------------------------------------------------------------
        Ti[0] = y
        if mpkt_i>=0:                                                           # Resultierender Massenstrom ↓           
            dTdt = (0
                    -mpkt_sek*prop['wasser']['cp']*(Ti[0] - T_sek_in)           # Sek. mpkt-bedingte Änderung
                    -(param['WS']['k'] * A_ext * (Ti[0] - Tamb))                # Wärmeverluste an Umgebung
                    # + (mpkt_i * prop['wasser']['cp'] * (Ti[i+1] - Ti[i]))      # Stofftransport
                    + (mpkt_prim * prop['wasser']['cp'] * (Ti[0+1] - Ti[0]))    # Stofftransport
                    +((A_i*param['WS']['lambda_eff']/z_i) * (Ti[0+1] - Ti[0]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt

        else:                                                                   # Resultierender Massenstrom ↑    
            dTdt = ((
              ((A_i*param['WS']['lambda_eff']/z_i)*(Ti[i] - Ti[i+1])) \
            + (mpkt_i * prop['wasser']['cp'] * (Ti[i] - T_sek_in)) \
            - (param['WS']['k']*A_ext*(Ti[i] - Tamb))) / (m_i*prop['wasser']['cp']))       
        return dTdt
    # -------------------------------Solving-DGL-------------------------------------------------------------------
    Ti[0] = float(odeint(EB_BOTTOM, Ti[0], np.linspace(0, 900))[-1]) # solve ODE
    # ------------------------------------------------------------------------------------------------------------- 
    
    
    # -------------------------------------------------------------------------------------------------------------
    # --------------------------------TOP-Layer--------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------
    def EB_TOP(y, t):
    # -----------------------------------DGL-----------------------------------------------------------------------
        Ti[-1] = y
        if mpkt_i>=0:                                                            # Resultierender Massenstrom ↓            
            dTdt = (0
                    +(mpkt_prim * prop['wasser']['cp'] * (T_prim_in - Ti[-1]))  # Prim. mpkt-bedingte Änderung
                    -(param['WS']['k']*A_ext*(Ti[-1] - Tamb))                   # Wärmeverluste an Umgebung 
                    +0 # keine Schicht oberhalb vorhanden                       # Stofftransport
                    +((A_i*param['WS']['lambda_eff']/z_i)*(Ti[-1-1] - Ti[-1]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt  
        
        else:                                                                   # Resultierender Massenstrom ↑   
            dTdt = (
              ((A_i*param['WS']['lambda_eff']/z_i)*(Ti[i] - Ti[i-1])) \
            + (mpkt_i * prop['wasser']['cp'] * (Ti[i] - Ti[i-1])) \
            - (param['WS']['k']*A_ext*(Ti[i] - Tamb))) / (m_i*prop['wasser']['cp'])
        return dTdt
    # -------------------------------Solving-DGL-------------------------------------------------------------------
    Ti[-1] = float(odeint(EB_TOP, Ti[-1], np.linspace(0, 900))[-1]) # solve ODE
    # -------------------------------------------------------------------------------------------------------------
    
    
    # -------------------------------------------------------------------------------------------------------------
    # ------------------------------MIDDLE-Layer-------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------
    def EB_MIDDLE(y, t): # Function for the middle layer
        Ti[i] = y
    # -----------------------------------DGL-----------------------------------------------------------------------
        if mpkt_i>=0:                                                           # Resultierender Massenstrom ↓ 
            dTdt = (0
                    +0 # kein Anschluss in diesen Layern                        # Sek/prim mpkt-bedingte Änderung
                    -(param['WS']['k']*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung
                    +(mpkt_i*prop['wasser']['cp']*(Ti[i+1] - Ti[i]))            # Stofftransport
                    +((A_i*param['WS']['lambda_eff'])/z_i * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt
   
        else:                                                                   # Resultierender Massenstrom ↑ 
            dTdt = ((((A_i*param['WS']['lambda_eff'])/z_i) * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) \
                    + (mpkt_i*prop['wasser']['cp']*(Ti[i] - Ti[i-1])) \
                    - (param['WS']['k']*A_ext*(Ti[i] - Tamb))) / (m_i*prop['wasser']['cp'])
        return dTdt
    # -------------------------------Solving-DGL-------------------------------------------------------------------
    for i in range(1, param['WS']['N'] - 1):
        Ti[i] = float(odeint(EB_MIDDLE, Ti[i], np.linspace(0, 900))[-1])
    # -------------------------------------------------------------------------------------------------------------    
     
    
    # ---------------------------Return-Temperatures---------------------------------------------------------------
    return Ti


#%%

# Sek. mpkt bedingte Änderung
# Wärmeverluste an Umgebung
# Stofftransport
# Wärmeleitung zw. Layer

# WS_theta = [[10, 12, 14, 16, 18]]
# ws = WS(Vpkt_prim       = 200/3600/1000, 
#         T_prim_in       = 25, 
#         Vpkt_sek        = 100/3600/1000, 
#         T_sek_in        = 12, 
#         T_vektor_prev   = WS_theta[-1],
#         Tamb            = 18)


#%%
WS_theta = [[10]*param['WS']['N']]


dauer = 45
for i in range(0,dauer):
    print('Loop: ' + str(i) )
    ws = WS(Vpkt_prim       = 600/3600/1000, 
            T_prim_in       = 25, 
            Vpkt_sek        = 500/3600/1000, 
            T_sek_in        = 12, 
            T_vektor_prev   = WS_theta[-1],
            Tamb            = 300)
    WS_theta.append(ws)
 

for i in range(param['WS']['N']):
    plt.plot(range(0, dauer+1), [p[i] for p in WS_theta], label=[str(a) for a in list(range(param['WS']['N']))][i], color = plt.cm.RdYlBu(np.linspace(0, 1, param['WS']['N']))[::-1][i])
plt.legend()
plt.grid()
plt.xticks([])
plt.ylim(8,26)
plt.xlim(0, 2000*20/900)








#%%












