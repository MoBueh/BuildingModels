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
    'k'             : 0.002,    # W/m²K     Wärmetransferkoeffizient
    'lambda_eff'    : 0.15    # W/mK      vertikale Wärmeleitfähigkeit von Wasser     
    }




#%%

# mpkt_i = 1000/3600/1000*980 - 500/3600/1000*980 # when positive, then TES is being charged by the HP. When negative the TES is being discharged. Charge and discharged dont happen at the same time. It is just the result of mdot what determines what is happening.

def WS(Vpkt_prim, T_prim_in, Vpkt_sek, T_sek_in, T_vektor_prev, Tamb):
    #                        _______________
    # Vpkt_prim, T_prim_in  |               | T_sek_out
    # --------------------->|_______________|--------------------->
    #                       |               |
    #                       |_______________|
    #                       |               |
    #                       |_______________|
    #            T_prim_out |               | Vpkt_sek, T_sek_in
    # <---------------------|_______________|<---------------------   
    
    """"------------------------------Parameter------------------------------"""
    z_i = param['WS']['H']/param['WS']['N']                                     # (25) Höhe einer Schicht
    A_ext = np.pi*param['WS']['D']*z_i                                          # (26) Außenfläche einer Schicht
    A_i = np.pi*(param['WS']['D']-2*param['WS']['th'])**2/4                     # (27) Fläche zwischen zwei Layern
    m_i = A_i*z_i*prop['wasser']['rho']                                         # (28) Masse pro Schicht
    
    """"------------------------------Massenströme---------------------------"""
    mpkt_prim = Vpkt_prim*prop['wasser']['rho']
    mpkt_sek = Vpkt_sek*prop['wasser']['rho']
    mpkt_i = mpkt_prim-mpkt_sek    
            
    Ti = copy.deepcopy(T_vektor_prev)
    # Ti = T_vektor_prev

    """ calculate the bottom layer first--------------------------------------------------"""
    def EB_BOTTOM(y, t): # Function for the bottom layer
        Ti[i] = y
        
        if mpkt_i>0:                
            dTdt = (0
                    -mpkt_sek*prop['wasser']['cp']*(Ti[i] - T_sek_in)           # Sek. mpkt-bedingte Änderung
                    -(param['WS']['k'] * A_ext * (Ti[i] - Tamb))                # Wärmeverluste an Umgebung
                    # + (mpkt_i * prop['wasser']['cp'] * (Ti[i+1] - Ti[i]))       # Stofftransport
                    + (mpkt_prim * prop['wasser']['cp'] * (Ti[i+1] - Ti[i]))       # Stofftransport
                    +((A_i*param['WS']['lambda_eff']/z_i) * (Ti[i+1] - Ti[i]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt
            # print(((A_i*param['WS']['lambda_eff']/z_i) * (Ti[i+1] - Ti[i])))
            
                                               
        else: #The model is still not working properly. If a +X is not added to T[i] then the bottom temperature is stucked in one value and does not allow to discharge the tank. Leroy will check it.           
            dTdt = ((
              ((A_i*param['WS']['lambda_eff']/z_i)*(Ti[i] - Ti[i+1])) \
            + (mpkt_i * prop['wasser']['cp'] * (Ti[i] - T_sek_in)) \
            - (param['WS']['k']*A_ext*(Ti[i] - Tamb))) / (m_i*prop['wasser']['cp']))
                
        return dTdt

    i = 0 #Bottom layer index
    y0 = Ti[i]  # initial condition  
    Ti[i] = float(odeint(EB_BOTTOM, y0, np.linspace(0, 20))[-1]) # solve ODE

    """ calculate and iterate through the middle layers------------------------------------- """
    
    def EB_MIDDLE(y, t): # Function for the middle layer
        Ti[i] = y
        
        
        if mpkt_i>0:
            dTdt = (0
                    +0 # kein Anschluss in diesen Layern                        # Sek/prim mpkt-bedingte Änderung
                    -(param['WS']['k']*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung
                    +(mpkt_i*prop['wasser']['cp']*(Ti[i+1] - Ti[i]))            # Stofftransport
                    +((A_i*param['WS']['lambda_eff'])/z_i * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt
            # print(((A_i*param['WS']['lambda_eff'])/z_i * (Ti[i+1] - 2*Ti[i] + Ti[i-1])))
            
        else:
            dTdt = ((((A_i*param['WS']['lambda_eff'])/z_i) * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) \
                    + (mpkt_i*prop['wasser']['cp']*(Ti[i] - Ti[i-1])) \
                    - (param['WS']['k']*A_ext*(Ti[i] - Tamb))) / (m_i*prop['wasser']['cp'])
        
        return dTdt
    
    for i in range(1, param['WS']['N'] - 1): #All middle layers index
        # if  mpkt_i>0:
        #     y0 = Ti[i+1] # take layer above the current layer as initial condition          
        # else:
        #     y0 = Ti[i-1] # take layer under the current layer as initial condition 
        y0 =  Ti[i]        

        Ti[i] = float(odeint(EB_MIDDLE, y0, np.linspace(0, 20))[-1]) # solve ODE
    
    """ calculate the top layer --------------------------------------------------------------"""
    def EB_TOP(y, t): # Function for the top layer
        Ti[i] = y
        
        if mpkt_i>0:                
            dTdt = (0
                    +(mpkt_prim * prop['wasser']['cp'] * (T_prim_in - Ti[i]))   # Prim. mpkt-bedingte Änderung
                    -(param['WS']['k']*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung 
                    +0 # keine Schicht oberhalb vorhanden                       # Stofftransport
                    +((A_i*param['WS']['lambda_eff']/z_i)*(Ti[i-1] - Ti[i]))    # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt  
            # print(((A_i*param['WS']['lambda_eff']/z_i)*(Ti[i-1] - Ti[i])))
        else:       
            dTdt = (
              ((A_i*param['WS']['lambda_eff']/z_i)*(Ti[i] - Ti[i-1])) \
            + (mpkt_i * prop['wasser']['cp'] * (Ti[i] - Ti[i-1])) \
            - (param['WS']['k']*A_ext*(Ti[i] - Tamb))) / (m_i*prop['wasser']['cp'])
                            
        return dTdt
    
    i = param['WS']['N'] - 1 # Top layer index
    y0 = Ti[i] #initial condition 
    Ti[i] = float(odeint(EB_TOP, y0, np.linspace(0, 20))[-1]) # solve ODE
    
    # T_iteration = [T_vektor_prev[0]] + Ti[1:-1] + [T_vektor_prev[-1]]    
    # return Ti, T_iteration
    return Ti
    # return [round(x, 2) for x in Ti]  
    # return [round(x, 3) for x in Ti], [round(x, 3) for x in T_iteration]
    # return [f"{x:.3f}" for x in Ti], [f"{x:.3f}" for x in T_iteration]


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


dauer = 2000
for i in range(0,dauer):
    print('Loop: ' + str(i) )
    ws = WS(Vpkt_prim       = 600/3600/1000, 
            T_prim_in       = 25, 
            Vpkt_sek        = 500/3600/1000, 
            T_sek_in        = 12, 
            T_vektor_prev   = WS_theta[-1],
            Tamb            = 18)
    WS_theta.append(ws)
 

for i in range(param['WS']['N']):
    # print(i)
    plt.plot(range(0, dauer+1), [p[i] for p in WS_theta], label=[str(a) for a in list(range(param['WS']['N']))][i], color = plt.cm.RdYlBu(np.linspace(0, 1, param['WS']['N']))[::-1][i])
    # plt.plot(range(0, 11), [p[i] for p in T], label=['Bot', 'Top'][i])
plt.legend()
plt.grid()
plt.xlim(0,2000)
plt.xticks([])

# plt.grid()







#%%












