# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 16:03:45 2023

@author: 49157
"""

#%% Directory
import os
os.chdir('C:/Users/49157/Desktop/Modellvergleich')

#%% Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Vergleich.B1R1C import B1R1C
from Vergleich.B5R1C_simplified import B5R1C
from Modelle.RadProc import RadProc
plt.rc('axes', axisbelow=True)
plt.rcParams['figure.dpi'] = 800

#%% input
inp = pd.read_csv('Input.csv', index_col = 0)

#%% Parameter-Set
param = {}

# 1R1C
param['B1R1C'] = {
    'UA' : 179.2,
    'C_m' : 26000000.0,
    }

# 5R1C
param['B5R1C'] = {
        'A_f': 100,                                                          # m²    conditioned area
        'Coe_Am': 2.5,                                                         # -     coefficient for the determination of mass-related area according to table 12
        'Coe_Cm': 260000,                                                      # J/K*m² coefficiant for the determination of Internal heat storage capacity according to table 12
        'V': 750,                                                            # m²    Volume of building
    'transparent': {
        'A'          : [12, 12, 12, 12, 0],                 # m²    Area of the component
        'U'          : 0.8,                                                    # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'g_tot'      : 0.4
        },
    'opaque': {
        'A'          : [63, 63, 63, 63, 100],         # m²    Area of the component
        'U'          : 0.4,                                                     # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'R'          : 0.05,                                                    # Wärmeübergang Strahlung (Absorptionskoeffizient *Oberflächenwärme-Durchlasswiderstand)
        }    
    }
#%% Simulation
days = 1

# 1R1C
GB_theta_i = [6.5]
# 5R1C
GB_theta_m = [6.5]
GB_theta_s = [6.5]
GB_theta_air = [6.5]
GB_theta_operative = [6.5]

for time in range(0,24*days,1):
    i_sol = []
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 0))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 90))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 180))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 270))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 0, pos_az = 270))
    # 1R1C
    GB_1R1C = B1R1C(phi_int = 0,
                    theta_e = inp.loc[time,'t'],
                    phi_HC = 0,
                    theta_prev = GB_theta_i[-1],
                    i_sol = sum(i_sol),
                    param = param['B1R1C'])
    GB_theta_i.append(GB_1R1C)
    # 5R1C
    GB_5R1C = B5R1C(phi_int = 0,
                    theta_e = inp.loc[time,'t'],
                    V_pkt = 20,
                    theta_sup = inp.loc[time,'t'],
                    phi_HC = 0,
                    theta_m_prev = GB_theta_m[-1],
                    i_sol = i_sol,
                    shading = [1]*5,
                    param = param['B5R1C'])
    GB_theta_m.append(GB_5R1C[0])
    GB_theta_s.append(GB_5R1C[1])
    GB_theta_air.append(GB_5R1C[2])
    GB_theta_operative.append(GB_5R1C[3])

df_1R1C = pd.DataFrame({'air':GB_theta_i[1:]})
df_5R1C = pd.DataFrame({'air' : GB_theta_air[1:],'surface' : GB_theta_s[1:], 'mass' : GB_theta_m[1:], 'operative' : GB_theta_operative[1:]})
del GB_theta_m, GB_theta_s, GB_theta_air, GB_theta_operative, time, GB_1R1C, GB_5R1C, GB_theta_i, i_sol

#%%
plt.plot(df_1R1C.index, df_1R1C['air'], label = '1R1C')
plt.plot(df_5R1C.index, df_5R1C['air'], label = '5R1C')
plt.plot(inp.index[0:24*days], inp['t'][0:24*days])
plt.grid()
plt.legend()


#%%





# plt.plot(Y)




















